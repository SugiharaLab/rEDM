#' Perform generalized forecasting using simplex projection or s-map
#'
#' \code{block_lnlp} uses multiple time series given as input to generate an 
#' attractor reconstruction, and then applies the simplex projection or s-map 
#' algorithm to make forecasts. This method generalizes the \code{simplex} and 
#' \code{s-map} routines, and allows for "mixed" embeddings, where multiple time 
#' series can be used as different dimensions of an attractor reconstruction.
#' 
#' The default parameters are set so that passing a vector as the only 
#' argument will use that vector to predict itself one time step ahead. If a 
#' matrix or data.frame is given as the only argument, the first column will be 
#' predicted, using the remaining columns as the embedding. Rownames will be 
#' converted to numeric if possible to be used as the time index, otherwise 
#' 1:NROW will be used instead. The default lib and pred are for leave-one-out 
#' cross-validation over the whole time series, and returning just the forecast 
#' statistics.
#' 
#' norm_type "L2 norm" (default) uses the typical Euclidean distance:
#' \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}}{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' norm_type "L1 norm" uses the Manhattan distance:
#' \deqn{distance(a,b) := \sum_i{|a_i - b_i|}}{distance(a, b) := \sum|a_i - b_i|}
#' 
#' method "simplex" (default) uses the simplex projection forecasting algorithm
#' 
#' method "s-map" uses the s-map forecasting algorithm
#' 
#' @param block either a vector to be used as the time series, or a 
#'   data.frame or matrix where each column is a time series
#' @param lib a 2-column matrix (or 2-element vector) where each row specifes the 
#'   first and last *rows* of the time series to use for attractor reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param norm_type the distance function to use. see 'Details'
#' @param method the prediction method to use. see 'Details'
#' @param tp the prediction horizon (how far ahead to forecast)
#' @param num_neighbors the number of nearest neighbors to use (any of "e+1", 
#'   "E+1", "e + 1", "E + 1" will peg this parameter to E+1 for each run, any
#'   value < 1 will use all possible neighbors.)
#' @param columns either a vector with the columns to use (indices or names), 
#'   or a list of such columns
#' @param target_column the index (or name) of the column to forecast
#' @param stats_only specify whether to output just the forecast statistics or 
#'   the raw predictions for each run
#' @param first_column_time indicates whether the first column of the given 
#'   block is a time column (and therefore excluded when indexing)
#' @param exclusion_radius excludes vectors from the search space of nearest 
#'   neighbors if their *time index* is within exclusion_radius (NULL turns 
#'   this option off)
#' @param epsilon excludes vectors from the search space of nearest neighbors 
#'   if their *distance* is farther away than epsilon (NULL turns this option 
#'   off)
#' @param theta the nonlinear tuning parameter (theta is only relevant if 
#'   method == "s-map")
#' @param silent prevents warning messages from being printed to the R console
#' @return If stats_only, then a data.frame with components for the parameters 
#'   and forecast statistics:
#' \tabular{ll}{
#'   cols \tab embedding\cr
#'   tp \tab prediction horizon\cr
#'   nn \tab number of neighbors\cr
#'   num_pred \tab number of predictions\cr
#'   rho \tab correlation coefficient between observations and predictions\cr
#'   mae \tab mean absolute error\cr
#'   rmse \tab root mean square error
#' }
#' Otherwise, a list where the number of elements is equal to the number of runs 
#'   (unique parameter combinations). Each element is a list with the following 
#'   components:
#' \tabular{ll}{
#'   params \tab data.frame of parameters (E, tau, tp, nn)\cr
#'   model_output \tab data.frame with columns for the time index, observations, 
#'     and predictions\cr
#'   stats \tab data.frame of forecast statistics (num_pred, rho, mae, rmse)\cr
#' }
#' @export 

# KNOWN BUGS: NaN flags are set based on first embedding, so order of embeddings can influence results

block_lnlp <- function(block, lib = c(1, NROW(block)), pred = c(1, NROW(block)), 
                       norm_type = c("L2 norm", "L1 norm"), 
                       method = c("simplex", "s-map"), 
                       tp = 1, num_neighbors = "e+1", columns = NULL, 
                       target_column = 1, stats_only = TRUE, first_column_time = FALSE, 
                       exclusion_radius = NULL, epsilon = NULL, theta = NULL, silent = FALSE)
{
    convert_to_column_indices <- function(columns)
    {
        if(is.numeric(columns))
        {
            if(any(columns > NCOL(block)))
                message("Warning: some column indices exceed the number of columns and were ignored.")
            return(columns[columns <= NCOL(block)])
        }
        # else
        indices <- match(columns, col_names)
        if(any(is.na(indices)))
            message("Warning: some column names could not be matched and were ignored.")
        return(indices[is.finite(indices)])
    }
    # make new model object
    model <- new(BlockLNLP)
    
    # setup data
    if(first_column_time)
    {
        if(is.vector(block))
            time <- block
        else
        {
            time <- block[,1]
            block <- block[,-1]
        }
    }
    else
    {
        time <- rownames(block)
    }
    if (is.null(time))
    {
        time <- 1:NROW(block)
    } else {
        time <- as.numeric(time)
        if(any(is.na(time)))
            time <- 1:NROW(block)
    }
    col_names <- colnames(block)
    model$set_time(time)
    model$set_block(data.matrix(block))
    model$set_target_column(convert_to_column_indices(target_column))
    
    # setup norm and pred types
    model$set_norm_type(switch(match.arg(norm_type), "L2 norm" = 2, "L1 norm" = 1))
    model$set_pred_type(switch(match.arg(method), "simplex" = 2, "s-map" = 1))
    if(match.arg(method) == "s-map")
    {
        if(is.null(theta))
            theta <- 0        
    }
    
    # setup lib and pred ranges
    if (is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if (is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    model$set_lib(lib)
    model$set_pred(pred)
    
    # handle exclusion radius
    if (is.null(exclusion_radius))
        exclusion_radius = -1;
    model$set_exclusion_radius(exclusion_radius)
    
    # handle epsilon
    if (is.null(epsilon))
        epsilon = -1;
    model$set_epsilon(epsilon)
    
    # handle silent flag
    if (silent)
        model$suppress_warnings()
    
    # convert embeddings to column indices
    if(is.null(col_names)) {
        col_names <- paste("ts_", seq_len(NCOL(block)))
    }
    if(is.null(columns)) {
        columns <- list(1:NCOL(block))
    } else if(is.list(columns)) {
        columns <- lapply(columns, function(embedding) {
            convert_to_column_indices(embedding)
        })
    } else if(is.vector(columns)) {
        columns <- list(convert_to_column_indices(columns))
    } else if(is.matrix(columns)) {
        columns <- lapply(1:NROW(columns), function(i) convert_to_column_indices(columns[i,]))
    }
    embedding_index <- seq_along(columns)
    
    # setup other params in data.frame
    if(match.arg(method) == "s-map")
    {
        params <- expand.grid(tp, num_neighbors, theta, embedding_index)
        names(params) <- c("tp", "nn", "theta", "embedding")
        params <- params[,c("embedding", "tp", "nn", "theta")]
        e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
        if (any(e_plus_1_index, na.rm = TRUE))
            params$nn <- 1 + sapply(columns, length)
        # apply model prediction function to params
        if (stats_only)
        {
            stats <- lapply(1:NROW(params), function(i) {
                model$set_params(params$tp[i], params$nn[i])
                model$set_theta(params$theta[i])
                model$set_embedding(columns[[params$embedding[i]]])
                model$run()
                return(model$get_stats())
            })
            params$embedding <- sapply(params$embedding, function(i) {
                paste(columns[[i]], sep = "", collapse = ", ")})
            output <- cbind(params, do.call(rbind, stats))
        } else {
            output <- lapply(1:NROW(params), function(i) {
                model$set_params(params$tp[i], params$nn[i])
                model$set_theta(params$theta[i])
                model$set_embedding(columns[[params$embedding[i]]])
                model$run()
                
                return(list(params = params[i,], 
                            embedding = paste(columns[[params$embedding[i]]], sep = "", collapse = ", "), 
                            model_output = model$get_output(), 
                            stats = model$get_stats()))
            })
        }
    }
    else
    {
        params <- expand.grid(tp, num_neighbors, embedding_index)
        names(params) <- c("tp", "nn", "embedding")
        params <- params[,c("embedding", "tp", "nn")]
        e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
        if (any(e_plus_1_index, na.rm = TRUE))
            params$nn <- 1 + sapply(columns, length)
        # apply model prediction function to params
        if (stats_only)
        {
            stats <- lapply(1:NROW(params), function(i) {
                model$set_params(params$tp[i], params$nn[i])
                model$set_embedding(columns[[params$embedding[i]]])
                model$run()
                return(model$get_stats())
            })
            params$embedding <- sapply(params$embedding, function(i) {
                paste(columns[[i]], sep = "", collapse = ", ")})
            output <- cbind(params, do.call(rbind, stats))
        } else {
            output <- lapply(1:NROW(params), function(i) {
                model$set_params(params$tp[i], params$nn[i])
                model$set_embedding(columns[[params$embedding[i]]])
                model$run()
                
                return(list(params = params[i,], 
                            embedding = paste(columns[[params$embedding[i]]], sep = "", collapse = ", "), 
                            model_output = model$get_output(), 
                            stats = model$get_stats()))
            })
        }
    }
    return(output)
}