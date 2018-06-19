#' Perform generalized forecasting using simplex projection or s-map
#'
#' \code{\link{block_lnlp}} uses multiple time series given as input to generate 
#'   an attractor reconstruction, and then applies the simplex projection or 
#'   s-map algorithm to make forecasts. This method generalizes the 
#'   \code{\link{simplex}} and \code{\link{s_map}} routines, and allows for 
#'   "mixed" embeddings, where multiple time series can be used as different 
#'   dimensions of an attractor reconstruction.
#' 
#' The default parameters are set so that passing a vector as the only argument
#'   will use that vector to predict itself one time step ahead. If a matrix or 
#'   data.frame is given as the only argument, the first column will be 
#'   predicted (one time step ahead), using the remaining columns as the 
#'   embedding. Rownames will be converted to numeric if possible to be used as 
#'   the time index, otherwise 1:NROW will be used instead. The default lib and 
#'   pred are for leave-one-out cross-validation over the whole time series, 
#'   and returning just the forecast statistics.
#' 
#' norm_type "L2 norm" (default) uses the typical Euclidean distance:
#' \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
#' }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' norm_type "L1 norm" uses the Manhattan distance:
#' \deqn{distance(a,b) := \sum_i{|a_i - b_i|}
#' }{distance(a, b) := \sum|a_i - b_i|}
#' norm type "P norm" uses the P norm, generalizing the L1 and L2 norm to use 
#'   $P$ as the exponent:
#' \deqn{distance(a,b) := \sum_i{(a_i - b_i)^P}^{1/P}
#' }{distance(a, b) := (\sum(a_i - b_i)^P)^(1/P)}
#' 
#' method "simplex" (default) uses the simplex projection forecasting algorithm
#' 
#' method "s-map" uses the s-map forecasting algorithm
#' 
#' @param block either a vector to be used as the time series, or a 
#'   data.frame or matrix where each column is a time series
#' @param lib a 2-column matrix (or 2-element vector) where each row specifies 
#'   the first and last *rows* of the time series to use for attractor 
#'   reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param norm_type the distance function to use. see 'Details'
#' @param P the exponent for the P norm
#' @param method the prediction method to use. see 'Details'
#' @param tp the prediction horizon (how far ahead to forecast)
#' @param num_neighbors the number of nearest neighbors to use. Note that the 
#'   default value will change depending on the method selected. (any of "e+1", 
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
#' @param save_smap_coefficients specifies whether to include the s_map 
#'   coefficients with the output (and forces stats_only = FALSE, as well)
#' @return A data.frame with components for the parameters and forecast 
#'   statistics:
#' \tabular{ll}{
#'   \code{cols} \tab embedding\cr
#'   \code{tp} \tab prediction horizon\cr
#'   \code{nn} \tab number of neighbors\cr
#'   \code{num_pred} \tab number of predictions\cr
#'   \code{rho} \tab correlation coefficient between observations and 
#'     predictions\cr
#'   \code{mae} \tab mean absolute error\cr
#'   \code{rmse} \tab root mean square error\cr
#'   \code{perc} \tab percent correct sign\cr
#'   \code{p_val} \tab p-value that rho is significantly greater than 0 using 
#'     Fisher's z-transformation\cr
#'   \code{const_rho} \tab same as \code{rho}, but for the constant predictor\cr
#'   \code{const_mae} \tab same as \code{mae}, but for the constant predictor\cr
#'   \code{const_rmse} \tab same as \code{rmse}, but for the constant predictor\cr
#'   \code{const_perc} \tab same as \code{perc}, but for the constant predictor\cr
#'   \code{const_p_val} \tab same as \code{p_val}, but for the constant predictor\cr
#'   \code{model_output} \tab data.frame with columns for the time index, 
#'     observations, predictions, and estimated prediction variance
#'     (if \code{stats_only == FALSE})\cr
#' }
#' If "s-map" is the method, then the same, but with additional columns:
#' \tabular{ll}{
#'   \code{theta} \tab the nonlinear tuning parameter\cr
#'   \code{smap_coefficients} \tab data.frame with columns for the s-map 
#'   coefficients (if \code{save_smap_coefficients == TRUE})\cr
#'   \code{smap_coefficient_covariances} \tab list of covariance matrices for 
#'   the s-map coefficients (if \code{save_smap_coefficients == TRUE})\cr
#' }
#' @examples 
#' data("two_species_model")
#' block <- two_species_model[1:200,]
#' block_lnlp(block, columns = c("x", "y"), first_column_time = TRUE)
#' 
block_lnlp <- function(block, lib = c(1, NROW(block)), pred = lib, 
                       norm_type = c("L2 norm", "L1 norm", "P norm"), P = 0.5, 
                       method = c("simplex", "s-map"), 
                       tp = 1, 
                       num_neighbors = switch(match.arg(method), 
                                              "simplex" = "e+1", "s-map" = 0), 
                       columns = NULL, 
                       target_column = 1, stats_only = TRUE, 
                       first_column_time = FALSE, 
                       exclusion_radius = NULL, epsilon = NULL, theta = NULL, 
                       silent = FALSE, save_smap_coefficients = FALSE)
{
    # make new model object
    model <- new(BlockLNLP)
    
    # setup data
    block <- setup_time_and_data_block(model, first_column_time, block)
    model$set_target_column(convert_to_column_indices(target_column, block, 
                                                      silent = silent))
    
    # setup norm and pred types
    model$set_norm_type(switch(match.arg(norm_type), 
                               "P norm" = 3, "L2 norm" = 2, "L1 norm" = 1))
    model$set_p(P)
    model$set_pred_type(switch(match.arg(method), "simplex" = 2, "s-map" = 1))
    if (match.arg(method) == "s-map")
    {
        if (is.null(theta))
            theta <- 0
        if (save_smap_coefficients)
        {
            stats_only <- FALSE
        }
        model$save_smap_coefficients()
    }
    
    # setup lib and pred ranges
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)
    model$set_lib(lib)
    model$set_pred(pred)
    
    # handle remaining arguments and flags
    setup_model_flags(model, exclusion_radius, epsilon, silent)
    
    # convert embeddings to column indices
    if (is.null(names(block)))
    {
        names(block) <- paste("ts_", seq_len(NCOL(block)))
    }
    if (is.null(columns))
    {
        columns <- list(1:NCOL(block))
    } else if (is.list(columns)) {
        columns <- lapply(columns, function(embedding) {
            convert_to_column_indices(embedding, block, silent = silent)
        })
    } else if (is.vector(columns)) {
        columns <- list(convert_to_column_indices(columns, block, silent = silent))
    } else if (is.matrix(columns)) {
        columns <- lapply(1:NROW(columns), function(i) {
            convert_to_column_indices(columns[i,], block, silent = silent)})
    }
    embedding_index <- seq_along(columns)
    
    # setup other params in data.frame
    if (match.arg(method) == "s-map")
    {
        params <- expand.grid(tp, num_neighbors, theta, embedding_index)
        names(params) <- c("tp", "nn", "theta", "embedding")
        params <- params[,c("embedding", "tp", "nn", "theta")]
        e_plus_1_index <- match(num_neighbors, 
                                c("e+1", "E+1", "e + 1", "E + 1"))
        if (any(e_plus_1_index, na.rm = TRUE))
            params$nn <- 1 + sapply(columns, length)
        params$nn <- as.numeric(params$nn)
        
        # check params
        idx <- sapply(seq(NROW(params)), function(i) {
            check_params_against_lib(1, 0, params$tp[i], lib, silent = silent)})
        if (!any(idx))
        {
            stop("No valid parameter combinations to run, stopping.")
        }
        params <- params[idx, ]
        
        # apply model prediction function to params
        output <- lapply(1:NROW(params), function(i) {
            model$set_embedding(columns[[params$embedding[i]]])
            model$set_params(params$tp[i], params$nn[i])
            model$set_theta(params$theta[i])
            model$run()
            if (stats_only)
            {
                df <- model$get_stats()
            } else {
                df <- model$get_stats()
                df$model_output <- I(list(model$get_output()))
                if (save_smap_coefficients)
                {
                    df$smap_coefficients <- 
                        I(list(model$get_smap_coefficients()))
                    df$smap_coefficient_covariances <- 
                        I(list(model$get_smap_coefficient_covariances()))
                }
            }
            return(df)
        })
    } else {
        # simplex
        params <- expand.grid(tp, num_neighbors, embedding_index)
        names(params) <- c("tp", "nn", "embedding")
        params <- params[, c("embedding", "tp", "nn")]
        e_plus_1_index <- match(num_neighbors, 
                                c("e+1", "E+1", "e + 1", "E + 1"))
        if (any(e_plus_1_index, na.rm = TRUE))
            params$nn <- 1 + sapply(columns, length)
        params$nn <- as.numeric(params$nn)
        
        # check params
        idx <- sapply(seq(NROW(params)), function(i) {
            check_params_against_lib(1, 0, params$tp[i], lib, silent = silent)})
        if (!any(idx))
        {
            stop("No valid parameter combinations to run, stopping.")
        }
        params <- params[idx, ]
        
        # apply model prediction function to params
        output <- lapply(1:NROW(params), function(i) {
            model$set_embedding(columns[[params$embedding[i]]])
            model$set_params(params$tp[i], params$nn[i])
            model$run()
            if (stats_only)
            {
                df <- model$get_stats()
            } else {
                df <- model$get_stats()
                df$model_output <- I(list(model$get_output()))
            }
            return(df)
        })
    }
    
    # create embedding column in params
    params$embedding <- sapply(params$embedding, function(i) {
        paste(columns[[i]], sep = "", collapse = ", ")})
    return(cbind(params, do.call(rbind, output), row.names = NULL))
}