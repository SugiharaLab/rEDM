#' Perform generalized forecasting using Gaussian processes
#'
#' \code{block_gp} uses multiple time series given as input to generate an 
#' attractor reconstruction, and then applies Gaussian process regression to 
#' approximate the dynamics and make forecasts. This method is the 
#' generalized version of \code{tde_gp}, which constructs the block from time 
#' lags of a time series to pass into this function.
#' 
#' The default parameters are set so that passing a vector as the only 
#' argument will use that vector to predict itself one time step ahead. If a 
#' matrix or data.frame is given as the only argument, the first column will be 
#' predicted (one time step ahead), using the remaining columns as the 
#' embedding. Rownames will be converted to numeric if possible to be used as 
#' the time index, otherwise 1:NROW will be used instead. The default lib and 
#' pred are for leave-one-out cross-validation over the whole time series, and 
#' returning just the forecast statistics.
#' 
#' [[DETAILS ABOUT IMPLEMENTATION]]
#' 
#' @param block either a vector to be used as the time series, or a 
#'   data.frame or matrix where each column is a time series
#' @param lib a 2-column matrix (or 2-element vector) where each row specifes the 
#'   first and last *rows* of the time series to use for attractor reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param tp the prediction horizon (how far ahead to forecast)
#' @param phi length-scale parameter. see 'Details'
#' @param v_e noise-variance parameter. see 'Details'
#' @param tau signal-variance parameter. see 'Details'
#' @param columns either a vector with the columns to use (indices or names), 
#'   or a list of such columns
#' @param target_column the index (or name) of the column to forecast
#' @param stats_only specify whether to output just the forecast statistics or 
#'   the raw predictions for each run
#' @param first_column_time indicates whether the first column of the given 
#'   block is a time column (and therefore excluded when indexing)
#' @param fit_params specify whether to use MLE to estimate params over the lib
#' @param silent prevents warning messages from being printed to the R console
#' @param save_covariance_matrix specifies whether to include the full 
#'   covariance matrix with the output (and forces the full output as if 
#'   stats_only were set to FALSE)
#' @param ... other parameters. see 'Details'
#' @return If stats_only, then a data.frame with components for the parameters 
#'   and forecast statistics:
#' \tabular{ll}{
#'   cols \tab embedding\cr
#'   tp \tab prediction horizon\cr
#'   phi \tab length-scale parameter\cr
#'   v_e \tab noise-variance parameter\cr
#'   tau \tab signal-variance parameter\cr
#'   num_pred \tab number of predictions\cr
#'   rho \tab correlation coefficient between observations and predictions\cr
#'   mae \tab mean absolute error\cr
#'   rmse \tab root mean square error\cr
#'   perc \tab percent correct sign\cr
#'   p_val \tab p-value that rho is significantly greater than 0 using Fisher's 
#'   z-transformation\cr
#' }
#' Otherwise, a list where the number of elements is equal to the number of runs 
#'   (unique parameter combinations). Each element is a list with the following 
#'   components:
#' \tabular{ll}{
#'   params \tab data.frame of parameters (embedding, tp, nn)\cr
#'   
#'   model_output \tab data.frame with columns for the time index, observations, 
#'     and predictions\cr
#'   stats \tab data.frame of forecast statistics\cr
#' }
#' @examples 
#' data("two_species_model")
#' block <- two_species_model[1:200,]
#' block_gp(block, columns = c("x", "y"), first_column_time = TRUE)
#' @export
block_lnlp <- function(block, lib = c(1, NROW(block)), pred = lib, 
                       tp = 1, phi = 0, v_e = 0, tau = 0, 
                       columns = NULL, target_column = 1, 
                       stats_only = TRUE, first_column_time = FALSE, 
                       fit_params = TRUE, silent = FALSE, 
                       save_covariance_matrix = FALSE, ...)
{
    convert_to_column_indices <- function(columns)
    {
        if(is.numeric(columns))
        {
            if(any(columns > NCOL(block)))
                warning("Some column indices exceed the number of columns and were ignored.")
            return(columns[columns <= NCOL(block)])
        }
        # else
        indices <- match(columns, col_names)
        if(any(is.na(indices)))
            warning("Some column names could not be matched and were ignored.")
        return(indices[is.finite(indices)])
    }
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

    # setup lib and pred ranges
    if(is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if(is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    
    if(!all(lib[,2] >= lib[,1]))
        warning("Some library rows look incorrectly formatted, please check the lib argument.")
    if(!all(pred[,2] >= pred[,1]))
        warning("Some library rows look incorrectly formatted, please check the pred argument.")
    
    # set indices for lib and pred
    model$set_lib(lib)
    model$set_pred(pred)
    
    
    # define x_lib
    model$set_block(data.matrix(block))

    # define y_lib
    model$set_target_column(convert_to_column_indices(target_column))

    best_params <- c(phi = phi, v_e = v_e, tau = tau)

    # check if fitting params
    if(fit_params)
    {
        best_params <- get_mle_params_gp(x_lib = x_lib, y_lib = y_lib, 
                                         params_init = best_params, 
                                         ...) 
    } 
    
    # compute mean and covariance for pred
    out <- compute_gp(x_lib = x_lib, y_lib = y_lib,
                      params = best_params, cov_matrix = save_covariance_matrix, 
                      x_pred = x_pred, ...)
    
    # massage gp output into correct return format
    if (stats_only)
    {
        stats <- lapply(1:NROW(params), function(i) {
            model$set_embedding(columns[[params$embedding[i]]])
            model$set_params(params$tp[i], params$nn[i])
            model$set_theta(params$theta[i])
            model$run()
            return(model$get_stats())
        })
        params$embedding <- sapply(params$embedding, function(i) {
            paste(columns[[i]], sep = "", collapse = ", ")})
        output <- cbind(params, do.call(rbind, stats), row.names = NULL)
    } else {
        output <- lapply(1:NROW(params), function(i) {
            model$set_embedding(columns[[params$embedding[i]]])
            model$set_params(params$tp[i], params$nn[i])
            model$set_theta(params$theta[i])
            model$run()
            
            if(short_output)
            {
                to_return <- list(params = params[i,],
                                  embedding = paste(columns[[params$embedding[i]]], sep = "", collapse = ", "), 
                                  model_output = model$get_short_output(), 
                                  stats = model$get_stats())
                if(save_smap_coefficients)
                    to_return$smap_coefficients <- do.call(rbind, model$get_short_smap_coefficients())
            } else {
                to_return <- list(params = params[i,],
                                  embedding = paste(columns[[params$embedding[i]]], sep = "", collapse = ", "), 
                                  model_output = model$get_output(), 
                                  stats = model$get_stats())
                if(save_smap_coefficients)
                    to_return$smap_coefficients <- do.call(rbind, model$get_smap_coefficients())
            }
            return(to_return)
        })
    }
    return(output)
}