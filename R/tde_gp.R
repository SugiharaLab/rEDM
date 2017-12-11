#' Perform univariate forecasting using Gaussian processes
#' 
#' \code{\link{tde_gp}} is used in the same vein as \code{\link{simplex}} or \code{\link{s_map}} to 
#' do time series forecasting using Gaussian processes. Here, the default 
#' parameters are set so that passing a time series as the only argument will 
#' run over E = 1:10 (embedding dimension) to created a lagged block, and 
#' passing in that block and all remaining arguments into \code{\link{block_gp}}.
#' 
#' See \code{\link{block_gp}} for implementation details of the Gaussian process 
#' regression.
#' 
#' @param time_series either a vector to be used as the time series, or a 
#'   data.frame or matrix with at least 2 columns (in which case the first 
#'   column will be used as the time index, and the second column as the time 
#'   series)
#' @param lib a 2-column matrix (or 2-element vector) where each row specifes 
#'   the first and last *rows* of the time series to use for attractor 
#'   reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param E the embedding dimensions to use for time delay embedding
#' @param tau the lag to use for time delay embedding
#' @param tp the prediction horizon (how far ahead to forecast)
#' @param phi length-scale parameter. see 'Details'
#' @param v_e noise-variance parameter. see 'Details'
#' @param eta signal-variance parameter. see 'Details'
#' @param fit_params specify whether to use MLE to estimate params over the lib
#' @param stats_only specify whether to output just the forecast statistics or 
#'   the raw predictions for each run
#' @param save_covariance_matrix specifies whether to include the full 
#'   covariance matrix with the output (and forces the full output as if 
#'   stats_only were set to FALSE)
#' @param silent prevents warning messages from being printed to the R console
#' @param ... other parameters. see 'Details'
#' @return If stats_only, then a data.frame with components for the parameters 
#'   and forecast statistics:
#' \tabular{ll}{
#'   E \tab embedding dimension\cr
#'   tp \tab prediction horizon\cr
#'   phi \tab length-scale parameter\cr
#'   v_e \tab noise-variance parameter\cr
#'   eta \tab signal-variance parameter\cr
#'   fit_params \tab whether params were fitted or not\cr
#'   num_pred \tab number of predictions\cr
#'   rho \tab correlation coefficient between observations and predictions\cr
#'   mae \tab mean absolute error\cr
#'   rmse \tab root mean square error\cr
#'   perc \tab percent correct sign\cr
#'   p_val \tab p-value that rho is significantly greater than 0 using Fisher's 
#'   z-transformation\cr
#' }
#' If stats_only is FALSE or save_covariance_matrix is TRUE, then there is an 
#' additional list-column variable:
#' \tabular{ll}{
#'   model_output \tab data.frame with columns for the time index, 
#'   observations, and mean-value for predictions\cr
#' }
#' If save_covariance_matrix is TRUE, then there is an additional list-column 
#'   variable:
#' \tabular{ll}{
#'   covariance_matrix \tab covariance matrix for predictions\cr
#' }
#' @examples 
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' tde_gp(ts, lib = c(1, 100), pred = c(101, 200), E = 5)
#' @export
tde_gp <- function(time_series, lib = c(1, NROW(time_series)), pred = lib, 
                   E = 1:10, tau = 1, tp = 1, 
                   phi = 0, v_e = 0, eta = 0, fit_params = TRUE, 
                   stats_only = TRUE, save_covariance_matrix = FALSE, 
                   silent = FALSE, ...)
{
    # restructure lib and pred if necessary
    if (is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if (is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    
    # setup data
    if (is.vector(time_series)) {
        if (!is.null(names(time_series))) {
            time <- as.numeric(names(time_series))
            if (any(is.na(time)))
                time <- seq_along(time_series)
        } else {
            time <- seq_along(time_series)
        }
    } else if ((is.matrix(time_series) || is.data.frame(time_series)) && 
               NCOL(time_series) >= 2) {
        time <- time_series[, 1]
        time_series <- time_series[, 2]
    }
    n <- length(time_series)
    
    params <- expand.grid(E = E, tau = tau)
    output <- do.call(rbind, lapply(1:NROW(params), function(i) {
        E <- params$E[i]
        tau <- params$tau[i]
        
        # make block
        block <- matrix(NA, nrow = n, ncol = E + 1)
        block[, 1] <- time
        block[, 2] <- time_series
        if (E > 1)
        {
            for (lag in 2:E)
            {
                # add lag of previous lag
                block[, lag + 1] <- c(rep.int(NA, tau), block[1:(n - tau), lag])
                
                # set NAs for beginning of lib sections
                block[as.vector(outer(lib[, 1], 1:tau, "+")) - 1, lag + 1] <- NA
            }
        }

        # pass along args to block_gp
        out_df <- block_gp(block, lib, pred, tp = tp, 
                           phi = phi, v_e = v_e, eta = eta, 
                           fit_params = fit_params, 
                           columns = 1 + (1:E), target_column = 2, 
                           stats_only = stats_only, 
                           save_covariance_matrix = save_covariance_matrix, 
                           silent = TRUE, ...)
        
    }))
    
    return(data.frame(E = params$E, 
                      tau = params$tau, 
                      output))
}
