#' (univariate) Time-Delay Embedding forecasting using Gaussian Processes
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
#' @inheritParams block_lnlp
#' @inheritParams simplex
#' @inheritParams tde_gp
#' @param ... other parameters. see 'Details'
#' @return If stats_only, then a data.frame with components for the parameters 
#'   and forecast statistics:
#' \tabular{ll}{
#'   \code{E} \tab embedding dimension\cr
#'   \code{tau} \tab time lag\cr
#'   \code{tp} \tab prediction horizon\cr
#'   \code{phi} \tab length-scale parameter\cr
#'   \code{v_e} \tab noise-variance parameter\cr
#'   \code{eta} \tab signal-variance parameter\cr
#'   \code{fit_params} \tab whether params were fitted or not\cr
#'   \code{num_pred} \tab number of predictions\cr
#'   \code{rho} \tab correlation coefficient between observations and 
#'     predictions\cr
#'   \code{mae} \tab mean absolute error\cr
#'   \code{rmse} \tab root mean square error\cr
#'   \code{perc} \tab percent correct sign\cr
#'   \code{p_val} \tab p-value that rho is significantly greater than 0 using 
#'     Fisher's z-transformation\cr
#'   \code{model_output} \tab data.frame with columns for the time index, 
#'     observations, mean-value for predictions, and independent variance for 
#'     predictions (if \code{stats_only == FALSE} or 
#'     \code{save_covariance_matrix == TRUE})\cr
#'   \code{covariance_matrix} \tab the full covariance matrix for predictions 
#'     (if \code{save_covariance_matrix == TRUE})\cr
#' }
#' @examples 
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' tde_gp(ts, lib = c(1, 100), pred = c(101, 200), E = 5)
#' 
tde_gp <- function(time_series, lib = c(1, NROW(time_series)), pred = lib, 
                   E = 1:10, tau = 1, tp = 1, 
                   phi = 0, v_e = 0, eta = 0, fit_params = TRUE, 
                   stats_only = TRUE, save_covariance_matrix = FALSE, 
                   silent = FALSE, ...)
{
    # setup lib and pred ranges
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)

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

    params <- expand.grid(E = E, tau = tau)
    output <- do.call(rbind, lapply(1:NROW(params), function(i) {
        E <- params$E[i]
        tau <- params$tau[i]
        
        # make block
        block <- make_block(block = data.frame(ts = time_series),
                            t = time, max_lag = E, tau = tau,
                            lib = lib, restrict_to_lib = FALSE)

        # pass along args to block_gp
        out_df <- block_gp(block, lib, pred, tp = tp, 
                           phi = phi, v_e = v_e, eta = eta, 
                           fit_params = fit_params, 
                           columns = 1:E, target_column = 1, 
                           stats_only = stats_only, 
                           save_covariance_matrix = save_covariance_matrix, 
                           first_column_time = TRUE, 
                           silent = TRUE, ...)
        
    }))
    
    return(data.frame(E = params$E, 
                      tau = params$tau, 
                      output))
}
