#' @name simplex
#' @aliases s_map
#' 
#' @title Perform univariate forecasting

#' @details \code{\link{simplex}} is typically applied, and the embedding dimension 
#' varied, to find an optimal embedding dimension for the data. Thus, the 
#' default parameters are set so that passing a time series as the only 
#' argument will run over E = 1:10 (embedding dimension), using leave-one-out 
#' cross-validation over the whole time series, and returning just the forecast 
#' statistics.
#' 
#' \code{\link{s_map}} is typically applied, with fixed embedding dimension, and theta 
#' varied, to test for nonlinear dynamics in the data. Thus, the default 
#' parameters are set so that passing a time series as the only  argument will 
#' run over a default list of thetas (0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 
#' 0.03, 0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6, and 8), using E = 1, 
#' leave-one-out cross-validation over the whole time series, and returning 
#' just the forecast statistics.
#'
#' \code{norm = 2} (default) uses the "L2 norm", Euclidean distance:
#'   \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
#'     }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' \code{norm = 1} uses the "L1 norm", Manhattan distance:
#'   \deqn{distance(a,b) := \sum_i{|a_i - b_i|}
#'     }{distance(a, b) := \sum|a_i - b_i|}
#' Other values generalize the L1 and L2 norm to use the given argument as the 
#'   exponent, P, as:
#'   \deqn{distance(a,b) := \sum_i{(a_i - b_i)^P}^{1/P}
#'     }{distance(a, b) := (\sum(a_i - b_i)^P)^(1/P)}
#' 
#' @inheritParams block_lnlp
#' @param time_series either a vector to be used as the time series, or a 
#'   data.frame or matrix with at least 2 columns (in which case the first 
#'   column will be used as the time index, and the second column as the time 
#'   series)
#' @param E the embedding dimensions to use for time delay embedding
#' @param tau the lag to use for time delay embedding
#' @param tp the prediction horizon (how far ahead to forecast)

#' @rdname simplex
#' 
#' @description \code{\link{simplex}} uses time delay embedding on a single time 
#'   series to generate an attractor reconstruction, and then applies the 
#'   simplex projection algorithm to make forecasts.
#' 
#' @return For \code{\link{simplex}}, a data.frame with components for the 
#'   parameters and forecast statistics:
#' \tabular{ll}{
#'   \code{E} \tab embedding dimension\cr
#'   \code{tau} \tab time lag\cr
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
#' @examples 
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' simplex(ts, lib = c(1, 100), pred = c(101, 200))
#' 
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' simplex(ts, stats_only = FALSE)
#'  
simplex <- function(time_series, lib = c(1, NROW(time_series)), pred = lib, 
                    norm = 2, 
                    E = 1:10, tau = 1, tp = 1, num_neighbors = "e+1", 
                    stats_only = TRUE, exclusion_radius = NULL, epsilon = NULL, 
                    silent = FALSE)
{
    # make new model object
    model <- new(LNLP)
    
    # setup data
    dat <- setup_time_and_time_series(time_series)
    time <- dat$time
    time_series <- dat$time_series
    model$set_time(time)
    model$set_time_series(time_series)
    
    # setup norm and pred types
    model$set_norm(norm)
    model$set_pred_type(2)
    
    # setup lib and pred ranges
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)
    model$set_lib(lib)
    model$set_pred(pred)

    # handle remaining arguments and flags
    setup_model_flags(model, exclusion_radius, epsilon, silent)

    # setup other params in data.frame
    params <- expand.grid(tp, num_neighbors, tau, E)
    names(params) <- c("tp", "nn", "tau", "E")
    params <- params[, c("E", "tau", "tp", "nn")]
    e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if (any(e_plus_1_index, na.rm = TRUE))
        params$nn <- params$E + 1
    params$nn <- as.numeric(params$nn)
    
    # check params
    idx <- sapply(seq(NROW(params)), function(i) {
        check_params_against_lib(params$E[i], params$tau[i], params$tp[i], lib, 
                                 silent = silent)})
    if (!any(idx))
    {
        stop("No valid parameter combinations to run, stopping.")
    }
    params <- params[idx, ]
    
    # apply model prediction function to params
    output <- lapply(seq(NROW(params)), function(i) {
        model$set_params(params$E[i], params$tau[i], params$tp[i], params$nn[i])
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
    
    return(cbind(params, do.call(rbind, output), row.names = NULL))
}

#' @rdname simplex
#' 
#' @description \code{\link{s_map}} is similar to \code{\link{simplex}}, but uses the S-map 
#'   algorithm to make forecasts.
#' @return For \code{\link{s_map}}, the same as for \code{\link{simplex}}, but 
#'   with additional columns:
#' \tabular{ll}{
#'   \code{theta} \tab the nonlinear tuning parameter\cr
#'   \code{smap_coefficients} \tab data.frame with columns for the s-map 
#'   coefficients (if \code{save_smap_coefficients == TRUE})\cr
#'   \code{smap_coefficient_covariances} \tab list of covariance matrices for 
#'   the s-map coefficients (if \code{save_smap_coefficients == TRUE})\cr
#' }
#' @examples 
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' s_map(ts, E = 2)
#' 
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' s_map(ts, E = 2, theta = 1, save_smap_coefficients = TRUE)
#' 
s_map <- function(time_series, lib = c(1, NROW(time_series)), pred = lib, 
                  norm = 2, 
                  E = 1, tau = 1, tp = 1, num_neighbors = 0, 
                  theta = c(0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 
                            0.3, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6, 8), 
                  stats_only = TRUE, exclusion_radius = NULL, epsilon = NULL, 
                  silent = FALSE, save_smap_coefficients = FALSE)
{
    # check inputs?
    
    # make new model object
    model <- new(LNLP)
    
    # setup data
    dat <- setup_time_and_time_series(time_series)
    time <- dat$time
    time_series <- dat$time_series
    model$set_time(time)
    model$set_time_series(time_series)
    
    # setup norm and pred types
    model$set_norm(norm)
    model$set_pred_type(1)
    
    # setup lib and pred ranges
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)
    model$set_lib(lib)
    model$set_pred(pred)
        
    # handle remaining arguments and flags
    setup_model_flags(model, exclusion_radius, epsilon, silent)
    
    # handle smap coefficients flag
    if (save_smap_coefficients)
    {
        stats_only <- FALSE
        model$save_smap_coefficients()
    }
    # setup other params in data.frame
    params <- expand.grid(theta, tp, num_neighbors, tau, E)
    names(params) <- c("theta", "tp", "nn", "tau", "E")
    params <- params[, c("E", "tau", "tp", "nn", "theta")]
    e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if (any(e_plus_1_index, na.rm = TRUE))
        params$nn <- params$E + 1
    params$nn <- as.numeric(params$nn)
    
    # check params
    idx <- sapply(seq(NROW(params)), function(i) {
        check_params_against_lib(params$E[i], params$tau[i], params$tp[i], lib, 
                                 silent = silent)})
    if (!any(idx))
    {
        stop("No valid parameter combinations to run, stopping.")
    }
    params <- params[idx, ]
    
    # apply model prediction function to params
    output <- lapply(1:NROW(params), function(i) {
        model$set_params(params$E[i], params$tau[i], params$tp[i], params$nn[i])
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
                df$smap_coefficients <- I(list(model$get_smap_coefficients()))
                df$smap_coefficient_covariances <- I(list(model$get_smap_coefficient_covariances()))
            }
        }
        return(df)
    })
    
    return(cbind(params, do.call(rbind, output), row.names = NULL))
}
