#' (generalized) Block forecasting using Gaussian Processes
#'
#' \code{\link{block_gp}} uses multiple time series given as input to generate an 
#' attractor reconstruction, and then applies Gaussian process regression to 
#' approximate the dynamics and make forecasts. This method is the 
#' generalized version of \code{\link{tde_gp}}, which constructs the block from 
#' lags of a time series to pass into this function.
#' 
#' The default parameters are set so that passing a vector as the only 
#' argument will use that vector to predict itself one time step ahead. If a 
#' matrix or data.frame is given as the only argument, the first column will be 
#' predicted (one time step ahead), using the remaining columns as the 
#' embedding. Rownames will be converted to numeric if possible to be used as 
#' the time index, otherwise 1:NROW will be used instead. The default lib and 
#' pred are to perform maximum likelihood estimation of the phi, v_e, and eta 
#' parameters over the whole time series, and return just the forecast 
#' statistics.
#' 
#' If phi, v_e, and eta parameters are given, all combinations of their values 
#' will be tried. If fit_params is also set to TRUE, these values will be the 
#' initial values for subsequent optimization of likelihood.
#' 
#' The basic model is:
#'   \deqn{y = f(x) + \textnormal{noise}
#'   }{y = f(x) + noise}
#' in which the function f(x) is modeled using a Gaussian process prior:
#'   \deqn{f \sim \textnormal{GP}(0, C)
#' }{f ~ GP(0, C)}
#' with mean = 0,
#' and covariance function, C, which is given by the squared-exponential kernel:
#'   \deqn{C_{ij} = eta * \exp(-phi^2 * ||x_i - x_j||^2)
#' }{C_{ij} = eta * exp(-phi^2 * ||x_i - x_j||^2)}
#' 
#' y is a realization from process f with normally-distributed i.i.d. process 
#' noise,
#'   \deqn{noise \sim \mathcal(N)(0, v_e)
#'   }{noise ~ N(0, v_e)}
#' such that the covariance of observations y_i and y_j is
#'   \deqn{K_{ij} = C_{ij} + v_e * \delta_{ij}}
#' where \eqn{\delta_{ij}} is the kronecker delta (i.e. it is 1 if \eqn{i = j} 
#' and 0 otherwise)
#'
#' From the model definition, the variance in y, after marginalizing over f, is 
#' given by eta + v_e. Thus to simplify specification of priors for the 
#' hyperparameters eta and v_e, the outputs y are normalized to zero mean and 
#' unit variance prior to fitting. This allows us to set (0, 1) bounds on 
#' eta and v_e which facilitates parameter estimation.  We set Beta(2, 2) 
#' priors for both eta and v_e to partition prior uncertainty equally across 
#' structural and process uncertainty.
#'
#' For a scalar input, the length-scale parameter phi controls the expected 
#' number of zero crossings on the unit interval as 
#'   \deqn{E(crossings) = \frac{\sqrt(2)}{\pi} \phi \approx 0.45 \phi
#'   }{E(crossings) = \sqrt(2)/\pi \phi \approx 0.45 \phi}
#' 
#' Thus to facilitate interpretation and prior specification, the distances in 
#' C are scaled by the max distance so that a model with \eqn{\phi = 2} would 
#' have roughly one zero crossing over the range of the data.  We assign phi a 
#' half-Normal prior with variance \eqn{\pi / 2} so that the prior mean phi is 
#' 1, which tends to avoid overfitting. 

#' To fit the GP we estimate eta, v_e, and phi by maximizing the posterior after
#' marginalizing over f(x). This is given by the multivariate normal likelihood 
#'   \deqn{logL = -1/2 \log{|K_d|}^{-1/2} y_d^T [K_d]^{-1} y_d
#' }{logL = -1/2 log(|K_d|)^(-1/2) y_d^T [K_d]^(-1) y_d}
#' where \eqn{K_d} is the matrix obtained by evaluating the covariance function 
#' at all pairs of inputs and \eqn{y_d} is the column vector of outputs. 
#' Predictions for new values of x are obtained by setting eta, v_e, and phi to 
#' the Maximum a Posteriori (MAP) estimates and using the GP conditional on the 
#' observed data.  Specifically, given \eqn{x_d} and \eqn{y_d}, 
#' the mean and variance for y evaluated at a new value of x are
#'   \deqn{E(y) = C(x, x_d) [K_d]^{-1} y_d
#'   }{E(y) = C(x, x_d) [K_d]^(-1) y_d}
#'   \deqn{V(y) = eta + v_e - C(x, x_d)[K_d]^{-1} C(x_d, x)
#'   }{V(y) = eta + v_e - C(x, x_d)[K_d]^(-1) C(x_d, x)}
#' where the vector \eqn{C(x, x_d)} is obtained by evaluating C at x and each 
#' of the observed inputs while holding eta, phi, and v_e at the MAP estimates.
#' @inheritParams block_lnlp
#' @param phi length-scale parameter. see 'Details'
#' @param v_e noise-variance parameter. see 'Details'
#' @param eta signal-variance parameter. see 'Details'
#' @param fit_params specify whether to use MLE to estimate params over the lib
#' @param save_covariance_matrix specifies whether to include the full 
#'   covariance matrix with the output (and forces the full output as if 
#'   stats_only were set to FALSE)
#' @param ... other parameters. see 'Details'
#' @return If stats_only, then a data.frame with components for the parameters 
#'   and forecast statistics:
#' \tabular{ll}{
#'   \code{embedding} \tab embedding\cr
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
#' block <- two_species_model[1:200,]
#' block_gp(block, columns = c("x", "y"), first_column_time = TRUE)
#' 
block_gp <- function(block, lib = c(1, NROW(block)), pred = lib, 
                     tp = 1, phi = 0, v_e = 0, eta = 0, 
                     fit_params = TRUE, 
                     columns = NULL, target_column = 1, 
                     stats_only = TRUE, save_covariance_matrix = FALSE, 
                     first_column_time = FALSE, silent = FALSE, ...)
{
    # setup data
    if (first_column_time)
    {
        if (is.vector(block))
            time <- block
        else
        {
            time <- block[, 1]
            block <- block[, -1, drop = FALSE]
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
        if (any(is.na(time)))
            time <- 1:NROW(block)
    }
    
    # setup embeddings
    col_names <- colnames(block)
    if (is.null(col_names)) {
        col_names <- paste("ts_", seq_len(NCOL(block)))
    }
    if (is.null(columns)) {
        columns <- list(1:NCOL(block))
    } else if (is.list(columns)) {
        columns <- lapply(columns, function(embedding) {
            convert_to_column_indices(embedding, block)
        })
    } else if (is.vector(columns)) {
        columns <- list(convert_to_column_indices(columns, block))
    } else if (is.matrix(columns)) {
        columns <- lapply(1:NROW(columns), function(i) {
            convert_to_column_indices(columns[i, ], block)})
    }
    
    # setup target    
    target_column <- convert_to_column_indices(target_column, block)
    
    # setup lib and pred ranges
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)

    params <- expand.grid(tp = tp, 
                          phi = phi, 
                          v_e = v_e, 
                          eta = eta, 
                          embedding_index = seq_along(columns))
    
    output <- do.call(rbind, lapply(1:NROW(params), function(i) {
        tp <- params$tp[i]
        phi <- params$phi[i]
        v_e <- params$v_e[i]
        eta <- params$eta[i]
        embedding <- columns[[params$embedding_index[i]]]
        
        # correct lib and pred for tp
        if (tp > 0)
        {
            lib[, 2] <- pmax(1, lib[, 2] - tp)
            pred[, 2] <- pmax(1, pred[, 2] - tp)
        } else if (tp < 0) {
            lib[, 1] <- pmin(NROW(block), lib[, 1] - tp)
            pred[, 1] <- pmax(NROW(block), pred[, 1] - tp)
        }
        
        # correct lib and pred if tp makes time series segments unusable
        lib <- lib[lib[, 2] >= lib[, 1], , drop = FALSE]
        pred <- pred[pred[, 2] >= pred[, 1], , drop = FALSE]
        if (NROW(lib) == 0)
            stop("No valid time series segments in lib ", 
                 "(after correcting for tp)")
        if (NROW(pred) == 0)
            stop("No valid time series segments in pred ", 
                 "(after correcting for tp)")
        
        # set indices for lib and pred
        lib_idx <- sort(unique(do.call(c, lapply(1:NROW(lib), function(i) {
            seq(from = lib[i, 1], to = lib[i, 2])}))))
        pred_idx <- sort(unique(do.call(c, lapply(1:NROW(pred), function(i) {
            seq(from = pred[i, 1], to = pred[i, 2])}))))
        
        # define inputs to fitting of GP (data, and params)
        x_lib <- as.matrix(block[lib_idx, embedding, drop = FALSE])
        y_lib <- block[lib_idx + tp, target_column]
        x_pred <- as.matrix(block[pred_idx, embedding, drop = FALSE])
        y_pred <- block[pred_idx + tp, target_column]
        time_pred <- time[pred_idx + tp]
        
        # filter x_lib and y_lib to finite values
        valid_lib_idx <- apply(is.finite(x_lib), 1, all) & is.finite(y_lib)
        if (sum(valid_lib_idx) < length(valid_lib_idx))
        {
            rEDM_warning("Trimmed ", length(valid_lib_idx), " lib points down to ", 
                         sum(valid_lib_idx), " valid ones.", silent = silent)
            x_lib <- x_lib[valid_lib_idx, , drop = FALSE]
            y_lib <- y_lib[valid_lib_idx]
        }
        if (NROW(x_lib) < 2)
            stop("Not enough data points in lib: n = ", NROW(x_lib))
        
        # filter x_pred and y_pred to finite values
        valid_pred_idx <- apply(is.finite(x_pred), 1, all)
        if (sum(valid_pred_idx) < length(valid_pred_idx))
        {
            rEDM_warning("Trimmed ", length(valid_pred_idx), " pred points down to ",
                    sum(valid_pred_idx), " valid ones.", silent = silent)
            x_pred <- x_pred[valid_pred_idx, , drop = FALSE]
            y_pred <- y_pred[valid_pred_idx]
        }


        # fit params if option is set, otherwise use as given
        best_params <- c(phi = phi, v_e = v_e, eta = eta)
        if (fit_params)
        {
            best_params <- get_mle_params_gp(x_lib = x_lib, y_lib = y_lib, 
                                             params_init = best_params, 
                                             ...) 
        } 
        
        # compute mean and covariance for pred
        out_gp <- compute_gp(x_lib = x_lib, y_lib = y_lib,
                             params = best_params, 
                             cov_matrix = save_covariance_matrix, 
                             x_pred = x_pred, ...)
        
        # compute stats for mean predictions
        stats <- compute_stats(y_pred, out_gp$mean_pred)
        
        # prepare output (default is param settings and stats)
        out_df <- data.frame(embedding = paste(embedding, sep = "", 
                                               collapse = ", "), 
                             tp = tp, 
                             phi = best_params["phi"], 
                             v_e = best_params["v_e"], 
                             eta = best_params["eta"], 
                             fit_params = fit_params, 
                             stats)
        
        # add in full output if requested
        if (!stats_only || save_covariance_matrix)
        {
            out_df$model_output <- list(data.frame(time = time_pred[valid_pred_idx], 
                                                   obs = y_pred, 
                                                   pred = out_gp$mean_pred, 
                                                   pred_var = out_gp$pred_var))
            if (save_covariance_matrix)
            {
                out_df$covariance_matrix <- list(out_gp$covariance_pred)
            }
        }
        row.names(out_df) <- NULL
        return(out_df)
    }))
    
    return(output)
}

get_mle_params_gp <- function(x_lib, y_lib, 
                              params_init = c(phi = 0, v_e = 0, eta = 0), 
                              mean_y = 0, var_y_lib = var(y_lib), 
                              param_rescaling_tol = 1e-3)
{
    likelihood_func <- function(x)
    {
        out <- compute_gp(x_lib, y_lib, params = x, 
                          mean_y = mean_y, var_y_lib = var_y_lib, 
                          gradient = TRUE, 
                          param_rescaling_tol = param_rescaling_tol)
        return(list(val = out$neg_log_likelihood, 
                    grad = out$gradient_neg_log_likelihood))
    }
    
    return(optim_rprop(likelihood_func, 
                       params_init))
}

optim_rprop <- function(func, x_init, 
                        max_iter = 200, delta_init = 0.1,
                        delta_min = 1e-6, delta_max = 50, 
                        eta_minus = 0.5, eta_plus = 1.2)
{
    ### Inputs:
    # func        function to minimize
    # x_init      initial value of x
    # max_iter    maximum # of iterations
    # delta_init  initial value for step size
    # delta_min   minimum value for step size
    # delta_max   maximum value for step size
    # eta_minus   scaling for decreasing step size
    # eta_plus    scaling for increasing step size
    
    ### Outputs:
    # x           "optimized" value of x (that minimizes func)
    
    ### Description
    # Perform function minimization using resilient backpropagation gradient 
    # descent
    
    # initial calulation of likelihood
    x <- x_init
    out <- func(x)
    s <- sqrt(sum(out$grad * out$grad))
    
    iter_count <- 0
    delta <- rep(delta_init, length(x_init))
    delta_f <- 10
    
    # begin iterations
    while ((s > 1e-4) &&                # STOP  if gradient is 0
           (iter_count < max_iter) &&   #    or if max iterations reached
           (delta_f > 1e-7))            #    or if no change in func output
    {
        # step 1: move
        x_new <- x - sign(out$grad) * delta
        out_new <- func(x_new)
        s <- sqrt(sum(out_new$grad * out_new$grad))
        delta_f <- abs(out_new$val / out$val - 1)
        
        # step 2: update step size
        grad_c <- sign(out$grad) * sign(out_new$grad)
        delta <- delta * (1 + 
                              (eta_plus - 1) * (grad_c > 0) + 
                              (eta_minus - 1) * (grad_c < 0))
        
        delta <- pmin(delta_max, pmax(delta_min, delta))
        
        # step 3: reset
        x <- x_new
        out <- out_new
        iter_count <- iter_count + 1
    }
    
    return(x)
}

compute_gp <- function(x_lib, y_lib, 
                       params = c(phi = 0, v_e = 0, eta = 0), 
                       mean_y = 0, 
                       max_x_lib = max(abs(x_lib)), 
                       var_y_lib = var(y_lib), 
                       x_pred = NULL, gradient = FALSE, cov_matrix = TRUE, 
                       param_rescaling_tol = 1e-3)
{
    ### Inputs: 
    # params     vector of parameters: [phi, v_e, eta]
    #              corresponding to length scale, process noise, pointwise 
    #              variance
    # x_lib      the n x E matrix of input states
    # y_lib      the n x 1 vector of corresponding outputs
    # mean_y     mean value of y (default is 0)
    # var_y_lib  variance for y (used because v_E and eta are proportional;
    #              default is variance of y_lib)
    # x_pred     the m x E matrix of input states to predict from (default is 
    #              NULL)
    # gradient   whether to compute the gradient of likelihood for params, used 
    #              for optimizing param values (default is FALSE)
    # param_rescaling_tol   tolerance for parameter rescaling (default is 0.001)
    
    ### Output
    # out        list with the following components:
    #   neg_log_likelihood   negative log likelihood (of lib vals and params)
    #
    #   (if gradient == TRUE)
    #     gradient_neg_log_likelihood   gradient of the negative log likelihood
    #                                     w.r.t. params
    #
    #   (if x_pred != NULL)
    #     mean_pred         mean value of predictions
    #     covariance_pred   covariance of predictions
    out <- list()
    
    ### Description
    # The basic model is:
    #     y = f(x) + noise
    # which we approximate using Gaussian Processes:
    #     f ~ GP(0, C)
    # with mean = 0,
    # and covariance C which is given by the squared-exponential kernel,
    #     C_ij = eta * exp(-phi^2 * ||x_i - x_j||^2)
    # y is a realization from process f with normally-distributed i.i.d. 
    #   process noise,
    #     e ~ N(0, v_e)
    # such that the covariance of observations y_i and y_j is
    #     K_ij = C_ij + v_e I_ij
    # with I the identity matrix or a kronecker delta
    
    ### Usage
    # (1)
    # create a function that computes likelihood | params, x_lib, y_lib
    #   this can be sent to a function optimizer to select best params value, 
    #   or sample from the posterior if we want to create distributions for the 
    #   params optionally, a gradient can be computed 
    #
    # (2)
    # compute predictions, y_pred | params, x_lib, y_lib, x_pred
    # note that the return value gives mean and variance
    
    ### Transform parameters
    # We do a transformation so that the parameters are constrained to (0, 1)
    v_e_min <- param_rescaling_tol
    v_e_max <- 1 - param_rescaling_tol
    eta_min <- param_rescaling_tol
    eta_max <- 1 - param_rescaling_tol
    
    phi <- exp(params["phi"]) / max_x_lib
    v_e <- v_e_min + (v_e_max - v_e_min) / (1 + exp(-params["v_e"]))
    eta <- eta_min + (eta_max - eta_min) / (1 + exp(-params["eta"]))
    d_params <- c(phi = phi,
                  v_e = v_e * (v_e_max - v_e_min - v_e) / (v_e_max - v_e_min),
                  eta = eta * (eta_max - eta_min - eta) / (eta_max - eta_min))
    
    ### Define priors
    # Gaussian prior with E(phi) = 1
    lambda_phi <- pi / 2
    log_likelihood_phi <- -0.5 * phi ^ 2 / lambda_phi
    d_log_likelihood_phi <- -phi / lambda_phi
    
    # Beta prior for v_e (alpha = beta = 2, E(v_e) = 1/2)
    a_v_e <- 2
    b_v_e <- 2
    log_likelihood_v_e <- (a_v_e - 1) * log(v_e) + (b_v_e - 1) * log(1 - v_e)
    d_log_likelihood_v_e <- (a_v_e - 1) / v_e - (b_v_e - 1) / (1 - v_e)
    
    # Beta prior for eta (alpha = beta = 2, E(eta) = 1/2)
    a_eta <- 2
    b_eta <- 2
    log_likelihood_eta <- (a_eta - 1) * log(eta) + (b_eta - 1) * log(1 - eta)
    d_log_likelihood_eta <- (a_eta - 1) / eta - (b_eta - 1) / (1 - eta)
    
    log_likelihood_params <- log_likelihood_phi + 
        log_likelihood_v_e + 
        log_likelihood_eta
    d_log_likelihood_params <- c(d_log_likelihood_phi, 
                                 d_log_likelihood_v_e, 
                                 d_log_likelihood_eta)
    
    ### start computing stuff
    # note that our parameter values for eta and v_e are between 0 and 1:
    #   i.e. they are relative to the variance in y
    #   the var_y_lib argument can be set on call to the this function (to 
    #     facilitate model comparisons with different settings), and with the 
    #     default to compute from y_lib directly
    
    eta_scaled <- eta * var_y_lib
    v_e_scaled <- v_e * var_y_lib
    
    # covariance calcs
    if (!is.null(x_pred)) # compute full distance matrix using lib and pred
    {
        dist_xy <- as.matrix(dist(rbind(x_lib, x_pred)))
        lib_idx <- 1:NROW(x_lib)
        pred_idx <- NROW(x_lib) + 1:NROW(x_pred)
        
        squared_dist_lib_lib <- dist_xy[lib_idx, lib_idx] ^ 2
        K_pred_pred <- eta_scaled * 
            exp(-phi ^ 2 * dist_xy[pred_idx, pred_idx] ^ 2)
        K_pred_lib <- eta_scaled * 
            exp(-phi ^ 2 * dist_xy[pred_idx, lib_idx] ^ 2)
    } else { # compute distance matrix using just lib
        squared_dist_lib_lib <- (as.matrix(dist(x_lib))) ^ 2
    }
    K_lib_lib <- eta_scaled * exp(-phi ^ 2 * squared_dist_lib_lib)
    Sigma <- K_lib_lib + v_e_scaled * diag(NROW(x_lib))
    
    # cholesky algorithm from Rasmussen & Williams (2006, algorithm 2.1)
    R <- chol(Sigma)
    alpha <- backsolve(R, forwardsolve(t(R), y_lib - mean_y))
    L_inv <- forwardsolve(t(R), diag(NROW(x_lib)))
    Sigma_inv <-  t(L_inv) %*% L_inv
    
    # marginal likelihood
    log_likelihood_lib <- (-0.5 * t(y_lib - mean_y) %*% alpha) - 
        sum(log(diag(R)))
    out$neg_log_likelihood <- -(log_likelihood_lib + log_likelihood_params)
    # out$neg_log_likelihood_L00 <- 
    #     0.5 * sum(log(diag(Sigma_inv))) -
    #     0.5 * sum(alpha ^ 2 / diag(Sigma_inv))
    
    ### Compute gradient if requested
    if (gradient)
    {
        # derivative of R
        # R = exp(-phi^2 * sq_dist)
        # d(R) = -2 * phi * sq_dist   *   exp(-phi^2 * sq_dist)
        #      = log(R) * 2 / phi     *   R
        d_R_d_p <- -2 * phi * squared_dist_lib_lib * (K_lib_lib / eta_scaled)
        W <- eta_scaled * d_R_d_p
        vQ <- alpha %*% t(alpha) - Sigma_inv
        d_log_likelihood_lib <- c(phi = 0.5 * sum(vQ * W), 
                                  v_e = 0.5 * sum(diag(vQ)), 
                                  eta = 0.5 * sum(vQ * K_lib_lib / eta_scaled))
        
        # J is gradient in parameter space: 
        # transform to get gradient in transformed parameters
        J <- d_log_likelihood_lib + d_log_likelihood_params
        out$gradient_neg_log_likelihood <- -J * d_params
    }
    
    ### Compute mean and variance for x_pred if given
    if (!is.null(x_pred))
    {
        covariance_matrix <- K_pred_pred - 
            K_pred_lib %*% Sigma_inv %*% t(K_pred_lib) + 
            v_e_scaled
        out$mean_pred <- mean_y + K_pred_lib %*% alpha
        out$pred_var <- diag(covariance_matrix)
        if (cov_matrix)
        {
            out$covariance_pred <- covariance_matrix
        }
    }
    
    return(out)
}