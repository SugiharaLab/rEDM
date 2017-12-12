#' Randomization test for nonlinearity using S-maps and surrogate data
#' 
#' \code{\link{test_nonlinearity}} tests for nonlinearity using S-maps by 
#' comparing improvements in forecast skill (delta rho and delta mae) between 
#' linear and nonlinear models with a null distribution from surrogate data.
#' 
#' @param ts the original time series
#' @param method which algorithm to use to generate surrogate data
#' @param num_surr the number of null surrogates to generate
#' @param T_period the period of seasonality for seasonal surrogates 
#'   (ignored for other methods)
#' @param E the embedding dimension for s_map
#' @param ... optional arguments to s_map
#' @return A data.frame containing the following components:
#' \tabular{ll}{
#'   delta_rho \tab the value of the delta rho statistic\cr
#'   delta_mae \tab the value of the delat mae statistic\cr
#'   num_surr \tab the size of the null distribution\cr
#'   delta_rho_p_value \tab the p-value for delta rho\cr
#'   delta_mae_p_value \tab the p-value for delta mae\cr
#' }
#' 
test_nonlinearity <- function(ts, method = "ebisuzaki", num_surr = 200, 
                              T_period = 1, E = 1, ...)
{
    compute_stats <- function(ts, ...)
    {
        results <- s_map(ts, stats_only = TRUE, silent = TRUE, ...)
        delta_rho <- max(results$rho) - results$rho[results$theta == 0]
        delta_mae <- results$mae[results$theta == 0] - min(results$mae)
        return(c(delta_rho = delta_rho, delta_mae = delta_mae))
    }
    
    actual_stats <- compute_stats(ts, ...)
    delta_rho <- actual_stats["delta_rho"]
    delta_mae <- actual_stats["delta_mae"]
    names(delta_rho) <- NULL
    names(delta_mae) <- NULL
    surrogate_data <- make_surrogate_data(ts, method, num_surr, T_period)
    null_stats <- data.frame(t(apply(surrogate_data, 2, compute_stats, ...)))
    
    return(data.frame(
        delta_rho = delta_rho, 
        delta_mae = delta_mae, 
        num_surr = num_surr, 
        E = E, 
        delta_rho_p_value = (sum(null_stats$delta_rho > delta_rho) + 1) / 
            (num_surr + 1), 
        delta_mae_p_value = (sum(null_stats$delta_mae > delta_mae) + 1) / 
            (num_surr + 1)))
}


#' Generate surrogate data for permutation/randomization tests
#'
#' \code{\link{make_surrogate_data}} generates surrogate data under several different 
#' null models.
#' 
#' Method "random_shuffle" creates surrogates by randomly permuting the values 
#' of the original time series.
#' 
#' Method "Ebisuzaki" creates surrogates by randomizing the phases of a Fourier 
#' transform, preserving the power spectra of the null surrogates.
#' 
#' Method "seasonal" creates surrogates by computing a mean seasonal trend of 
#' the specified period and shuffling the residuals.
#' 
#' See \code{\link{test_nonlinearity}} for context.
#' 
#' @param ts the original time series
#' @param method which algorithm to use to generate surrogate data
#' @param num_surr the number of null surrogates to generate
#' @param T_period the period of seasonality for seasonal surrogates 
#'   (ignored for other methods)
#' @return A matrix where each column is a separate surrogate with the same 
#'   length as `ts`.
#' @examples
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' make_surrogate_data(ts, method = "ebisuzaki")
#' 
make_surrogate_data <- function(ts, method = c("random_shuffle", "ebisuzaki", 
                                               "seasonal"), 
                                num_surr = 100, T_period = 1)
{  
    method <- match.arg(method)
    if (method == "random_shuffle")
    {
        return(sapply(1:num_surr, function(i) {
            sample(ts, size = length(ts))
        }))
    }
    else if (method == "ebisuzaki")
    {
        if (any(!is.finite(ts)))
            stop("input time series contained invalid values")
        
        n <- length(ts)
        n2 <- floor(n / 2)
        
        sigma <- sd(ts)
        a <- fft(ts)
        amplitudes <- abs(a)
        amplitudes[1] <- 0
        
        return(sapply(1:num_surr, function(i) {
            if (n %% 2 == 0) # even length
            {
                thetas <- 2 * pi * runif(n2 - 1)
                angles <- c(0, thetas, 0, -rev(thetas))
                recf <- amplitudes * exp(complex(imaginary = angles))
                recf[n2] <- complex(real = sqrt(2) * amplitudes[n2] * 
                                        cos(runif(1) * 2 * pi))
            }
            else # odd length
            {
                thetas <- 2 * pi * runif(n2)
                angles <- c(0, thetas, -rev(thetas))
                recf <- amplitudes * exp(complex(imaginary = angles))
            }
            temp <- Re(fft(recf, inverse = T) / n)
            
            # adjust variance of the surrogate time series to match original
            return(temp / sd(temp) * sigma)
        }))
    }
    else # method = "seasonal"
    {
        if (any(!is.finite(ts)))
            stop("input time series contained invalid values")
                
        n <- length(ts)
        I_season <- suppressWarnings(matrix(1:T_period, nrow = n, ncol = 1))
        
        # Calculate seasonal cycle using smooth.spline
        seasonal_F <- smooth.spline(c(I_season - T_period, I_season, 
                                      I_season + T_period), 
                                    c(ts, ts, ts))
        seasonal_cyc <- predict(seasonal_F, I_season)$y
        seasonal_resid <- ts - seasonal_cyc
        
        return(sapply(1:num_surr, function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        }))
    }
}
