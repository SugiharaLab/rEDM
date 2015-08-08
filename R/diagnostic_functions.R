#' Generate surrogate data for permutation/randomization tests
#'
#' \code{make_surrogate data} generates surrogate data under several different 
#' null models.
#' 
#' Method "random_shuffle" creates surrogates by randomly permuting the values 
#' of the original time series.
#' 
#' Method "Ebisuzaki" creates surrogates by randomizing the phases of a Fourier 
#' transform, preserving the power spectra of the null surrogates.
#' 
#' See \code{test_nonlinearity} for context.
#' 
#' @param ts the original time series
#' @param method which algorithm to use to generate surrogate data
#' @param num_surr the number of null surrogates to generate
#' @param T_period the period of seasonality for seasonal surrogates (ignored for other methods)
#' @return A matrix where each column is a separate surrogate with the same length as \code{ts}.
#' @export 
#' 
make_surrogate_data <- function(ts, method = c("random_shuffle", "ebisuzaki", "seasonal"), 
                                num_surr = 100, T_period = 1)
{  
    method <- match.arg(method)
    if(method == "random_shuffle")
    {
        return(sapply(1:num_surr, function(i) {
            sample(ts, size = length(ts))
        }))
    }
    else if(method == "ebisuzaki")
    {
        if(any(!is.finite(ts)))
            stop("input time series contained invalid values")
        
        n <- length(ts)
        n2 <- floor(n/2)
        
        mu <- mean(ts)
        sigma <- sd(ts)
        a <- fft(ts)
        amplitudes <- abs(a)
        amplitudes[1] <- 0
        
        return(sapply(1:num_surr, function(i) {
            if(n %% 2 == 0) # even length
            {
                thetas <- 2*pi*runif(n2-1)
                angles <- c(0, thetas, 0, -rev(thetas))
                recf <- amplitudes * exp(complex(imaginary = angles))
                recf[n2] <- complex(real = sqrt(2) * amplitudes[n2] * cos(runif(1)*2*pi))
            }
            else # odd length
            {
                thetas <- 2*pi*runif(n2)
                angles <- c(0, thetas, -rev(thetas))
                recf <- amplitudes * exp(complex(imaginary = angles))
            }
            temp <- Re(fft(recf, inverse = T) / n)
            
            # adjust variance of the surrogate time series to match the original            
            return(temp / sd(temp) * sigma)
        }))
        
        return(output)
    }
    else
    {
        if(any(!is.finite(ts)))
            stop("input time series contained invalid values")
                
        n <- length(ts)
        I_season <- suppressWarnings(matrix(1:T_period, nrow=n, ncol=1))
        
        # Calculate seasonal cycle using smooth.spline
        seasonal_F <- smooth.spline(rbind(I_season - T_period, I_season, I_season + T_period), 
                                    rbind(ts, ts, ts))
        seasonal_cyc <- predict(seasonal_F,I_season)$y
        seasonal_resid <- ts - seasonal_cyc
        
        return(sapply(1:num_surr, function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        }))
    }
}


# test for cross map convergence with library size
# equivalent of ccmtest from multispatialCCM
test_convergence <- function(ccm_results)
{
    return()
}
