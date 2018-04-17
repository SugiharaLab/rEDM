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

#' Make a lagged block for multiview
#'
#' \code{\link{make_block}} generates a lagged block with the appropriate max_lag and 
#' tau, while respecting lib (by inserting NANs, when trying to lag past lib 
#' regions)
#' 
#' @param block a data.frame or matrix where each column is a time series
#' @param t the time index for the block
#' @param max_lag the total number of lags to include for each variable. So if 
#'   max_lag == 3, a variable X will appear with lags X[t], X[t - tau], 
#'   X[t - 2*tau]
#' @param tau the lag to use for time delay embedding
#' @param lib a 2-column matrix (or 2-element vector) where each row specifies 
#'   the first and last *rows* of the time series to use for attractor 
#'   reconstruction
#' @param restrict_to_lib whether to restrict the final lagged block to 
#'   just the rows specified in lib (if lib exists)
#' @return A data.frame with the lagged columns and a time column. If the 
#'   original block had columns X, Y, Z and max_lag = 3, then the returned 
#'   data.frame will have columns TIME, X, X_1, X_2, Y, Y_1, Y_2, Z, Z_1, Z_2.
#' @examples 
#' data("block_3sp")
#' make_block(block_3sp[, c(2, 5, 8)])
make_block <- function(block, t = NULL, max_lag = 3, tau = 1, lib = NULL, 
                       restrict_to_lib = TRUE)
{
    # be sure to convert block if input is a vector
    if(is.vector(block))
        block <- matrix(block, ncol = 1)
    num_vars <- NCOL(block)
    num_rows <- NROW(block)
    
    # coerce lib into matrix if necessary
    if (!is.null(lib))
    {
        if (is.vector(lib)) lib <- matrix(lib, ncol = 2, byrow = TRUE)
    }
    # output is the returned data frame
    output <- matrix(NA, nrow = num_rows, ncol = 1 + num_vars * max_lag)
    col_names <- character(1 + num_vars * max_lag)
    
    # create the time column
    if (is.null(t))
        output[, 1] <- 1:num_rows
    else
        output[, 1] <- t
    col_names[1] <- "time"
    
    # add max_lag lags for each column in block
    col_index <- 2
    if (is.null(colnames(block)))
        colnames(block) <- paste0("col", seq_len(num_vars))
    for (j in 1:num_vars)
    {
        ts <- block[, j]
        output[, col_index] <- ts
        col_names[col_index] <- colnames(block)[j]
        col_index <- col_index + 1
        
        ## add lags if required
        if (max_lag > 1)
        {
            for (i in 1:(max_lag - 1))
            {
                ts <- c(rep_len(NA, tau), ts[1:(num_rows - tau)])
                
                # make sure we pad beginning of lib segments with tau x NAs
                if (!is.null(lib))
                {
                    for (k in 1:NROW(lib))
                    {
                        ts[lib[k, 1] - 1 + (1:tau)] <- NA
                    }
                }
                output[, col_index] <- ts
                col_names[col_index] <- paste0(colnames(block)[j], "_", i * tau)
                col_index <- col_index + 1
            }
        }
    }
    
    # only take rows of lib if appropriate
    if (!is.null(lib) && restrict_to_lib)
    {
        row_idx <- sort(unique(
            do.call(c, mapply(seq, lib[, 1], lib[, 2], SIMPLIFY = FALSE))
        ))
        output <- output[row_idx, ]
    }
    
    output <- data.frame(output)
    names(output) <- col_names
    return(output)
}
