#' @name make_surrogate_data
#' @aliases make_surrogate_shuffle make_surrogate_ebisuzaki 
#'   make_surrogate_seasonal
#' 
#' @title Generate surrogate data for permutation/randomization tests
#'
#' @description This is a wrapper function for generating surrogate time series 
#'   using several different null models.
#' 
#' See \code{\link{test_nonlinearity}} for an example context for usage of this 
#' function.
#' 
#' @param ts the original time series
#' @param method which algorithm to use to generate surrogate data
#' @param num_surr the number of null surrogates to generate
#' @param ... remaining arguments are passed on to the specific function to 
#'   make surrogate time series of that type
#' @return A matrix where each column is a separate surrogate with the same 
#'   length as `ts`.
#' @examples
#' data("two_species_model")
#' ts <- two_species_model$x[1:200]
#' make_surrogate_data(ts, method = "ebisuzaki")
#' 
make_surrogate_data <- function(ts, method = c("random_shuffle", "ebisuzaki", 
                                               "seasonal", "twin"), 
                                num_surr = 100, ...)
{  
    stopifnot(num_surr > 0)
    
    method <- match.arg(method)
    switch(method, 
           random_shuffle = make_surrogate_shuffle(ts, num_surr, ...), 
           ebisuzaki = make_surrogate_ebisuzaki(ts, num_surr, ...), 
           seasonal = make_surrogate_seasonal(ts, num_surr, ...), 
           twin = make_surrogate_twin(ts, num_surr, ...))
}

#' @rdname make_surrogate_data
#'
#' @description \code{make_surrogate_shuffle()} creates surrogates by randomly 
#'   permuting the values of the original time series. 
#'
#' @inheritParams make_surrogate_data
#'
#' @examples
#' make_surrogate_shuffle(rnorm(100), 10)
#' 
make_surrogate_shuffle <- function(ts, num_surr = 100)
{
    matrix(unlist(
        lapply(seq(num_surr), function(i) {
            sample(ts, size = length(ts))
        })
    ), ncol = num_surr)
}

#' @rdname make_surrogate_data
#'
#' @description \code{make_surrogate_ebisuzaki()} creates surrogates by 
#'   randomizing the phases of a Fourier transform, preserving the power 
#'   spectra of the null surrogates
#'
#' @inheritParams make_surrogate_data
#'
#' @examples
#' make_surrogate_ebisuzaki(rnorm(100), 10)
#' 
make_surrogate_ebisuzaki <- function(ts, num_surr = 100)
{
    if (any(!is.finite(ts)))
        stop("input time series contained invalid values")
    
    n <- length(ts)
    n2 <- floor(n / 2)
    
    sigma <- sd(ts)
    a <- fft(ts)
    amplitudes <- abs(a)
    amplitudes[1] <- 0
    
    matrix(unlist(
        lapply(seq(num_surr), function(i) {
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
        temp <- Re(fft(recf, inverse = TRUE) / n)
        
        # adjust variance of the surrogate time series to match original
        return(temp / sd(temp) * sigma)
        })
    ), ncol = num_surr)
}

#' @rdname make_surrogate_data
#'
#' @description \code{make_surrogate_seasonal()} creates surrogates by 
#'   computing a mean seasonal trend of the specified period and shuffling the 
#'   residuals.
#'
#' @inheritParams make_surrogate_data
#' @param T_period the period of seasonality for seasonal surrogates 
#'   (ignored for other methods)
#'
#' @examples
#' make_surrogate_seasonal(rnorm(100) + sin(1:100 * pi / 6), 10)
#' 
make_surrogate_seasonal <- function(ts, num_surr = 100, T_period = 12)
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
    
    matrix(unlist(
        lapply(seq(num_surr), function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        })
    ), ncol = num_surr)
}

#' @title Generate a time index for a twin surrogate
#' @description Use the `twins` recurrence structure to construct a surrogate 
#'   time series. Note that we don't use the values of the actual time series, 
#'   but instead use the time series index (the surrogate is a reordered set of 
#'   points, with possible duplicates). This is easier for the code, as the 
#'   method is agnostic to the actual time series anyway, since it only cares 
#'   about the twins.
#' @details Suppose that the original time series is x, and the surrogate is s. 
#'   The algorithm supposes that the j-th value in s is the m-th value of x. 
#'   For the initial point (j = 1), we sample from either the twins of x, or 
#'   the same phase point in other cycles (depending on the `initial_point` 
#'   argument).
#'   Then, at each next value of j, we sample from the possible twins of x(m), 
#'   with the following conditions:
#'     (i) the selected twin of m cannot be the last point in the time series
#'     (ii) if `initial_point = "same_season"`, then the selected twin of m 
#'          cannot be j
#'   If these are met, then s(j+1) = x(m+1), and we continue.
#' @param twins a list of the twins. The list has length equal to the time 
#'   series, and each element is a vector of the candidate twins. See 
#'   \code{\link{identify_twins}}.
#' @inheritParams make_surrogate_twin
#' @return A vector of the same length as the original time series (here, 
#'   `length(twins)`) and containing the reordered time indices.
#'   
make_twin_idx <- function(twins, 
                          phase_lock = TRUE, 
                          T_period = 24, 
                          initial_point = "same_season")
{
    # sample from twins data structure to get next point
    get_next_idx <- function(idx, pos)
    {
        if (idx == 0) return(0) # stop if out of next points already
        candidates <- twins[[idx]]
        
        # make sure we don't choose the end of the time series
        candidates <- candidates[candidates < ts_length]
        
        # make sure we don't resync to the original time series
        if (phase_lock && 
            (initial_point == "same_season"))
        {
            candidates <- candidates[candidates != pos]
        }
        
        # return 0 if no valid candidates
        if (length(candidates) < 1) return(0)
        
        # otherwise, sample from one of the valid candidates
        return(candidates[sample.int(length(candidates), 1)] + 1)
    }
    
    ts_length <- length(twins)
    surr <- rep.int(0, ts_length)
    stopifnot(ts_length > 1)
    
    # select the initial point of the surrogate
    #   either from points occurring in the same season, 
    #   or twins of the first point
    if (phase_lock)
    {
        if (initial_point == "same_season")
        {
            candidates <- seq(1 + T_period, ts_length - 1, by = T_period)
        } else {
            candidates <- twins[[1]]
        }
    } else {
        candidates <- seq(1, ts_length - 1)
    }
    surr[1] <- candidates[sample.int(length(candidates), 1)]
    
    # build surrogate by sampling for the next point
    for (j in 1:(ts_length - 1))
    {
        surr[j + 1] <- get_next_idx(surr[j], j)
    }
    return(surr)
}

#' @title Construct twins based on the recurrence structure
#' @description Identify similar points in the `original_e` data structure. 
#'   Use quantiles to identify close enough points. Then 
#' @details The algorithm is as follows:
#'   (1) create recurrence matrix: values are 1 if their distance (using the 
#'      "maximum" measure) is lower than a specific quantile threshold
#'   (2) create twins if the columns of the recurrence matrix are identical 
#'       (and if the points are in the same phase, if `phase_lock = TRUE`)
#'   (3) check if the number of twins is satisfactor. If not, repeat and use 
#'       the next value of the quantile threshold
#' @param block the multivariate time series block. Each row is a data point, 
#'   and each column is a coordinate. We expect this to be either a lagged 
#'   block from a single time series or multiple time series.
#' @inheritParams make_surrogate_twin
#' @param quantile_vec the quantiles used to filter candidate twins.
#' @param min_num_twins how many twins are necessary to stop adjusting the 
#'   threshold for twins
#' @return A list of the twins. The list has length equal to the time 
#'   series, and each element is a vector of the candidate twins. 
#'   
identify_twins <- function(block, 
                           phase_lock = TRUE, 
                           T_period = 24, 
                           quantile_vec = c(0.125, 0.12, 0.11, 0.10, 0.09, 0.08, 
                                            0.07, 0.06, 0.05, 0.15, 0.16, 0.17, 
                                            0.18, 0.19, 0.20, 0.04), 
                           min_num_twins = 10)
{
    dist_mat <- as.matrix(dist(block, method = "maximum"))
    
    for (s in quantile_vec)
    {
        # make recurrence matrix
        threshold <- quantile(dist_mat, s)
        recurrence_matrix <- 0 + (dist_mat > threshold)
        
        # generate twins
        twins <- which(as.matrix(dist(t(recurrence_matrix), "maximum")) == 0,  
                       arr.ind = TRUE)
        if (phase_lock)
        {
            twins <- twins[(twins[, 1] - twins[, 2]) %% T_period == 0, ]
        }
        twins <- twins[order(twins[, 1]), ]
        num_twins <- NROW(twins) - NROW(recurrence_matrix)
        
        # check for enough twins and structure of output
        if (num_twins >= min_num_twins)
        {
            twins <- split(twins[, 2], twins[, 1])
            if (!isTRUE(all.equal(as.numeric(names(twins)), seq(length(twins)))))
            {
                stop("`twins` did not have the correct structure:\n", 
                     "length(twins) = ", length(twins), "\n", 
                     "names(twins) = ", names(twins))
            }
            break()
        }
    }
    
    if (num_twins < min_num_twins)
    {
        stop("Did not find enough twins after exhausting all quantile thresholds.\n", 
             "Wanted at least ", min_num_twins, " twins.")
    }
    return(twins)
}

#' @rdname make_surrogate_data
#'
#' @description \code{make_surrogate_twin()} creates surrogates using the twin-
#'   surrogate method, with the option to preserve the phase for seasonal/
#'   periodic data
#'
#' @inheritParams make_surrogate_seasonal
#' @param dim the embedding dimension for the state-space reconstruction, in 
#'   which twins are identified
#' @param tau the lag for the state-space reconstruction
#' @param phase_lock whether twins have to occur at the same phase
#' @param initial_point how to sample the initial point. If `"same_season"`, 
#'   then the initial point is chosen from the same phase in a different cycle, 
#'   and the surrogate is not allowed to line up in both phase and cycle with 
#'   the original time series.
#' @inheritParams identify_twins
#'
#' @examples
#' make_surrogate_twin(rnorm(100) + sin(1:100 * pi / 6), 10)
#' 
make_surrogate_twin <- function(ts,
                                 num_surr = 1,
                                 dim = 1, tau = 1, 
                                 phase_lock = TRUE, 
                                 T_period = 24, 
                                 initial_point = "same_season", 
                                 ...)
{
    # generate time-lag embedding matrix
    if (dim > 1) {
        block <- rEDM::make_block(ts, max_lag = dim, tau = tau)
        block <- block[-seq((dim - 1) * tau), -1]
        block <- block[, rev(seq(NCOL(block)))]
        block <- as.matrix(block)
    } else if (dim == 1) {
        block <- as.matrix(ts)
    } else {
        warning("Embedding dimension should be >= 1, given ", dim)
    }
    
    # find plausible twins for the surrogate method
    twins <- identify_twins(block, phase_lock = phase_lock, 
                            T_period = T_period, ...)
    
    # generate twin surrogates
    surrogates <- list()
    num_iter <- 0
    while ((length(surrogates) < num_surr) && # stop loop when we have enough 
           num_iter < (30 * num_surr))        # surrogates or reached max num iter
    {
        num_iter <- num_iter + 1
        surr <- make_twin_idx(twins, phase_lock = phase_lock, 
                              T_period = T_period, 
                              initial_point = initial_point)
        
        # if the surrogate is valid, not too short, then append it
        if (tail(surr, 1) != 0) {
            ts <- block[surr, dim]
            if (dim >= 2)
                ts <- c(block[surr[1], 1:(dim - 1)], ts)
            surrogates <- c(surrogates, list(ts))
        }
    }
    surrogates <- do.call(cbind, surrogates)
    rownames(surrogates) <- NULL
    return(surrogates)
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
        if(is.list(ts))
        {
            ts <- unlist(ts)
        }
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
