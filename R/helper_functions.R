rEDM_warning <- function(..., silent = FALSE)
{
    if(!silent)
        warning(..., call. = FALSE)
    return()
}

check_params_against_lib <- function(E, tau, tp, lib, silent = FALSE)
{
    vector_start <- max(-(E - 1) * tau, 0, tp)
    vector_end <- min(-(E - 1) * tau, 0, tp)
    vector_length <- abs(vector_start - vector_end) + 1
    
    max_lib_segment <- max(lib[, 2] - lib[, 1] + 1)
    if(vector_length > max_lib_segment)
    {
        rEDM_warning("Invalid parameter combination: E = ", E, ", tau = ", tau, 
                     ", tp = ", tp, silent = silent)
        return(FALSE)
    }
    return(TRUE)
}

convert_to_column_indices <- function(columns, block, silent = FALSE)
{
    if (is.numeric(columns))
    {
        if (any(columns > NCOL(block) | columns < 1))
            rEDM_warning("Some column indices exceeded the possible range ", 
                         "and were ignored.", silent = silent)
        out <- columns[columns <= NCOL(block) & columns >= 1]
    } else {
        indices <- match(columns, colnames(block))
        if (any(is.na(indices)))
            rEDM_warning("Some column names could not be matched and were ignored.", 
                         silent = silent)
        out <- indices[is.finite(indices)]
    }
    if(length(out) < 1)
    {
        stop("Columns requested were invalid, please check the ", 
             match.call()$columns, " argument.")
    }
    return(out)
}

coerce_lib <- function(lib, silent = FALSE) {
    if (is.vector(lib)) lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if (!all(lib[, 2] >= lib[, 1]))
        rEDM_warning("Some library rows look incorrectly formatted, please ", 
                "check the ", match.call()$lib, " argument.", silent = silent)
    return(lib)
}

setup_time_and_data_block <- function(model, first_column_time, block)
{
    if (first_column_time)
    {
        if (is.vector(block))
            time <- block
        else
        {
            time <- block[, 1]
            block <- block[, -1]
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
    model$set_time(time)
    model$set_block(data.matrix(block))
    return(block)
}

setup_model_flags <- function(model, exclusion_radius, epsilon, silent)
{
    # handle exclusion radius
    if (is.null(exclusion_radius))
    {
        exclusion_radius <- -1
    }
    model$set_exclusion_radius(exclusion_radius)
    
    # handle epsilon
    if (is.null(epsilon))
    {
        epsilon <- -1
    }
    model$set_epsilon(epsilon)
    
    # handle silent flag
    if (silent)
    {
        model$suppress_warnings()
    }
    return()
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