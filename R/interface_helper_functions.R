rEDM_warning <- function(..., silent = FALSE)
{
    if (!silent)
        warning(..., call. = FALSE)
    return()
}

check_params_against_lib <- function(E, tau, tp, lib, silent = FALSE)
{
    vector_start <- max(-(E - 1) * tau, 0, tp)
    vector_end <- min(-(E - 1) * tau, 0, tp)
    vector_length <- abs(vector_start - vector_end) + 1
    
    max_lib_segment <- max(lib[, 2] - lib[, 1] + 1)
    if (vector_length > max_lib_segment)
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
    if (length(out) < 1)
    {
        stop("in ", deparse(sys.call(-1)), ", Columns requested were invalid, please check the ", 
             match.call()$columns, " argument.", call. = FALSE)
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

setup_time_and_time_series <- function(time_series)
{
    if (is.mts(time_series)) {
        time_series <- time_series[, 1]
    }
    if (is.ts(time_series)) {
        time <- as.numeric(time(time_series))
        time_series <- as.numeric(time_series)
    } else if (is.vector(time_series)) {
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
    if (any(is.na(time)))
        time <- seq(length(time_series))
    return(list(time = time, 
                time_series = time_series))
}

setup_time_and_block <- function(block, first_column_time = FALSE)
{
    if (first_column_time)
    {
        if (is.mts(block)) # convert multivariate time series into matrix
        {
            class(block) <- "matrix"
        }
        
        if (is.vector(block) || is.ts(block))
        {
            stop("in ", deparse(sys.call(-1)), ", No data columns to work with if `first_column_time = TRUE`.", call. = FALSE)
        } else {# matrix or data.frame
            time <- block[, 1]
            block <- block[, -1, drop = FALSE]
        }
    }
    else
    {
        if (is.mts(block)) # use time index from multivariate time series
        {
            time <- time(block)
            class(block) <- "matrix"
            attr(block, "tsp") <- NULL
        } else {
            time <- rownames(block)
        }
    }
    if (is.null(time))
    {
        time <- seq(NROW(block))
    } else {
        time <- as.numeric(time)
    }
    if (any(is.na(time)))
        time <- seq(NROW(block))
    return(list(time = time, 
                block = data.matrix(block)))
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