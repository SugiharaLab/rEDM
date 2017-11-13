convert_to_column_indices <- function(columns, block)
{
    if (is.numeric(columns))
    {
        if (any(columns > NCOL(block)))
            warning("Some column indices exceed the number of columns ", 
                    "and were ignored.")
        return(columns[columns <= NCOL(block)])
    }
    # else
    indices <- match(columns, colnames(block))
    if (any(is.na(indices)))
        warning("Some column names could not be matched and were ignored.")
    return(indices[is.finite(indices)])
}

setup_lib_and_pred <- function(model, lib, pred)
{
    if (is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if (is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    
    if (!all(lib[, 2] >= lib[, 1]))
        warning("Some library rows look incorrectly formatted, please check ", 
                "the lib argument.")
    if (!all(pred[, 2] >= pred[, 1]))
        warning("Some library rows look incorrectly formatted, please check ", 
                "the pred argument.")
    
    model$set_lib(lib)
    model$set_pred(pred)
    return()
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