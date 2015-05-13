rm(list = ls())
library(rEDM)
load("vostok_data.Rdata")
vostok_data <- rev(vostok_data)

n <- NROW(vostok_data)
vostok_data$temp <- vostok_data$temp + 30
vostok_data$sqrt_temp <- sqrt(vostok_data$temp)
lib <- c(1, n)
pred <- c(1, n)
if(TRUE)
{
    ccm_results <- ccm(vostok_data, lib = lib, pred = pred,  E = 3,
                       lib_sizes = n, random_libs = FALSE, num_samples = 1, 
                       lib_column = "CO2", target_column = c("temp", "sqrt_temp"))
}

# ccm params
if(FALSE)
{
    block <- vostok_data
    norm_type <- "L2 norm"
    E <- 3
    tau <- 1
    tp <- 0
    num_neighbors <- "e+1"
    lib_sizes <- n
    random_libs <- FALSE
    num_samples <- 1
    replace <- FALSE
    lib_column <- "CO2"
    target_column <- c("temp", "sqrt_temp")
    first_column_time <- FALSE
    exclusion_radius <- NULL
    epsilon <- NULL
    silent <- FALSE
    
    convert_to_column_indices <- function(columns)
    {
        if(is.numeric(columns))
        {
            if(any(columns > NCOL(block)))
                message("Warning: some column indices exceed the number of columns and were ignored.")
            return(columns[columns <= NCOL(block)])
        }
        # else
        indices <- match(columns, col_names)
        if(any(is.na(indices)))
            message("Warning: some column names could not be matched and were ignored.")
        return(indices[is.finite(indices)])
    }
    
    # make new model object
    model <- new(Xmap)
    
    # setup data
    if(first_column_time)
    {
        if(is.vector(block))
            time <- block
        else {
            time <- block[,1]
            block <- block[,-1]
        }
    } else {
        time <- rownames(block)
    }
    if (is.null(time))
    {
        time <- 1:NROW(block)
    } else {
        time <- as.numeric(time)
        if(any(is.na(time)))
            time <- 1:NROW(block)
    }
    col_names <- colnames(block)
    model$set_time(time)
    model$set_block(data.matrix(block))
    model$set_lib_column(convert_to_column_indices(lib_column))
    model$set_target_columns(convert_to_column_indices(target_column))
    
    # setup norm type
    model$set_norm_type(2)
    
    # setup lib and pred ranges
    if (is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if (is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    
    if(!all(lib[,2] >= lib[,1]))
        warning("Some library rows look incorrectly formatted, please check the lib argument.")
    if(!all(pred[,2] >= pred[,1]))
        warning("Some library rows look incorrectly formatted, please check the pred argument.")
    
    model$set_lib(lib)
    model$set_pred(pred)
    model$set_lib_sizes(lib_sizes)
    
    # handle exclusion radius
    if (is.null(exclusion_radius))
        exclusion_radius = -1;
    model$set_exclusion_radius(exclusion_radius)
    
    # TODO: handle epsilon
    
    # handle silent flag
    if (silent)
        model$suppress_warnings()
    
    # check inputs?
    
    params <- data.frame(E, tau, tp, num_neighbors, lib_column)
    e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if (any(e_plus_1_index, na.rm = TRUE))
        params$num_neighbors <- params$E+1
    
    model$set_params(params$E, params$tau, params$tp, params$num_neighbors, 
                     random_libs, num_samples, replace)
    message("AAAA")
    model$run()
    message("BBBB")
    stats <- model$get_output()
    results <- cbind(params, stats)
}
