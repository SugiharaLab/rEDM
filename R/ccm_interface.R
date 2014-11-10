#' Perform convergent cross mapping using simplex projection
#'
#' \code{ccm} uses time delay embedding on one time series to generate an 
#' attractor reconstruction, and then applies the simplex projection algorithm 
#' to estimate concurrent values of another time series. This method is 
#' typically applied, varying the library sizes, to determine if one time series
#' contains the necessary dynamic information to recover the influence of 
#' another, causal variable.
#' 
#' The default parameters are set so that passing a matrix as the only argument
#' will use E = 1 (embedding dimension), and leave-one-out cross-validation over
#' the whole time series to compute cross-mapping from the first column to the 
#' second column, letting the library size vary from 10 to 100 in increments of 
#' 10.
#' 
#' norm_type "L2 norm" (default) uses the typical Euclidean distance:
#' \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}}{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' norm_type "L1 norm" uses the Manhattan distance:
#' \deqn{distance(a,b) := \sum_i{|a_i - b_i|}}{distance(a, b) := \sum|a_i - b_i|}
#' 
#' @param block either a vector to be used as the time series, or a 
#'   data.frame or matrix where each column is a time series
#' @param lib a 2-column matrix (or 2-element vector) where each row specifes the 
#'   first and last *rows* of the time series to use for attractor reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param norm_type the distance function to use. see 'Details'
#' @param E the embedding dimensions to use for time delay embedding
#' @param tau the lag to use for time delay embedding
#' @param tp the prediction horizon (how far ahead to forecast)
#' @param num_neighbors the number of nearest neighbors to use (any of "e+1", 
#'   "E+1", "e + 1", "E + 1" will peg this parameter to E+1 for each run, any
#'   value < 1 will use all possible neighbors.)
#' @param lib_sizes the vector of library sizes to try
#' @param random_libs indicates whether to use randomly sampled libs
#' @param num_samples is the number of random samples at each lib size (this 
#'   parameter is ignored if random_libs is FALSE)
#' @param replace indicates whether to sample vectors with replacement
#' @param lib_column the index (or name) of the column to cross map from
#' @param target_column the index (or name) of the column to cross map to
#' @param first_column_time indicates whether the first column of the given 
#'   block is a time column (and therefore excluded when indexing)
#' @param exclusion_radius excludes vectors from the search space of nearest 
#'   neighbors if their *time index* is within exclusion_radius (NULL turns 
#'   this option off)
#' @param epsilon excludes vectors from the search space of nearest neighbors 
#'   if their *distance* is farther away than epsilon (NULL turns this option 
#'   off)
#' @param silent prevents warning messages from being printed to the R console
#' @return A data.frame with forecast statistics for the different parameter 
#'   settings:
#' \tabular{ll}{
#'   L \tab library length (number of vectors)\cr
#'   num_pred \tab number of predictions\cr
#'   rho \tab correlation coefficient between observations and predictions\cr
#'   mae \tab mean absolute error\cr
#'   rmse \tab root mean square error
#' }
#' @export 

ccm <- function(block, lib = c(1, NROW(block)), pred = c(1, NROW(block)), 
                norm_type = c("L2 norm", "L1 norm"), E = 1, tau = 1, 
                tp = 0, num_neighbors = "e+1", lib_sizes = seq(10, 100, by = 10), 
                random_libs = TRUE, num_samples = 100, replace = TRUE, 
                lib_column = 1, target_column = 2, first_column_time = FALSE, 
                exclusion_radius = NULL, epsilon = NULL, silent = FALSE)
{
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
        else
        {
            time <- block[,1]
            block <- block[,-1]
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
        if(any(is.na(time)))
            time <- 1:NROW(block)
    }
    model$set_time(time)
    model$set_block(as.matrix(block))
    model$set_lib_column(convert_to_column_indices(lib_column))
    model$set_target_column(convert_to_column_indices(target_column))
    
    # setup norm type
    model$set_norm_type(switch(match.arg(norm_type), "L2 norm" = 2, "L1 norm" = 1))
    
    # setup lib and pred ranges
    if (is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if (is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
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
    
    params <- data.frame(E, tau, tp, num_neighbors, lib_column, target_column)
    e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if (any(e_plus_1_index, na.rm = TRUE))
        params$num_neighbors <- params$E+1

    model$set_params(params$E, params$tau, params$tp, params$num_neighbors, 
                     random_libs, num_samples, replace)
    model$run()
    stats <- model$get_output()
    return(cbind(params, stats))
}

ccm_means <- function(ccm_df)
{
    ccm_means <- aggregate(ccm_df, by = list(ccm_df$lib_size), function(x) {mean(x, na.rm = TRUE)})
    return(ccm_means[,-1]) # drop Group.1 column
}