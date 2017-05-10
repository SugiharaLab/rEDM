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
#' norm type "P norm" uses the LP norm, generalizing the L1 and L2 norm to use $p$ as the exponent:
#' \deqn{distance(a,b) := \sum_i{(a_i - b_i)^p}^{1/p}}{distance(a, b) := (\sum(a_i - b_i)^p)^(1/p)}
#' 
#' @param block either a vector to be used as the time series, or a 
#'   data.frame or matrix where each column is a time series
#' @param lib a 2-column matrix (or 2-element vector) where each row specifes the 
#'   first and last *rows* of the time series to use for attractor reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param norm_type the distance function to use. see 'Details'
#' @param P the exponent for the P norm
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
#' @param RNGseed will set a seed for the random number generator, enabling 
#'   reproducible runs of ccm with randomly generated libraries
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
#' @examples
#' data("sardine_anchovy_sst")
#' anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
#'   lib_column = "anchovy", target_column = "np_sst", 
#'   lib_sizes = seq(10, 80, by = 10), num_samples = 100)
#' @export 
ccm <- function(block, lib = c(1, NROW(block)), pred = lib, 
                norm_type = c("L2 norm", "L1 norm", "LP norm"), P = 0.5, E = 1, 
                tau = 1, tp = 0, num_neighbors = "e+1", 
                lib_sizes = seq(10, 100, by = 10), random_libs = TRUE, 
                num_samples = 100, replace = TRUE, lib_column = 1, 
                target_column = 2, first_column_time = FALSE, RNGseed = NULL, 
                exclusion_radius = NULL, epsilon = NULL, silent = FALSE)
{
    convert_to_column_indices <- function(columns)
    {
        if(is.numeric(columns))
        {
            if(any(columns > NCOL(block)))
                warning("Some column indices exceed the number of columns and were ignored.")
            return(columns[columns <= NCOL(block)])
        }
        # else
        indices <- match(columns, col_names)
        if(any(is.na(indices)))
            warning("Some column names could not be matched and were ignored.")
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
    col_names <- colnames(block)
    model$set_time(time)
    model$set_block(data.matrix(block))
    model$set_lib_column(convert_to_column_indices(lib_column))
    model$set_target_columns(convert_to_column_indices(target_column))
    
    # setup norm type
    model$set_norm_type(switch(match.arg(norm_type), "P norm" = 3, "L2 norm" = 2, "L1 norm" = 1))
    model$set_p(P)
    
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
    else
        warning("Note: CCM results are typically interpreted in the opposite direction of causation. Please see 'Detecting causality in complex ecosystems' (Sugihara et al. 2012) for more details.")
    
    # check inputs?
    
    params <- data.frame(E, tau, tp, num_neighbors, lib_column)
    e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if (any(e_plus_1_index, na.rm = TRUE))
        params$num_neighbors <- params$E+1

    model$set_params(params$E, params$tau, params$tp, params$num_neighbors, 
                     random_libs, num_samples, replace)
    if(!is.null(RNGseed))
        model$set_seed(RNGseed)
    model$run()
    stats <- model$get_output()
    return(cbind(params, stats, row.names = NULL))
}

#' Take output from ccm and compute means as a function of library size.
#'
#' \code{ccm_means} is a utility function to summarize output from the \code{\link{ccm}} 
#' function
#' 
#' @param ccm_df a data.frame, usually output from the \code{\link{ccm}} function
#' @param FUN a function that aggregates the numerical statistics (by default, uses the mean)
#' @param ... optional arguments to FUN
#' @return A data.frame with forecast statistics aggregated at each unique library
#'   size
#' @examples 
#' data("sardine_anchovy_sst")
#' anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
#'   lib_column = "anchovy", target_column = "np_sst", 
#'   lib_sizes = seq(10, 80, by = 10), num_samples = 100)
#' a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
#' @export 
ccm_means <- function(ccm_df, FUN = mean, ...)
{
    lib <- ccm_df$lib_column[!duplicated(ccm_df$lib_size)]
    target <- ccm_df$target_column[!duplicated(ccm_df$lib_size)]
    ccm_df$lib_column <- NULL
    ccm_df$target_column <- NULL
    ccm_means <- aggregate(ccm_df, by = list(ccm_df$lib_size), FUN, ...)
    col_idx <- which(names(ccm_means) == "lib_size")
    ccm_means <- cbind(ccm_means[,1:(col_idx-1)], 
                       lib_column = lib, target_column = target, 
                       ccm_means[,col_idx:NCOL(ccm_means)])
    return(ccm_means[,-1]) # drop Group.1 column
}