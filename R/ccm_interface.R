#' Perform convergent cross mapping using simplex projection
#'
#' \code{\link{ccm}} uses time delay embedding on one time series to generate an 
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
#' \code{norm = 2} (default) uses the "L2 norm", Euclidean distance:
#'   \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
#'     }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' \code{norm = 1} uses the "L1 norm", Manhattan distance:
#'   \deqn{distance(a,b) := \sum_i{|a_i - b_i|}
#'     }{distance(a, b) := \sum|a_i - b_i|}
#' Other values generalize the L1 and L2 norm to use the given argument as the 
#'   exponent, P, as:
#'   \deqn{distance(a,b) := \sum_i{(a_i - b_i)^P}^{1/P}
#'     }{distance(a, b) := (\sum(a_i - b_i)^P)^(1/P)}
#' 
#' @inheritParams block_lnlp
#' @inheritParams simplex
#' @param lib_sizes the vector of library sizes to try
#' @param random_libs indicates whether to use randomly sampled libs
#' @param num_samples is the number of random samples at each lib size (this 
#'   parameter is ignored if random_libs is FALSE)
#' @param replace indicates whether to sample vectors with replacement
#' @param lib_column the index (or name) of the column to cross map from
#' @param RNGseed will set a seed for the random number generator, enabling 
#'   reproducible runs of ccm with randomly generated libraries
#' @return A data.frame with forecast statistics for the different parameter 
#'   settings:
#' \tabular{ll}{
#'   \code{L} \tab library length (number of vectors)\cr
#'   \code{num_pred} \tab number of predictions\cr
#'   \code{rho} \tab correlation coefficient between observations and predictions\cr
#'   \code{mae} \tab mean absolute error\cr
#'   \code{rmse} \tab root mean square error
#' }
#' @examples
#' data("sardine_anchovy_sst")
#' anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
#'   lib_column = "anchovy", target_column = "np_sst", 
#'   lib_sizes = seq(10, 80, by = 10), num_samples = 100)
#'  
ccm <- function(block, lib = c(1, NROW(block)), pred = lib, 
                norm = 2, E = 1, 
                tau = 1, tp = 0, num_neighbors = "e+1", 
                lib_sizes = seq(10, 100, by = 10), random_libs = TRUE, 
                num_samples = 100, replace = TRUE, lib_column = 1, 
                target_column = 2, first_column_time = FALSE, RNGseed = NULL, 
                exclusion_radius = NULL, epsilon = NULL, silent = FALSE)
{
    # make new model object
    model <- new(Xmap)
    
    # setup data
    block <- setup_time_and_data_block(model, first_column_time, block)
    my_lib_column <- convert_to_column_indices(lib_column, block, 
                                               silent = silent)
    model$set_lib_column(my_lib_column)
    my_target_column <- convert_to_column_indices(target_column, block, 
                                                  silent = silent)
    model$set_target_column(my_target_column)
    
    # setup norm type
    model$set_norm(norm)
    
    # setup lib and pred ranges
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)
    model$set_lib(lib)
    model$set_pred(pred)
    
    # check lib_sizes
    prev_num_lib_sizes <- length(lib_sizes)
    lib_sizes <- lib_sizes[lib_sizes >= 0]
    lib_sizes <- unique(sort(lib_sizes))
    if (length(lib_sizes) < 1)
        stop("No valid lib sizes found among input", lib_sizes)
    if (length(lib_sizes) < prev_num_lib_sizes)
        rEDM_warning("Some requested lib sizes were redundant or bad and ignored.", 
                     silent = silent)
    model$set_lib_sizes(lib_sizes)
    
    # handle exclusion radius
    if (is.null(exclusion_radius))
    {
        exclusion_radius <- -1
    }
    model$set_exclusion_radius(exclusion_radius)
    
    # TODO: handle epsilon
    
    # handle silent flag
    if (silent)
    {
        model$suppress_warnings()
    }
    rEDM_warning("Note: CCM results are typically interpreted in the opposite ", 
                 "direction of causation. Please see 'Detecting causality in ", 
                 "complex ecosystems' (Sugihara et al. 2012) for more details.", 
                 silent = silent)
    
    params <- data.frame(E, tau, tp, nn = num_neighbors, lib_column, target_column)
    e_plus_1_index <- match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1"))
    if (any(e_plus_1_index, na.rm = TRUE))
        params$nn <- params$E + 1
    params$nn <- as.numeric(params$nn)
    
    if (!check_params_against_lib(params$E, params$tau, params$tp, lib, 
                                  silent = silent))
    {
        stop("Parameter combination was invalid, stopping.")
    }
    
    model$set_params(params$E, params$tau, params$tp, params$nn, 
                     random_libs, num_samples, replace)
    if (!is.null(RNGseed))
        model$set_seed(RNGseed)
    model$run()
    stats <- model$get_output()
    return(cbind(params, stats, row.names = NULL))
}

#' Take output from ccm and compute means as a function of library size.
#'
#' \code{\link{ccm_means}} is a utility function to summarize output from the 
#'   \code{\link{ccm}} function
#' 
#' @param ccm_df a data.frame, usually output from the \code{\link{ccm}} 
#'   function
#' @param FUN a function that aggregates the numerical statistics (by default, 
#'   uses the mean)
#' @param ... optional arguments to FUN
#' @return A data.frame with forecast statistics aggregated at each unique 
#'   library size
#' @examples 
#' data("sardine_anchovy_sst")
#' anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
#'   lib_column = "anchovy", target_column = "np_sst", 
#'   lib_sizes = seq(10, 80, by = 10), num_samples = 100)
#' a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
#'  
ccm_means <- function(ccm_df, FUN = mean, ...)
{
    lib <- ccm_df$lib_column[!duplicated(ccm_df$lib_size)]
    target <- ccm_df$target_column[!duplicated(ccm_df$lib_size)]
    ccm_df$lib_column <- NULL
    ccm_df$target_column <- NULL
    ccm_means <- aggregate(ccm_df, by = list(ccm_df$lib_size), FUN, ...)
    col_idx <- which(names(ccm_means) == "lib_size")
    ccm_means <- cbind(ccm_means[, 1:(col_idx - 1)], 
                       lib_column = lib, target_column = target, 
                       ccm_means[, col_idx:NCOL(ccm_means)])
    return(ccm_means[, -1]) # drop Group.1 column
}