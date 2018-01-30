#' Perform forecasting using multiview embedding
#'
#' \code{\link{multiview}} applies the method described in Ye & Sugihara (2016) for 
#'   forecasting, wherein multiple attractor reconstructions are tested, and a 
#'   single nearest neighbor is selected from each of the top "k" 
#'   reconstructions to produce final forecasts.
#' 
#'   uses multiple time series given as input to generate an attractor 
#'   reconstruction, and then applies the simplex projection or s-map algorithm 
#'   to make forecasts. This method generalizes the \code{\link{simplex}} and 
#'   \code{\link{s-map}} routines, and allows for "mixed" embeddings, where multiple 
#'   time series can be used as different dimensions of an attractor 
#'   reconstruction.
#' 
#'   The default parameters are set so that, given a matrix of time series, 
#'   forecasts will be produced for the first column. By default, all possible 
#'   combinations of the columns are used for the attractor construction, the 
#'   k = sqrt(m) heuristic will be used, forecasts will be one time step ahead. 
#'   Rownames will be converted to numeric if possible to be used as the time 
#'   index, otherwise 1:NROW will be used instead. The default lib and pred are 
#'   to use the first half of the data for the "library" and to predict over the
#'   second half of the data. Unless otherwise set, the output will be just the 
#'   forecast statistics.
#' 
#' norm_type "L2 norm" (default) uses the typical Euclidean distance:
#' \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
#' }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' norm_type "L1 norm" uses the Manhattan distance:
#' \deqn{distance(a,b) := \sum_i{|a_i - b_i|}
#' }{distance(a, b) := \sum|a_i - b_i|}
#' norm type "P norm" uses the P norm, generalizing the L1 and L2 norm to use 
#'   $P$ as the exponent:
#' \deqn{distance(a,b) := \sum_i{(a_i - b_i)^P}^{1/P}
#' }{distance(a, b) := (\sum(a_i - b_i)^P)^(1/P)}
#' 
#' @param block either a vector to be used as the time series, or a 
#'   data.frame or matrix where each column is a time series
#' @param lib a 2-column matrix (or 2-element vector) where each row specifies 
#'   the first and last *rows* of the time series to use for attractor 
#'   reconstruction
#' @param pred (same format as lib), but specifying the sections of the time 
#'   series to forecast.
#' @param norm_type the distance function to use. see 'Details'
#' @param P the exponent for the P norm
#' @param E the embedding dimensions to use for time delay embedding
#' @param tau the lag to use for time delay embedding
#' @param tp the prediction horizon (how far ahead to forecast)
#' @param max_lag the maximum number of lags to use for variable combinations. 
#'   So if max_lag == 3, a variable X will appear with lags X[t], X[t - tau], 
#'   X[t - 2*tau]
#' @param num_neighbors the number of nearest neighbors to use for the 
#'   in-sample prediction (any of "e+1", "E+1", "e + 1", "E + 1" will peg this 
#'   parameter to E+1 for each run, any value < 1 will use all possible 
#'   neighbors.)
#' @param k the number of embeddings to use (any of "sqrt", "SQRT" will use 
#'   k = floor(sqrt(m)))
#' @param na.rm logical. Should missing values (including `NaN`` be omitted 
#'   from the calculations?)
#' @param target_column the index (or name) of the column to forecast
#' @param stats_only specify whether to output just the forecast statistics or 
#'   the raw predictions for each run
#' @param first_column_time indicates whether the first column of the given 
#'   block is a time column (and therefore excluded when indexing)
#' @param exclusion_radius excludes vectors from the search space of nearest 
#'   neighbors if their *time index* is within exclusion_radius (NULL turns 
#'   this option off)
#' @param silent prevents warning messages from being printed to the R console
#' @return A data.frame with components for the parameters and forecast 
#'   statistics:
#' \tabular{ll}{
#'   E \tab embedding dimension\cr
#'   tau \tab time lag\cr
#'   tp \tab prediction horizon\cr
#'   nn \tab number of neighbors\cr
#'   k \tab number of embeddings used\cr
#'   num_pred \tab number of predictions\cr
#'   rho \tab correlation coefficient between observations and predictions\cr
#'   mae \tab mean absolute error\cr
#'   rmse \tab root mean square error\cr
#'   perc \tab percent correct sign\cr
#'   p_val \tab p-value that rho is significantly greater than 0 using Fisher's 
#'   z-transformation\cr
#'   const_rho \tab same as rho, but for the constant predictor\cr
#'   const_mae \tab same as mae, but for the constant predictor\cr
#'   const_rmse \tab same as rmse, but for the constant predictor\cr
#'   const_perc \tab same as perc, but for the constant predictor\cr
#'   const_p_val \tab same as p_val, but for the constant predictor
#' }
#' If \code{stats_only == FALSE}, then additionally a list column:
#' \tabular{ll}{
#'   model_output \tab data.frame with columns for the time index, 
#'   observations, and predictions\cr
#' }
#' @examples 
#' data("block_3sp")
#' block <- block_3sp[, c(2, 5, 8)]
#' multiview(block, k = c(1, 3, "sqrt"))
#' 
multiview <- function(block, lib = c(1, floor(NROW(block) / 2)), 
                      pred = c(floor(NROW(block) / 2) + 1, NROW(block)), 
                      norm_type = c("L2 norm", "L1 norm", "P norm"), P = 0.5, 
                      E = 3, tau = 1, tp = 1, max_lag = 3, 
                      num_neighbors = "e+1", k = "sqrt", na.rm = FALSE, 
                      target_column = 1, 
                      stats_only = TRUE, first_column_time = FALSE, 
                      exclusion_radius = NULL, silent = FALSE)
{
    # setup params
    lib <- coerce_lib(lib)
    pred <- coerce_lib(pred)
    if (any(match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1")), 
           na.rm = TRUE))
        num_neighbors <- E + 1
    num_neighbors <- as.numeric(num_neighbors)
    
    # generate lagged block and list of embeddings
    if (max_lag < 1)
        warning("Maximum lag must be positive - setting to 1.")
    num_vars <- NCOL(block)
    if (first_column_time)
    {
        num_vars <- num_vars - 1
        lagged_block <- make_block(block[, 2:NCOL(block)], max_lag = max_lag, 
                                   t = block[, 1], lib = lib, tau = tau)
    } else {
        lagged_block <- make_block(block, max_lag = max_lag, 
                                   lib = lib, tau = tau)
    }
    
    embeddings_list <- t(combn(num_vars * max_lag, E, simplify = TRUE))
    valid_embeddings_idx <- apply(embeddings_list %% max_lag, 1, 
                                  function(x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings_idx, ]
    my_embeddings <- lapply(1:NROW(embeddings_list),
                            function(i) {embeddings_list[i, ]})
    
    ## make sure that if target_column is given as a column index, it
    ## is aligned with the lagged data frame.
    if (is.numeric(target_column))
        target_column <- 1 + max_lag * (target_column - 1)
    
    # make in-sample forecasts
    in_results <- block_lnlp(lagged_block, lib = lib, pred = lib, 
                             norm_type = norm_type, P = P, method = "simplex", 
                             tp = tp, num_neighbors = num_neighbors, 
                             columns = my_embeddings, 
                             target_column = target_column, 
                             stats_only = TRUE, first_column_time = TRUE, 
                             exclusion_radius = exclusion_radius, 
                             silent = silent)
    
    # rank embeddings
    num_embeddings <- NROW(in_results)
    k[k == "sqrt"] <- floor(sqrt(num_embeddings))
    k <- as.numeric(k)
    k_list <- sort(unique(pmin(k, num_embeddings)))
    
    in_sample_ranking <- order(in_results$rho, decreasing = TRUE)
    ranked_embeddings <- my_embeddings[in_sample_ranking]
    best_embeddings <- ranked_embeddings[1:max(k_list)]
    
    # make out-sample forecasts
    out_results <- block_lnlp(lagged_block, lib = lib, pred = pred, 
                              norm_type = norm_type, P = P, method = "simplex", 
                              tp = tp, num_neighbors = 1, 
                              columns = best_embeddings, 
                              target_column = target_column, 
                              stats_only = FALSE, first_column_time = TRUE, 
                              exclusion_radius = exclusion_radius, 
                              silent = silent)
    out_time <- out_results$model_output[[1]]$time
    out_obs <- out_results$model_output[[1]]$obs
    out_forecasts <- do.call(cbind, 
                             lapply(out_results$model_output, function(df) {
                                 df$pred})
    )
    
    # compute mve forecasts
    mve_forecasts <- lapply(k_list, function(k) {
        if (k > 1)
        {
            return(rowMeans(out_forecasts[, 1:k], na.rm = na.rm))
        } else {
            return(out_forecasts[, 1])
        }
    })
    
    params <- data.frame(E = E, tau = tau, tp = tp, 
                         nn = num_neighbors, k = k_list)
    
    # return output
    output <- lapply(mve_forecasts, function(pred) {
        if (stats_only)
        {
            df <- compute_stats(out_obs, pred)
        } else {
            df <- compute_stats(out_obs, pred)
            df$model_output <- I(list(data.frame(
                time = out_time,
                obs = out_obs, 
                pred = pred)))
        }
        return(df)
    })
    return(cbind(params, do.call(rbind, output), row.names = NULL))
}