#' Perform forecasting using multiview embedding
#'
#' \code{\link{multiview}} applies the method described in Ye & Sugihara (2016) for 
#'   forecasting, wherein multiple attractor reconstructions are tested, and a 
#'   single nearest neighbor is selected from each of the top \code{k} 
#'   reconstructions to produce final forecasts.
#' 
#'   uses multiple time series given as input to generate an attractor 
#'   reconstruction, and then applies the simplex projection or s-map algorithm 
#'   to make forecasts. This method generalizes the \code{\link{simplex}} and 
#'   \code{\link{s_map}} routines, and allows for "mixed" embeddings, where multiple 
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
#' @param max_lag the maximum number of lags to use for variable combinations. 
#'   So if max_lag == 3, a variable X will appear with lags X[t], X[t - tau], 
#'   X[t - 2*tau]
#' @param k the number of embeddings to use ("sqrt" will use k = floor(sqrt(m)), 
#'   "all" or values less than 1 will use k = m)
#' @param na.rm logical. Should missing values (including `NaN`` be omitted 
#'   from the calculations?)
#' @param save_lagged_block specify whether to output the lagged block that 
#'   is constructed as part of running \code{multiview}
#' @return A data.frame with components for the parameters and forecast 
#'   statistics:
#' \tabular{ll}{
#'   E \tab embedding dimension\cr
#'   tau \tab time lag\cr
#'   tp \tab prediction horizon\cr
#'   nn \tab number of neighbors\cr
#'   k \tab number of embeddings used\cr
#' }
#' \tabular{ll}{
#'   \code{E} \tab embedding dimension\cr
#'   \code{tau} \tab time lag\cr
#'   \code{tp} \tab prediction horizon\cr
#'   \code{nn} \tab number of neighbors\cr
#'   \code{k} \tab number of embeddings used\cr
#'   \code{num_pred} \tab number of predictions\cr
#'   \code{rho} \tab correlation coefficient between observations and 
#'     predictions\cr
#'   \code{mae} \tab mean absolute error\cr
#'   \code{rmse} \tab root mean square error\cr
#'   \code{perc} \tab percent correct sign\cr
#'   \code{p_val} \tab p-value that rho is significantly greater than 0 using 
#'     Fisher's z-transformation\cr
#'   \code{model_output} \tab data.frame with columns for the time index, 
#'     observations, predictions, and estimated prediction variance
#'     (if \code{stats_only == FALSE})\cr
#'   \code{embeddings} \tab list of the columns used in each of the embeddings 
#'     that comprise the model (if \code{stats_only == FALSE})\cr
#' }
#' @examples 
#' data("block_3sp")
#' block <- block_3sp[, c(2, 5, 8)]
#' multiview(block, k = c(1, 3, "sqrt"))
#' 
multiview <- function(block, lib = c(1, floor(NROW(block) / 2)), 
                      pred = c(floor(NROW(block) / 2) + 1, NROW(block)), 
                      norm = 2, E = 3, tau = 1, tp = 1, max_lag = 3, 
                      num_neighbors = "e+1", k = "sqrt", na.rm = FALSE, 
                      target_column = 1, 
                      stats_only = TRUE, save_lagged_block = FALSE, 
                      first_column_time = FALSE, 
                      exclusion_radius = NULL, silent = FALSE)
{
    # setup params
    lib <- coerce_lib(lib, silent = silent)
    pred <- coerce_lib(pred, silent = silent)
    if (any(match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1")), 
           na.rm = TRUE))
        num_neighbors <- E + 1
    num_neighbors <- as.numeric(num_neighbors)
    
    # generate lagged block and list of embeddings
    if (max_lag < 1)
        rEDM_warning("Maximum lag must be positive - setting to 1.", 
                     silent = silent)
    num_vars <- NCOL(block)
    if (first_column_time)
    {
        num_vars <- num_vars - 1
        lagged_block <- make_block(block[, 2:NCOL(block)], max_lag = max_lag, 
                                   t = block[, 1], lib = lib, tau = tau, 
                                   restrict_to_lib = FALSE)
    } else {
        lagged_block <- make_block(block, max_lag = max_lag, 
                                   lib = lib, tau = tau, 
                                   restrict_to_lib = FALSE)
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
                             norm = norm, method = "simplex", 
                             tp = tp, num_neighbors = num_neighbors, 
                             columns = my_embeddings, 
                             target_column = target_column, 
                             stats_only = TRUE, first_column_time = TRUE, 
                             exclusion_radius = exclusion_radius, 
                             silent = silent)
    
    # rank embeddings
    num_embeddings <- NROW(in_results)
    k[tolower(k) == "sqrt"] <- floor(sqrt(num_embeddings))
    k[k < 1] <- num_embeddings
    k[tolower(k) == "all"] <- num_embeddings
    k <- as.numeric(k)
    k_list <- sort(unique(pmin(k, num_embeddings)))
    
    in_sample_ranking <- order(in_results$rho, decreasing = TRUE)
    ranked_embeddings <- my_embeddings[in_sample_ranking]
    best_embeddings <- ranked_embeddings[1:max(k_list)]
    
    # make out-sample forecasts
    out_results <- block_lnlp(lagged_block, lib = lib, pred = pred, 
                              norm = norm, method = "simplex", 
                              tp = tp, num_neighbors = 1, 
                              columns = best_embeddings, 
                              target_column = target_column, 
                              stats_only = FALSE, first_column_time = TRUE, 
                              exclusion_radius = exclusion_radius, 
                              silent = silent)
    out_time <- out_results$model_output[[1]]$time
    out_obs <- out_results$model_output[[1]]$obs
    out_pred <- do.call(cbind, lapply(out_results$model_output, 
                                      function(df) {df$pred}))
    out_pred_var <- do.call(cbind, lapply(out_results$model_output, 
                                          function(df) {df$pred_var}))
    # compute mve forecasts
    mve_forecasts <- lapply(k_list, function(k) {
        if (k > 1) {
            return(data.frame(pred = rowMeans(out_pred[, 1:k], na.rm = na.rm), 
                              pred_var = rowMeans(out_pred_var[, 1:k], na.rm = na.rm) + 
                                  apply(out_pred_var[, 1:k], 1, var, na.rm = na.rm)))
        } else {
            return(data.frame(pred = out_pred[, 1], 
                              pred_var = out_pred_var[, 1]))
        }
    })

    # return output
    output <- lapply(mve_forecasts, function(mve_output) {
        if (stats_only)
        {
            df <- compute_stats(out_obs, mve_output$pred)
        } else {
            df <- compute_stats(out_obs, mve_output$pred)
            df$model_output <- I(list(data.frame(
                time = out_time,
                obs = out_obs, 
                pred = mve_output$pred, 
                pred_var = mve_output$pred_var)))
        }
        if (save_lagged_block)
        {
            df$lagged_block <- I(list(lagged_block))
        }
        return(df)
    })
    output <- do.call(rbind, output)
    
    if (!stats_only)
    {
        output$embeddings <- lapply(k_list, function(k) {
            best_embeddings[seq(k)]})
    }
    
    # setup params to append
    params <- data.frame(E = E, tau = tau, tp = tp, 
                         nn = num_neighbors, k = k_list)
    return(cbind(params, output, row.names = NULL))
}