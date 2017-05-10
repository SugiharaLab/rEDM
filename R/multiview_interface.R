#' Perform forecasting using multiview embedding
#'
#' \code{multiview} applies the method described in Ye & Sugihara (2016) for 
#' forecasting, wherein multiple attractor reconstructions are tested, and a 
#' single nearest neighbor is selected from each of the top "k" reconstructions 
#' to produce final forecasts.
#' 
#' uses multiple time series given as input to generate an 
#' attractor reconstruction, and then applies the simplex projection or s-map 
#' algorithm to make forecasts. This method generalizes the \code{simplex} and 
#' \code{s-map} routines, and allows for "mixed" embeddings, where multiple time 
#' series can be used as different dimensions of an attractor reconstruction.
#' 
#' The default parameters are set so that, given a matrix of time series, 
#' forecasts will be produced for the first column. By default, all possible 
#' combinations of the columns are used for the attractor construction, the 
#' k = sqrt(m) heuristic will be used, forecasts will be one time step ahead. 
#' Rownames will be converted to numeric if possible to be used as the time 
#' index, otherwise 1:NROW will be used instead. The default lib and pred are 
#' to use the first half of the data for the "library" and to predict over the 
#' second half of the data. Unless otherwise set, the output will be just the 
#' forecast statistics.
#' 
#' norm_type "L2 norm" (default) uses the typical Euclidean distance:
#' \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}}{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
#' norm_type "L1 norm" uses the Manhattan distance:
#' \deqn{distance(a,b) := \sum_i{|a_i - b_i|}}{distance(a, b) := \sum|a_i - b_i|}
#' norm type "P norm" uses the P norm, generalizing the L1 and L2 norm to use $p$ as the exponent:
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
#' @param max_lag the maximum number of lags to use for variable combinations
#' @param num_neighbors the number of nearest neighbors to use for the in-sample 
#'   prediction (any of "e+1", "E+1", "e + 1", "E + 1" will peg this parameter 
#'   to E+1 for each run, any value < 1 will use all possible neighbors.)
#' @param k the number of embeddings to use (any of "sqrt", "SQRT" will use 
#'   k = floor(sqrt(m)))
#' @param na.rm logical. Should missing values (including \code{NaN} be omitted 
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
#' @param short_output specifies whether to return a truncated output data.frame
#'   whose rows only include the predictions made and not the whole input block
#' @return If stats_only, then a data.frame with components for the parameters 
#'   and forecast statistics:
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
#' Otherwise, a list where the number of elements is equal to the number of runs 
#'   (unique parameter combinations). Each element is a list with the following 
#'   components:
#' \tabular{ll}{
#'   params \tab data.frame of parameters (E, tau, tp, nn, k)\cr
#'   lib_stats \tab data.frame of in-sample forecast statistics\cr
#'   model_output \tab data.frame with columns for the time index, observations, 
#'     and predictions\cr
#'   pred_stats \tab data.frame of forecast statistics\cr
#' }
#' @examples 
#' data("block_3sp")
#' block <- block_3sp[, c(2, 5, 8)]
#' multiview(block, k = c(1, 3, "sqrt"))
#' @export
multiview <- function(block, lib = c(1, floor(NROW(block)/2)), 
                      pred = c(floor(NROW(block)/2), NROW(block)), 
                      norm_type = c("L2 norm", "L1 norm", "P norm"), P = 0.5, 
                      E = 3, tau = 1, tp = 1, max_lag = 3, num_neighbors = "e+1", 
                      k = "sqrt", na.rm = FALSE, target_column = 1, 
                      stats_only = TRUE, first_column_time = FALSE, 
                      exclusion_radius = NULL, silent = FALSE, 
                      short_output = FALSE)
{
    # setup params
    if(is.vector(lib))
        lib <- matrix(lib, ncol = 2, byrow = TRUE)
    if(is.vector(pred))
        pred <- matrix(pred, ncol = 2, byrow = TRUE)
    if(any(match(num_neighbors, c("e+1", "E+1", "e + 1", "E + 1")), na.rm = TRUE))
        num_neighbors <- E + 1
    
    # generate lagged block and list of embeddings
    if(max_lag < 0)
        warning("Maximum lag must be non-negative - setting to 0.")
    num_vars <- NCOL(block)
    if(first_column_time)
    {
        num_vars <- num_vars - 1
        lagged_block <- make_block(block[, 2:NCOL(block)], max_lag = max_lag, 
                                   t = block[,1], lib = lib, tau = tau)
    } else {
        lagged_block <- make_block(block, max_lag = max_lag, 
                                   lib = lib, tau = tau)
    }
    
    embeddings_list <- t(combn(num_vars*max_lag, E, simplify = TRUE))
    valid_embeddings_idx <- apply(embeddings_list %% max_lag, 1, function (x) {1 %in% x})
    embeddings_list <- embeddings_list[valid_embeddings_idx,]
    my_embeddings <- lapply(1:NROW(embeddings_list), function(i) {embeddings_list[i,]})
    
    # make in-sample forecasts
    in_sample_results <- block_lnlp(lagged_block, lib = lib, pred = lib, 
                                    norm_type = norm_type, P = P, method = "simplex", 
                                    tp = tp, num_neighbors = num_neighbors, 
                                    columns = my_embeddings, target_column = target_column, 
                                    stats_only = TRUE, first_column_time = TRUE, 
                                    exclusion_radius = exclusion_radius, silent = silent)
    
    # rank embeddings
    num_embeddings <- NROW(in_sample_results)
    k[k == "sqrt"] <- floor(sqrt(num_embeddings))
    k <- as.numeric(k)
    k_list <- sort(unique(pmin(k, num_embeddings)))
    
    in_sample_ranking <- order(in_sample_results$rho, decreasing = TRUE)
    ranked_embeddings <- my_embeddings[in_sample_ranking]
    best_embeddings <- ranked_embeddings[1:max(k_list)]
    
    # make out-sample forecasts
    out_sample_results <- block_lnlp(lagged_block, lib = lib, pred = pred, 
                                     norm_type = norm_type, P = P, method = "simplex", 
                                     tp = tp, num_neighbors = 1, 
                                     columns = best_embeddings, target_column = target_column, 
                                     stats_only = FALSE, first_column_time = TRUE, 
                                     exclusion_radius = exclusion_radius, silent = silent, 
                                     short_output = short_output)
    out_sample_output <- out_sample_results[[1]]$model_output
    out_sample_forecasts <- do.call(cbind, lapply(out_sample_results, function(block_lnlp_output) {
        block_lnlp_output$model_output$pred}))
    
    # compute mve forecasts
    mve_forecasts <- lapply(k_list, function(k) {
        if(k > 1)
        {
            return(rowMeans(out_sample_forecasts[, 1:k], na.rm = na.rm))
        } else {
            return(mve_pred <- out_sample_forecasts[, 1])
        }
    })
    
    mve_stats <- do.call(rbind, lapply(mve_forecasts, function(pred) {
        compute_stats(out_sample_output$obs, pred)
    }))
    params <- data.frame(E = E, tau = tau, tp = tp, nn = num_neighbors, k = k_list)
    
    # return output
    if(stats_only)
    {
        output <- cbind(params, mve_stats)
    } else {
        output <- lapply(seq_along(k_list), function(i) {
            return(list(params = params[i,], 
                        lib_stats = in_sample_results[in_sample_ranking[1:k_list[i]],], 
                        model_output = data.frame(time = out_sample_output$time, 
                                                  obs = out_sample_output$obs, 
                                                  pred = mve_forecasts[[i]]), 
                        pred_stats = mve_stats[i, ]))
        })
    }
    return(output)
}

#' Make a lagged block for multiview
#'
#' \code{make_block} generates a lagged block with the appropriate max_lag and 
#' tau, while respecting lib (by inserting NANs, when trying to lag past lib 
#' regions)
#' 
#' @param block a data.frame or matrix where each column is a time series
#' @param max_lag the total number of lags to include for each variable
#' @param t the time index for the block
#' @param lib a 2-column matrix (or 2-element vector) where each row specifes the 
#'   first and last *rows* of the time series to use for attractor reconstruction
#' @param tau the lag to use for time delay embedding
#' @return A data.frame with the lagged columns and a time column
#' 
make_block <- function(block, max_lag = 3, t = NULL, lib = NULL, tau = 1)
{
    num_vars <- NCOL(block)
    num_rows <- NROW(block)
    output <- matrix(NA, nrow = num_rows, ncol = 1+num_vars*max_lag)
    if(is.null(t))
        output[, 1] <- 1:num_rows
    else
        output[, 1] <- t
    
    # add max_lag lags for each column in block
    col_index <- 2
    for (j in 1:num_vars)
    {
        ts <- block[,j]
        output[, col_index] <- ts
        col_index <- col_index + 1
        
        if(max_lag > 1)
        {
            for(i in 1:(max_lag-1))
            {
                ts <- c(rep_len(NA, tau), ts[1:(num_rows-tau)])
                if(!is.null(lib))
                {
                    for(k in 1:NROW(lib))
                    {
                        ts[lib[k,1]] <- NA
                    }
                }
                output[, col_index] <- ts
                col_index <- col_index + 1
            }
        }
    }
    
    return(output)
}