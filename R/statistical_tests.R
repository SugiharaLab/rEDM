#' Randomization test for nonlinearity using S-maps and surrogate data
#' 
#' \code{\link{test_nonlinearity}} tests for nonlinearity using S-maps by 
#' comparing improvements in forecast skill (delta rho and delta mae) between 
#' linear and nonlinear models with a null distribution from surrogate data.
#' 
#' @param ts the original time series
#' @param method which algorithm to use to generate surrogate data
#' @param num_surr the number of null surrogates to generate
#' @param T_period the period of seasonality for seasonal surrogates 
#'   (ignored for other methods)
#' @param E the embedding dimension for s_map
#' @param ... optional arguments to s_map
#' @return A data.frame containing the following components:
#' \tabular{ll}{
#'   delta_rho \tab the value of the delta rho statistic\cr
#'   delta_mae \tab the value of the delat mae statistic\cr
#'   num_surr \tab the size of the null distribution\cr
#'   delta_rho_p_value \tab the p-value for delta rho\cr
#'   delta_mae_p_value \tab the p-value for delta mae\cr
#' }
#' 
test_nonlinearity <- function(ts, method = "ebisuzaki", num_surr = 200, 
                              T_period = 1, E = 1, ...)
{
    compute_stats <- function(ts, ...)
    {
        results <- s_map(ts, stats_only = TRUE, silent = TRUE, ...)
        delta_rho <- max(results$rho) - results$rho[results$theta == 0]
        delta_mae <- results$mae[results$theta == 0] - min(results$mae)
        return(c(delta_rho = delta_rho, delta_mae = delta_mae))
    }
    
    actual_stats <- compute_stats(ts, ...)
    delta_rho <- actual_stats["delta_rho"]
    delta_mae <- actual_stats["delta_mae"]
    names(delta_rho) <- NULL
    names(delta_mae) <- NULL
    surrogate_data <- make_surrogate_data(ts, method, num_surr, T_period)
    null_stats <- data.frame(t(apply(surrogate_data, 2, compute_stats, ...)))
    
    return(data.frame(
        delta_rho = delta_rho, 
        delta_mae = delta_mae, 
        num_surr = num_surr, 
        E = E, 
        delta_rho_p_value = (sum(null_stats$delta_rho > delta_rho) + 1) / 
            (num_surr + 1), 
        delta_mae_p_value = (sum(null_stats$delta_mae > delta_mae) + 1) / 
            (num_surr + 1)))
}