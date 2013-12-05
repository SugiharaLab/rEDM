compute_stats <- function(out)
{
    num_pred <- sum(is.finite(out$obs) & is.finite(out$pred))
    rho <- cor(out$obs, out$pred, use = "pairwise.complete.obs")
    mae <- mean(abs(out$pred - out$obs), na.rm = TRUE)
    rmse <- sqrt(mean((out$pred - out$obs)^2, na.rm = TRUE))
    return(data.frame(num_pred = num_pred, rho = rho, mae = mae, rmse = rmse))
}
