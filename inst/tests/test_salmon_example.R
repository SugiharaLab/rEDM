# scale each cycle line to mean = 0, variance = 1
normalize_by_cycle_line <- function(ts)
{
    n <- length(ts)
    means <- rep.int(NA, times = 4)
    sds <- rep.int(NA, times = 4)
    mu <- rep.int(NA, times = n)
    sigma <- rep.int(NA, times = n)
    for(k in 1:4)
    {
        index <- seq(from = k, to = n, by = 4)
        means[k] <- mean(ts[index], na.rm = TRUE)
        sds[k] <- sd(ts[index], na.rm = TRUE)
        mu[index] <- means[k]
        sigma[index] <- sds[k]
    }
    ts <- (ts - mu) / sigma
    df <- data.frame(cbind(ts, mu, sigma))
    return(df)
}

normalize <- function(block)
{
    if(NCOL(block) > 1)
    {
        n <- NROW(block)
        means <- sapply(block, mean, na.rm = TRUE)
        sds <- sapply(block, sd, na.rm = TRUE)
        return((block - matrix(rep(means, each = n), nrow = n)) / 
                   matrix(rep(sds, each = n), nrow = n))
    }
    else
        return((block - mean(block, na.rm = TRUE)) / sd(block, na.rm = TRUE))
}

preprocess_data <- function(stock_name = "Early Stuart")
{
    # load data
    data("FR_sockeye")
    
    # process biological data
    stock_df <- subset(sockeye_data, stk == stock_name)
    n <- NROW(stock_df)
    stock_df$rec45 <- stock_df$rec4 + stock_df$rec5
    stock_df$ret <- stock_df$rec4 + c(NA, stock_df$rec5[1:(n-1)]) # age-4 and age-5 fish
    
    temp <- normalize_by_cycle_line(stock_df$rec4)
    stock_df$rec4_n <- temp$ts
    stock_df$rec4_mu <- temp$mu
    stock_df$rec4_sigma <- temp$sigma
    
    temp <- normalize_by_cycle_line(stock_df$rec5)
    stock_df$rec5_n <- temp$ts
    stock_df$rec5_mu <- temp$mu
    stock_df$rec5_sigma <- temp$sigma
    
    temp <- normalize_by_cycle_line(stock_df$eff)
    stock_df$eff_n <- temp$ts
    stock_df$eff_mu <- temp$mu
    stock_df$eff_sigma <- temp$sigma
    
    # include environmental data
    desired_years <- stock_df$yr + 2
    index_in_env_data <- match(desired_years, environmental_data$year)
    index_in_stock_df <- 1:length(desired_years)
    discharge_names <- c("D_max", "D_apr", "D_may", "D_jun")
    discharge_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = length(discharge_names)))
    discharge_cols[index_in_stock_df,] <- environmental_data[index_in_env_data, discharge_names]
    stock_df[, discharge_names] <- discharge_cols
    
    temperature_names <- c("ET_apr", "ET_may", "ET_jun", "PT_apr", "PT_may", "PT_jun", "PT_jul")
    temperature_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = length(temperature_names)))
    temperature_cols[index_in_stock_df,] <- environmental_data[index_in_env_data, temperature_names]
    stock_df[, temperature_names] <- temperature_cols
    
    desired_years <- stock_df$yr + 1
    index_in_env_data <- match(desired_years, environmental_data$year)
    pdo_names <- "PDO_win"
    pdo_cols <- data.frame(matrix(NA, nrow = length(desired_years), ncol = length(pdo_names)))
    pdo_cols[index_in_stock_df,] <- environmental_data[index_in_env_data, pdo_names]
    stock_df[, pdo_names] <- pdo_cols
    
    return(stock_df)
}

make_forecasts <- function(stock_df, columns = c("eff", "D_jun"))
{
    valid <- is.finite(stock_df$rec45) & is.finite(stock_df$eff)
    years <- stock_df$yr[valid]
    returns <- stock_df$ret[valid]
    spawners <- stock_df$eff_n[valid]
    recruits_4 <- stock_df$rec4_n[valid]
    mu_4 <- stock_df$rec4_mu[valid]
    sigma_4 <- stock_df$rec4_sigma[valid]
    recruits_5 <- stock_df$rec5_n[valid]
    mu_5 <- stock_df$rec5_mu[valid]
    sigma_5 <- stock_df$rec5_sigma[valid]
    env_names = c("D_max", "D_apr", "D_may", "D_jun", 
                  "ET_apr", "ET_may", "ET_jun", "PT_apr", "PT_may", "PT_jun", "PT_jul", 
                  "PDO_win")
    env <- normalize(stock_df[,env_names])
    
    # make block
    block <- data.frame(years = years, eff = spawners, 
                        rec4 = recruits_4, rec5 = recruits_5)
    block <- cbind(block, env[valid, ])
    
    # make lib and pred
    lib <- c(1, NROW(block)) 
    pred <- c(1, NROW(block))
    
    rec4_preds <- block_lnlp(block, lib = lib, pred = pred, 
                             target_column = 2, columns = columns, 
                             stats_only = FALSE)[[1]]$model_output
    rec5_preds <- block_lnlp(block, lib = lib, pred = pred, 
                             target_column = 3, columns = columns, 
                             stats_only = FALSE)[[1]]$model_output
    rec4 <- rec4_preds$pred*sigma_4 + mu_4
    rec5 <- rec5_preds$pred*sigma_5 + mu_5
    rec4_var <- rec4_preds$pred_var * sigma_4 * sigma_4
    rec5_var <- rec5_preds$pred_var * sigma_5 * sigma_5
    forecasts <- data.frame(pred = rec4 + c(NA, rec5[1:NROW(block)-1]), 
                            var = rec4_var + c(NA, rec5_var[1:NROW(block)-1]))
    forecasts$std_err <- sqrt(forecasts$var)
    output <- cbind(year = years, obs = returns, forecasts)
}

stock_df <- preprocess_data()
forecasts <- make_forecasts(stock_df)













