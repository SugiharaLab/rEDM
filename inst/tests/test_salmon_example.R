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

preprocess_data <- function(stock_name = "Early Stuart")
{
    # load data
    data("FR_sockeye")
    
    # process biological data
    stock_df <- subset(sockeye_data, stk == stock_name)
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

stock_df <- preprocess_data()
forecasts <- make_forecasts(stock_df)

make_forecasts <- function(stock_df, columns = c("eff, D_jun"), 
                           lib = c(1, NROW(stock_df)), 
                           pred = c(1, NROW(stock_df)))
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
    env <- normalize(stock_df[,env_names])
    
    # make block
    block <- data.frame(years = years, eff = spawners, 
                        rec4 = recruits_4, rec5 = recruits_5)
    block <- cbind(block, env[valid, ])
    
    if(length(returns) < 2) # check for enough data
        return(data.frame(year = NaN, obs = NaN, pred = NaN))
    
    rec4_preds <- do.call(cbind, block_lnlp_4(block, target_column = 2, columns = columns))
    rec5_preds <- do.call(cbind, block_lnlp_4(block, target_column = 3, columns = columns))
    rec4_preds <- rec4_preds*sigma_4 + mu_4
    rec5_preds <- rec5_preds*sigma_5 + mu_5
    forecasts <- data.frame(rec4_preds + rbind(NA, rec5_preds[1:NROW(block)-1,]))
    names(forecasts) <- lapply(columns, function(v) paste(v, sep = "", collapse = ", "))
    output <- cbind(year = years, obs = returns, forecasts)
}











