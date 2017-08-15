library(rEDM)

data("two_species_model")

ts <- two_species_model$x[1:200]

# univariate forecasting with Gaussian processes and E = 1:10, using maximum 
#　likelihood to estimate params over library points
out <- tde_gp(ts)

# univariate forecasting with Gaussian processes and E = 5, and first half to 
# predict second half
out <- tde_gp(ts, lib = c(1, 100), pred = c(101, 200), E = 5, stats_only = FALSE)