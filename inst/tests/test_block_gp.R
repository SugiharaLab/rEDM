library(rEDM)

data("two_species_model")

block <- two_species_model[1:200,]

# multivariate simplex projection using x and y to predict x
output <- block_gp(block, columns = c("x", "y"), 
                   first_column_time = TRUE, stats_only = FALSE, 
                   save_covariance_matrix = TRUE)
