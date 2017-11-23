library(rEDM)

data("two_species_model")

block <- two_species_model[1:200, ]

# try different param combinations (no fitting)
output <- block_gp(block, columns = c("x", "y"), 
                   phi = 0.5402931, 
                   v_e = seq(from = -300, to = 300, by = 50), 
                   eta = 7, 
                   fit_params = FALSE, 
                   first_column_time = TRUE, stats_only = FALSE)
plot(output$model_output[[4]]$obs, 
     output$model_output[[4]]$pred)

# multivariate simplex projection using x and y to predict x
output <- block_gp(block, columns = c("x", "y"), 
                   first_column_time = TRUE, stats_only = FALSE, 
                   save_covariance_matrix = TRUE)
