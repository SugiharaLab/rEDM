library(rEDM)

data("two_species_model")

block <- two_species_model[1:200,]

# multivariate simplex projection using x and y to predict x
output <- block_lnlp(block, columns = c("x", "y"), 
                     first_column_time = TRUE, stats_only = FALSE)

output <- output[[1]]$model_output
output$pred_err <- sqrt(output$pred_var)
t <- 1:50
plot(t, output$obs[t], type = "l")
points(t, output$pred[t], col = "blue")
for(i in t)
{
    lines(c(i,i), c(output$pred[i]-output$pred_err[i], 
                    output$pred[i]+output$pred_err[i]), col = "blue")
}

# cross mapping using x to predict y
block_lnlp(block, target_column = 2, columns = c("x"), first_column_time = TRUE)
