library(rEDM)

data("block_3sp")

block <- block_3sp[, c(2, 5, 8)]

# do multiview using top 1, top 3, top k = sqrt(m) embeddings
multiview(block, k = c(1, 3, "sqrt"), 
          silent = TRUE)

# same, but save output for plotting
output <- multiview(block, k = c(1, 3, "sqrt"), 
          stats_only = FALSE, 
          silent = TRUE)

plot(output$model_output[[3]]$obs, 
     output$model_output[[3]]$pred)
