data("two_species_model")

block <- two_species_model[1:200,]

# multivariate simplex projection using x and y to predict x
block_lnlp(block, columns = c("x", "y"), first_column_time = TRUE)

# cross mapping using x to predict y
block_lnlp(block, target_column = 2, columns = c("x"), first_column_time = TRUE)
