num_replicates <- 100

data("two_species_model")

block <- cbind(two_species_model$time, 
               two_species_model$x, 
               c(NA, two_species_model$x[1:(NROW(two_species_model)-1)]))
write.table(block, file = "block_lnlp_data.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cbind(two_species_model$time, two_species_model$x), file = "lnlp_data.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

lib <- c(1, 1000)
pred <- c(1, 1000)

simplex_new <- function() {
    simplex(two_species_model$x)
}
simplex_old <- function() {
    lnlp()
}

a <- simplex_new()
b <- simplex_old()
expect_equal(a$rho, b$rho, tolerance = 0.0001, scale = 1)
expect_equal(a$mae, b$MAE, tolerance = 0.0001, scale = 1)
expect_equal(a$rmse, b$RMSE, tolerance = 0.0001, scale = 1)            

test_that("old simplex is the same as new simplex", {
    for(i in seq_len(num_replicates))
    {
        b <- simplex_new()
        expect_equal(a$rho, b$rho)
        expect_equal(a$mae, b$mae)
        expect_equal(a$rmse, b$rmse)            
    }
})