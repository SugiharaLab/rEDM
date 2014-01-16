library(Redm)
library(microbenchmark)

setwd("~/code/Redm/inst/tests/")
source("lnlp_func.R")
source("block_lnlp_func.R")

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

smap_new <- function() {
    s_map(two_species_model$x, lib = c(1, 500), pred = c(501, 1000), E = 3)
}
smap_old <- function() {
    lnlp(lib_file = c(1, 500), pred_file = c(501, 1000), E = 3, pred_type = 1, num_neighbors = 0, 
         theta_epsilon = c(0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 
                           0.3, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6, 8))
}

block_new <- function () {
    block_lnlp(block, columns = list(c(1, 2), c(1)), first_column_time = TRUE)
}
block_old <- function() {
    block_lnlp_f(columns = c(1, 2), dump = 1)
}

# benchmarks
if(FALSE)
{
    microbenchmark(simplex_new(), simplex_old(), times = 20)
    microbenchmark(smap_new(), smap_old(), times = 20)
    microbenchmark(block_new(), block_old())
}