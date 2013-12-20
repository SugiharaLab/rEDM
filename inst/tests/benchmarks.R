library(microbenchmark)

setwd("~/code/Redm/inst/tests/")
source("lnlp_func.R")

data("two_species_model")

write.table(two_species_model[, c("time", "x")], file = "lnlp_data.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

lib <- c(1, 1000)
pred <- c(1, 1000)

f1 <- function () {
    s_map(two_species_model$x, lib = lib, pred = pred, E = 3)
}
f2 <- function() {
    lnlp(lib_file = lib, pred_file = pred, pred_type = 1, E = 3, num_neighbors = 0, 
         theta_epsilon = c(0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 
                           0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6, 8))
}

microbenchmark(f1(), f2(), times = 20)