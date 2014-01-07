library(Redm)
library(microbenchmark)

setwd("~/code/Redm/inst/tests/")
source("lnlp_func.R")
source("block_lnlp_func.R")

data("two_species_model")

block <- cbind(two_species_model$time, 
               two_species_model$x, 
               c(NA, two_species_model$x[1:(NROW(two_species_model)-1)]))

lib <- c(1, 1000)
pred <- c(1, 1000)

f1 <- function () {
    block_lnlp(block, lib, pred, columns = c(1, 2), first_column_time = TRUE)
}

write.table(block, file = "block_lnlp_data.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

f2 <- function() {
    block_lnlp_f(lib_file = lib, pred_file = pred, columns = c(1, 2), dump = 1)
}