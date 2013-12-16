library(microbenchmark)

setwd("~/code/Redm/inst/tests/")
source("lnlp_func.R")

data("two_species_model")

write.table(two_species_model[, c("time", "x")], file = "lnlp_data.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

lib <- c(1, 1000)
pred <- c(1, 1000)

out_Redm <- simplex(two_species_model$x, lib = lib, pred = pred)
out_bin <- lnlp(lib_file = lib, pred_file = pred)
