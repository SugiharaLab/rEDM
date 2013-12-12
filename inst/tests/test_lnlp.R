data("two_species_model")
output <- simplex(two_species_model$x, lib = c(1, 500), pred = c(501, 1000), 
                  stats_only = TRUE)