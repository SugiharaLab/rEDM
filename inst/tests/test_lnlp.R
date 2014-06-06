data("two_species_model")

ts <- two_species_model$x[1:200]

# univariate simplex projection using E = 1:10, and leave-one-out cross-validation
x <- simplex(ts, stats_only = FALSE)

# univariate simplex projection using E = 1:10, and first half to predict second half
simplex(ts, lib = c(1, 100), pred = c(101, 200))

# univariate s-map using E = 2
s_map(ts, E = 2)
