context("Check structure of included datasets")

test_that("two_species_model is correct", {
    expect_error(data("two_species_model"), NA)
    expect_true(exists("two_species_model"))
    expect_equal(dim(two_species_model), c(1000, 3))
    expect_equal(names(two_species_model), c("time", "x", "y"))
})