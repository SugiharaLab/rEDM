context("Check simplex function")

data("two_species_model")
ts <- two_species_model$x[1:200]

test_that("simplex works", {
    expect_error(simplex_out <- simplex(ts, lib = c(1, 100), 
                                        pred = c(101, 200)), 
                 NA)
    expect_s3_class(simplex_out, "data.frame")
    expect_true("E" %in% names(simplex_out))
    expect_true("tau" %in% names(simplex_out))
    expect_true("tp" %in% names(simplex_out))
    expect_true("nn" %in% names(simplex_out))
    expect_true("num_pred" %in% names(simplex_out))
    expect_true("rho" %in% names(simplex_out))
    expect_true("mae" %in% names(simplex_out))
    expect_true("rmse" %in% names(simplex_out))
    expect_equal(NROW(simplex_out), 10)
})

test_that("simplex model_output works", {
    expect_warning(simplex_out <- simplex(ts, E = 3, stats_only = FALSE))
    expect_s3_class(simplex_out, "data.frame")
    expect_true("model_output" %in% names(simplex_out))
    expect_true(is.list(simplex_out$model_output))
    expect_error(model_output <- simplex_out$model_output[[1]], NA)
    expect_s3_class(model_output, "data.frame")
    expect_true("time" %in% names(model_output))
    expect_true("obs" %in% names(model_output))
    expect_true("pred" %in% names(model_output))
    expect_true("pred_var" %in% names(model_output))
    expect_equal(dim(model_output), c(197, 4))
})

test_that("simplex error checking works", {
    expect_warning(simplex(1:10))
    expect_error(suppressWarnings(simplex(1:5, E = 5, silent = TRUE)))
    expect_error(suppressWarnings(simplex(1:5, E = 2, tau = 4, silent = TRUE)))
    expect_error(suppressWarnings(simplex(1:5, E = 1, tp = 5, silent = TRUE)))
    expect_error(suppressWarnings(simplex(1:5, E = 1, tp = -5, silent = TRUE)))
})
