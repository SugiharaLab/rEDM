context("Check s-map function")

data("two_species_model")
ts <- two_species_model$x[1:200]
theta_list <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 
                0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)

test_that("s-map works", {
    expect_warning(smap_out <- s_map(ts, E = 2, theta = theta_list))
    expect_s3_class(smap_out, "data.frame")
    expect_true("E" %in% names(smap_out))
    expect_true("tau" %in% names(smap_out))
    expect_true("tp" %in% names(smap_out))
    expect_true("nn" %in% names(smap_out))
    expect_true("theta" %in% names(smap_out))
    expect_true("num_pred" %in% names(smap_out))
    expect_true("rho" %in% names(smap_out))
    expect_true("mae" %in% names(smap_out))
    expect_true("rmse" %in% names(smap_out))
    expect_equal(NROW(smap_out), length(theta_list))
})

test_that("s-map model_output works", {
    expect_warning(smap_out <- s_map(ts, E = 2, theta = 1, 
                                     stats_only = FALSE, 
                                     silent = TRUE), NA)
    expect_s3_class(smap_out, "data.frame")
    expect_true("model_output" %in% names(smap_out))
    expect_true(is.list(smap_out$model_output))
    expect_error(model_output <- smap_out$model_output[[1]], NA)
    expect_s3_class(model_output, "data.frame")
    expect_true("time" %in% names(model_output))
    expect_true("obs" %in% names(model_output))
    expect_true("pred" %in% names(model_output))
    expect_true("pred_var" %in% names(model_output))
    expect_equal(dim(model_output), c(198, 4))
})

test_that("s-map smap_coefficients works", {
    expect_warning(smap_out <- s_map(ts, E = 2, theta = 1, 
                                     save_smap_coefficients = TRUE, 
                                     silent = TRUE), NA)
    expect_s3_class(smap_out, "data.frame")
    expect_true("smap_coefficients" %in% names(smap_out))
    expect_true(is.list(smap_out$smap_coefficients))
    expect_error(smap_coefficients <- smap_out$smap_coefficients[[1]], NA)
    expect_s3_class(smap_coefficients, "data.frame")
    expect_true("c_1" %in% names(smap_coefficients))
    expect_true("c_2" %in% names(smap_coefficients))
    expect_true("c_0" %in% names(smap_coefficients))
    expect_equal(dim(smap_coefficients), c(198, 3))
})

test_that("s-map error checking works", {
    expect_warning(s_map(1:10))
    expect_error(suppressWarnings(s_map(1:5, E = 5, silent = TRUE)))
    expect_error(suppressWarnings(s_map(1:5, E = 2, tau = 4, silent = TRUE)))
    expect_error(suppressWarnings(s_map(1:5, E = 1, tp = 5, silent = TRUE)))
    expect_error(suppressWarnings(s_map(1:5, E = 1, tp = -5, silent = TRUE)))
})