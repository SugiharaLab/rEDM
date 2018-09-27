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
    expect_equal(digest::digest(round(smap_out$rho, 4)), 
                 "51c159a4ab5d37fe1b0cf3b1c07bc509")
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
    expect_equal(dim(model_output), c(200, 4))
    expect_equal(digest::digest(round(model_output, 4)), 
                 "d180a19cc4629e64c36712153226e8e0")
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
    expect_equal(dim(smap_coefficients), c(200, 3))
    expect_equal(digest::digest(round(smap_coefficients, 4)), 
                 "d63e5b1c25e79b2c65352fa3ec118e99")
})

test_that("s-map smap_coefficient_covariances works", {
    expect_warning(smap_out <- s_map(ts, E = 2, theta = 1, 
                                     save_smap_coefficients = TRUE, 
                                     silent = TRUE), NA)
    expect_s3_class(smap_out, "data.frame")
    expect_true("smap_coefficient_covariances" %in% names(smap_out))
    expect_true(is.list(smap_out$smap_coefficient_covariances))
    expect_error(smap_coeff_covariances <- smap_out$smap_coefficient_covariances[[1]], NA)
    expect_true(is.list(smap_coeff_covariances))
    expect_equal(length(smap_coeff_covariances), 200)
    expect_null(smap_coeff_covariances[[1]])
    expect_null(smap_coeff_covariances[[200]])
    expect_equal(sapply(smap_coeff_covariances[2:199], dim), 
                 matrix(3, nrow = 2, ncol = 198))
    expect_error(covariance_mat <- do.call(rbind, smap_coeff_covariances[2:199]), NA)
    expect_equal(digest::digest(round(covariance_mat, 4)), 
                 "4c4643578927d95cc473c9bf873afcf6")
})

test_that("s-map works on time series", {
    expect_warning(output <- s_map(AirPassengers, 
                                   E = 7, theta = 1, stats_only = FALSE))
    model_output <- round(output$model_output[[1]], 4)
    expect_equal(digest::digest(model_output), "2e8bddf61e78cef8493b47046c00d071")
    
    output <- output[, !(names(output) %in% "model_output")]
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "b31f94cc3d234f45e4b26aa4a4444159")
})


test_that("s-map works on multivariate time series", {
    expect_warning(output <- s_map(EuStockMarkets, 
                                     E = 6, theta = 1, stats_only = FALSE))
    model_output <- round(output$model_output[[1]], 4)
    expect_equal(digest::digest(model_output), "b72187b45cc3bc56fae75daeec740e11")
    
    output <- output[, !(names(output) %in% "model_output")]
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "db27e5eeec3de0f898ee82f9837558f4")
})

test_that("s-map error checking works", {
    expect_warning(s_map(1:10))
    expect_error(s_map(1:5, E = 5, silent = TRUE))
    expect_error(s_map(1:5, E = 2, tau = 4, silent = TRUE))
    expect_error(s_map(1:5, E = 1, tp = 5, silent = TRUE))
    expect_error(s_map(1:5, E = 1, tp = -5, silent = TRUE))
})