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
    simplex_out <- data.frame(lapply(simplex_out, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(simplex_out) <- attributes(simplex_out)[sort(names(attributes(simplex_out)))]
    expect_equal(digest::digest(simplex_out), 
                 "19b7ca9c40138adf04ba91e6be68fa82")
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
    expect_equal(dim(model_output), c(200, 4))
    expect_equal(digest::digest(round(model_output, 4)), 
                 "3ca91650810e31a0f84eccc8d35a513c")
})

test_that("simplex works on time series", {
    expect_warning(output <- simplex(AirPassengers, 
                                     E = 7, stats_only = FALSE))
    model_output <- round(output$model_output[[1]], 4)
    expect_equal(digest::digest(model_output), "22d06a4d868f1e3070aef99f270f1dee")
    
    output <- output[, !(names(output) %in% "model_output")]
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "e3d42ba41f5d9e7527af7e31a2cad576")
})

test_that("simplex works on multivariate time series", {
    expect_warning(output <- simplex(EuStockMarkets, 
                                     E = 6, stats_only = FALSE))
    model_output <- round(output$model_output[[1]], 4)
    expect_equal(digest::digest(model_output), "a5bb66a31ce1816b11599acb59f3b55b")
    
    output <- output[, !(names(output) %in% "model_output")]
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "53a3720878b52780ddf89b9a604d49b2")
})

test_that("simplex error checking works", {
    expect_warning(simplex(1:10))
    expect_error(simplex(1:5, E = 5, silent = TRUE))
    expect_error(simplex(1:5, E = 2, tau = 4, silent = TRUE))
    expect_error(simplex(1:5, E = 1, tp = 5, silent = TRUE))
    expect_error(simplex(1:5, E = 1, tp = -5, silent = TRUE))
})
