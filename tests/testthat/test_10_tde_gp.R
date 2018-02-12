context("Check tde_gp function")

data("two_species_model")
ts <- two_species_model[1:200, "x"]

test_that("tde_gp works", {
    expect_error(output <- tde_gp(ts), #, silent = TRUE),
                 NA)
    expect_s3_class(output, "data.frame")
    expect_true("E" %in% names(output))
    expect_true("tau" %in% names(output))
    expect_true("embedding" %in% names(output))
    expect_true("tp" %in% names(output))
    expect_true("phi" %in% names(output))
    expect_true("v_e" %in% names(output))
    expect_true("eta" %in% names(output))
    expect_true("fit_params" %in% names(output))
    expect_true("num_pred" %in% names(output))
    expect_true("rho" %in% names(output))
    expect_true("mae" %in% names(output))
    expect_true("rmse" %in% names(output))
    expect_equal(NROW(output), 10)
    expect_equal(digest::digest(round(output$rho, 4)),
                 "effa17e5fc295299b049aa59b9d0a1bf")
})

test_that("tde_gp model_output works", {
    expect_error(output <- tde_gp(ts, lib = c(1, 100), pred = c(101, 200), 
                                  E = 5, stats_only = FALSE),
                 NA)
    expect_s3_class(output, "data.frame")
    expect_true("model_output" %in% names(output))
    expect_true(is.list(output$model_output))
    expect_error(model_output <- output$model_output[[1]], NA)
    expect_s3_class(model_output, "data.frame")
    expect_true("time" %in% names(model_output))
    expect_true("obs" %in% names(model_output))
    expect_true("pred" %in% names(model_output))
    expect_true("pred_var" %in% names(model_output))
    expect_equal(dim(model_output), c(99, 4))
    expect_equal(digest::digest(round(model_output, 4)),
                 "8f765da4f14bb15fd30981aff41c221f")
})