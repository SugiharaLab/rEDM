context("Check multiview function")

data("block_3sp")
block <- block_3sp[, c(2, 5, 8)]

test_that("multiview works", {
    expect_error(output <- multiview(block, k = c(1, 3, "sqrt"),
                                     silent = TRUE),
                 NA)
    expect_s3_class(output, "data.frame")
    expect_true("E" %in% names(output))
    expect_true("tau" %in% names(output))
    expect_true("tp" %in% names(output))
    expect_true("nn" %in% names(output))
    expect_true("k" %in% names(output))
    expect_true("num_pred" %in% names(output))
    expect_true("rho" %in% names(output))
    expect_true("mae" %in% names(output))
    expect_true("rmse" %in% names(output))
    expect_equal(NROW(output), 3)
    expect_equal(digest::digest(round(output$rho, 4)),
                 "605e8e3ede358655d719b453223f9c9a")
})

test_that("multiview model_output works", {
    expect_warning(output <- multiview(block, k = c(1, 3, "sqrt"),
                                       stats_only = FALSE, 
                                       silent = TRUE),
                   NA)
    expect_s3_class(output, "data.frame")
    expect_true("model_output" %in% names(output))
    expect_true(is.list(output$model_output))
    expect_error(model_output <- output$model_output[[3]], NA)
    expect_s3_class(model_output, "data.frame")
    expect_true("time" %in% names(model_output))
    expect_true("obs" %in% names(model_output))
    expect_true("pred" %in% names(model_output))
    expect_true("pred_var" %in% names(model_output))
    expect_equal(dim(model_output), c(100, 4))
    expect_equal(digest::digest(round(model_output, 4)),
                 "ab30b9bd5f3c69d2a867f7aceb92b203")
})