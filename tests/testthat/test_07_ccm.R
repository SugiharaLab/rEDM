context("Check ccm function")

data("sardine_anchovy_sst")

test_that("ccm works", {
    expect_error(ccm_out <- ccm(sardine_anchovy_sst, E = 3, 
                                lib_sizes = seq(10, 80, by = 10), 
                                lib_column = "anchovy", target_column = "np_sst",
                                random_libs = FALSE, silent = TRUE),
                 NA)
    expect_s3_class(ccm_out, "data.frame")
    expect_true("E" %in% names(ccm_out))
    expect_true("tau" %in% names(ccm_out))
    expect_true("tp" %in% names(ccm_out))
    expect_true("nn" %in% names(ccm_out))
    expect_true("lib_column" %in% names(ccm_out))
    expect_true("target_column" %in% names(ccm_out))
    expect_true("lib_size" %in% names(ccm_out))
    expect_true("num_pred" %in% names(ccm_out))
    expect_true("rho" %in% names(ccm_out))
    expect_true("mae" %in% names(ccm_out))
    expect_true("rmse" %in% names(ccm_out))
    expect_equal(NROW(ccm_out), 533)
    expect_equal(digest::digest(round(ccm_out$rho, 4)), 
                 "405d54f5c1f6bd46131339134c00a539")
})

test_that("ccm_means works", {
    ccm_out <- ccm(sardine_anchovy_sst, E = 3, 
                   lib_sizes = seq(10, 80, by = 10), 
                   lib_column = "anchovy", target_column = "np_sst",
                   random_libs = FALSE, silent = TRUE)
    expect_error(ccm_results <- ccm_means(ccm_out), NA)
    expect_s3_class(ccm_results, "data.frame")
    expect_true("E" %in% names(ccm_results))
    expect_true("tau" %in% names(ccm_results))
    expect_true("tp" %in% names(ccm_results))
    expect_true("nn" %in% names(ccm_results))
    expect_true("lib_column" %in% names(ccm_results))
    expect_true("target_column" %in% names(ccm_results))
    expect_true("lib_size" %in% names(ccm_results))
    expect_equal(anyDuplicated(ccm_results$lib_size), 0)
    expect_true("num_pred" %in% names(ccm_results))
    expect_true("rho" %in% names(ccm_results))
    expect_true("mae" %in% names(ccm_results))
    expect_true("rmse" %in% names(ccm_results))
    expect_equal(NROW(ccm_results), 8)
    expect_equal(digest::digest(round(ccm_results$rho, 4)), 
                  "c966dd8ec76cb29d1fc0dc34b64638ea")
})

test_that("ccm model_output works", {
    expect_error(ccm_out <- ccm(sardine_anchovy_sst, E = 3, 
                                lib_sizes = seq(10, 80, by = 10), 
                                lib_column = "anchovy", target_column = "np_sst",
                                stats_only = FALSE, 
                                random_libs = FALSE, silent = TRUE),
                 NA)
    expect_s3_class(ccm_out, "data.frame")
    expect_true("model_output" %in% names(ccm_out))
    expect_true(is.list(ccm_out$model_output))
    expect_error(model_output <- ccm_out$model_output[[1]], NA)
    expect_s3_class(model_output, "data.frame")
    expect_true("time" %in% names(model_output))
    expect_true("obs" %in% names(model_output))
    expect_true("pred" %in% names(model_output))
    expect_true("pred_var" %in% names(model_output))
    
    # check that number of rows of model output match up with number of predictions
    expect_equal(vapply(ccm_out$model_output, NROW, 1), 
                 rep.int(78, times = 533))
    expect_equal(vapply(ccm_out$model_output, NCOL, 1), 
                 rep.int(4, times = 533))
    
    # check a few digests
    idx <- 1
    model_output <- ccm_out$model_output[[idx]]
    model_stats <- compute_stats(model_output$obs, model_output$pred)
    expect_equal(model_stats[, c("num_pred", "rho", "mae", "rmse")], 
                 ccm_out[idx, c("num_pred", "rho", "mae", "rmse")], 
                 check.attributes = FALSE)
    expect_equal(digest::digest(round(model_output, 4)), 
                 "c4725803e1b8974909c8128475b5d463")
    
    idx <- 533
    model_output <- ccm_out$model_output[[idx]]
    model_stats <- compute_stats(model_output$obs, model_output$pred)
    expect_equal(model_stats[, c("num_pred", "rho", "mae", "rmse")], 
                 ccm_out[idx, c("num_pred", "rho", "mae", "rmse")], 
                 check.attributes = FALSE)
    expect_equal(digest::digest(round(model_output, 4)), 
                 "7608d92d62c38edf583730e720635730")
    
    ### add test for ccm_means on ccm_output with model_output
    expect_error(ccm_results <- ccm_means(ccm_out), NA)
    expect_s3_class(ccm_results, "data.frame")
    expect_true("E" %in% names(ccm_results))
    expect_true("tau" %in% names(ccm_results))
    expect_true("tp" %in% names(ccm_results))
    expect_true("nn" %in% names(ccm_results))
    expect_true("lib_column" %in% names(ccm_results))
    expect_true("target_column" %in% names(ccm_results))
    expect_true("lib_size" %in% names(ccm_results))
    expect_equal(anyDuplicated(ccm_results$lib_size), 0)
    expect_true("num_pred" %in% names(ccm_results))
    expect_true("rho" %in% names(ccm_results))
    expect_true("mae" %in% names(ccm_results))
    expect_true("rmse" %in% names(ccm_results))
    expect_equal(NROW(ccm_results), 8)
    expect_equal(digest::digest(round(ccm_results$rho, 4)), 
                 "c966dd8ec76cb29d1fc0dc34b64638ea")
})

test_that("ccm works on multivariate time series", {
    expect_warning(output <- ccm(EuStockMarkets[1:300, ], 
                                 lib_column = "DAX",
                                 target_column = "CAC", 
                                 random_libs = FALSE))
    
    output <- data.frame(lapply(output, function(y) 
        if (is.numeric(y)) round(y, 4) else y))
    attributes(output) <- attributes(output)[sort(names(attributes(output)))]
    expect_equal(digest::digest(output), "004d6d6aa7e57fed21ec9385b02dac03")
})

test_that("ccm error checking works", {
    df <- data.frame(a = 1:5, b = 0:4)
    expect_warning(ccm(df))
    expect_warning(ccm(df, lib_sizes = c(2, 2, 3)))
    expect_error(ccm(df, E = 6, silent = TRUE))
    expect_error(ccm(df, E = 2, tau = 5, silent = TRUE))
    expect_error(ccm(df, E = 1, tp = 5, silent = TRUE))
    expect_error(ccm(df, E = 1, tp = -5, silent = TRUE))
    expect_error(ccm(df, lib_column = "c", silent = TRUE))
    expect_error(ccm(df, lib_column = 0, silent = TRUE))
    expect_error(ccm(df, lib_column = 3, silent = TRUE))
    expect_error(ccm(df, target_column = "c", silent = TRUE))
    expect_error(ccm(df, target_column = 0, silent = TRUE))
    expect_error(ccm(df, target_column = 3, silent = TRUE))
    expect_error(ccm(df, lib_sizes = -1, silent = TRUE))
})