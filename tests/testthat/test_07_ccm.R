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
    expect_equal(digest::digest(ccm_out$rho), 
                 "072e0544e8765687ba7c56608e651357")
})

test_that("ccm_means works", {
    ccm_out <- ccm(sardine_anchovy_sst, E = 3, 
                   lib_sizes = seq(10, 80, by = 10), num_samples = 100, 
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
    # expect_equal(digest::digest(ccm_results), 
    #              "483eca2eef8b8baef8cc181d56a83ce8")
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