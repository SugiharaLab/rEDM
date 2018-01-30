context("Check ccm function")

data("sardine_anchovy_sst")

test_that("ccm works", {
    expect_error(ccm_out <- ccm(sardine_anchovy_sst, E = 3, 
                                lib_sizes = seq(10, 80, by = 10), num_samples = 100, 
                                lib_column = "anchovy", target_column = "np_sst",
                                RNGseed = 42, silent = TRUE),
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
    expect_equal(NROW(ccm_out), 800)
})

test_that("ccm_means works", {
    expect_error(ccm_means <- ccm_means(ccm_out), NA)
    expect_s3_class(ccm_means, "data.frame")
    expect_true("E" %in% names(ccm_means))
    expect_true("tau" %in% names(ccm_means))
    expect_true("tp" %in% names(ccm_means))
    expect_true("nn" %in% names(ccm_means))
    expect_true("lib_column" %in% names(ccm_means))
    expect_true("target_column" %in% names(ccm_means))
    expect_true("lib_size" %in% names(ccm_means))
    expect_equal(anyDuplicated(ccm_means$lib_size), 0)
    expect_true("num_pred" %in% names(ccm_means))
    expect_true("rho" %in% names(ccm_means))
    expect_true("mae" %in% names(ccm_means))
    expect_true("rmse" %in% names(ccm_means))
    expect_equal(NROW(ccm_means), 8)
})

test_that("ccm error checking works", {
    df <- data.frame(a = rnorm(5), b = rnorm(5))
    expect_warning(ccm(df))
    expect_error(suppressWarnings(ccm(df, E = 6, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, E = 2, tau = 5, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, E = 1, tp = 5, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, E = 1, tp = -5, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, lib_column = "c", silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, lib_column = 0, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, lib_column = 3, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, target_column = "c", silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, target_column = 0, silent = TRUE)))
    expect_error(suppressWarnings(ccm(df, target_column = 3, silent = TRUE)))
})