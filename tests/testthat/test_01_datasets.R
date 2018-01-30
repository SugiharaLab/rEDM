context("Check structure of included datasets")

test_that("two_species_model is correct", {
    expect_error(data("two_species_model"), NA)
    expect_true(exists("two_species_model"))
    expect_equal(dim(two_species_model), c(1000, 3))
    expect_equal(names(two_species_model), c("time", "x", "y"))
    expect_equal(digest::digest(two_species_model), 
                 "787da8c5a012f2a75e7c89b7fcc992b0")
})

test_that("block_3sp is correct", {
    expect_error(data("block_3sp"), NA)
    expect_true(exists("block_3sp"))
    expect_equal(dim(block_3sp), c(200, 10))
    expect_equal(names(block_3sp), c("time", "x_t", "x_t-1", "x_t-2", 
                                     "y_t", "y_t-1", "y_t-2", 
                                     "z_t", "z_t-1", "z_t-2"))
    expect_equal(digest::digest(block_3sp), 
                 "819609bd735c8ef02173b868c13055d4")
})

test_that("paramecium_didinium is correct", {
    expect_error(data("paramecium_didinium"), NA)
    expect_true(exists("paramecium_didinium"))
    expect_equal(dim(paramecium_didinium), c(71, 3))
    expect_equal(names(paramecium_didinium), c("time", "paramecium", 
                                               "didinium"))
    expect_equal(digest::digest(paramecium_didinium), 
                 "17e4ef1b9af8ce86fd8c8cc6df3c9a0d")
})

test_that("sardine_anchovy_sst is correct", {
    expect_error(data("sardine_anchovy_sst"), NA)
    expect_true(exists("sardine_anchovy_sst"))
    expect_equal(dim(sardine_anchovy_sst), c(78, 5))
    expect_equal(names(sardine_anchovy_sst), c("year", "anchovy", "sardine", 
                                               "sio_sst", "np_sst"))
    expect_equal(digest::digest(sardine_anchovy_sst), 
                 "5500a49be04dcdf771ddb362796d8211")
})

test_that("tentmap_del is correct", {
    expect_error(data("tentmap_del"), NA)
    expect_true(exists("tentmap_del"))
    expect_equal(length(tentmap_del), 999)
    expect_equal(digest::digest(tentmap_del), 
                 "bfe18a43a6f5cdcc38433a0761397b0b")
})

test_that("thrips_block is correct", {
    expect_error(data("thrips_block"), NA)
    expect_true(exists("thrips_block"))
    expect_equal(dim(thrips_block), c(81, 6))
    expect_equal(names(thrips_block), c("Year", "Month", "Thrips_imaginis", 
                                        "maxT_degC", "Rain_mm", "Season"))
    expect_equal(digest::digest(thrips_block), 
                 "7810f6159af2b81483d1146c3a89d1aa")
})

test_that("e120_invnit16 is correct", {
    expect_error(data("e120_invnit16"), NA)
    expect_true(exists("e120_invnit16"))
    expect_equal(dim(e120_invnit16), c(238, 7))
    expect_equal(names(e120_invnit16), c("Exp", "Year", "Plot", 
                                         "AbvBioAnnProd", "noh020tot", 
                                         "invrichness", "SummerPrecip.mm."))
    expect_equal(digest::digest(e120_invnit16), 
                 "c854a723a312ac2c5657199e9b6265d6")
})

test_that("sockeye_returns is correct", {
    expect_error(data("sockeye_returns"), NA)
    expect_true(exists("sockeye_returns"))
    expect_equal(dim(sockeye_returns), c(55, 10))
    expect_equal(names(sockeye_returns), c("year", "Early_Stuart", 
                                           "Late_Stuart", "Stellako", 
                                           "Quesnel", "Chilko", "Seymour", 
                                           "Late_Shuswap", "Birkenhead", 
                                           "Weaver"))
    expect_equal(digest::digest(sockeye_returns), 
                 "fb910d884cc0bca744bd2cc9979aad0a")
})