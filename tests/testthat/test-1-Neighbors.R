# The exclusion / compaction logic is backend-agnostic; exercise it with the
# exact brute-force backend so tests run without RANN. When RANN is present,
# the final block confirms the two backends agree.

embed <- matrix(0:9, ncol = 1)   # rows 1..10 hold values 0..9
BR    <- "brute"

test_that("disjoint lib/pred returns correct nearest neighbours (unsquared)", {
  fn <- FindNeighbors(embed, 1:5, 6:10, knn = 2, libOverlap = FALSE, backend = BR)
  expect_identical(fn$neighbors[1, ], c(5L, 4L))
  expect_equal(fn$distances[1, ], c(1, 2))
  expect_identical(fn$neighbors[5, ], c(5L, 4L))
})

test_that("libOverlap excludes the self-match", {
  fn <- FindNeighbors(embed, 1:10, 1:10, knn = 2, libOverlap = TRUE, backend = BR)
  expect_false(5L %in% fn$neighbors[5, ])
  expect_identical(sort(fn$neighbors[5, ]), c(4L, 6L))
  expect_equal(fn$distances[5, ], c(1, 1))
  expect_identical(fn$neighbors[1, ], c(2L, 3L))
})

test_that("exclusionRadius applies the Theiler window", {
  fn <- FindNeighbors(embed, 1:10, 1:10, knn = 2, exclusionRadius = 1,
                      libOverlap = TRUE, backend = BR)
  expect_false(any(c(4L, 5L, 6L) %in% fn$neighbors[5, ]))
  expect_identical(sort(fn$neighbors[5, ]), c(3L, 7L))
  expect_equal(fn$distances[5, ], c(2, 2))
  expect_identical(fn$neighbors[1, ], c(3L, 4L))
})

test_that("validLib restricts the library, order preserved", {
  vlib <- rep(FALSE, 10); vlib[c(1, 2, 3, 8, 9, 10)] <- TRUE
  fn <- FindNeighbors(embed, 1:10, 5, knn = 2, validLib = vlib,
                      libOverlap = TRUE, backend = BR)
  expect_true(all(fn$neighbors[1, ] %in% which(vlib)))
  expect_identical(fn$neighbors[1, ], c(3L, 2L))
  expect_equal(fn$distances[1, ], c(2, 3))
})

test_that("deficiency yields Inf-padded neighbours without leakage", {
  # pred = 2, lib = 1:4, exclusionRadius = 5 excludes the whole library:
  # the row is fully excluded and must be Inf-padded (issue #74), not
  # back-filled with the raw nearest (which would reinstate the self-match).
  expect_warning(
    fn <- FindNeighbors(embed, 1:4, 2, knn = 3, exclusionRadius = 5,
                        libOverlap = TRUE, backend = BR, verbose = TRUE))
  expect_equal(ncol(fn$neighbors), 3L)                 # shape preserved
  expect_true(all(is.infinite(fn$distances[1, ])))     # fully excluded -> Inf
  # no finite-distance neighbour lies within the exclusion radius (no leak)
  fin <- is.finite(fn$distances[1, ])
  expect_false(any(fin & abs(2 - fn$neighbors[1, ]) <= 5))
})

test_that("partial deficiency keeps valid neighbours, Inf-pads the rest", {
  # pred = 2, lib = 1:6, exclusionRadius = 2 leaves only rows 5,6 valid for
  # knn = 5 : two finite neighbours followed by three Inf-padded slots.
  fn <- FindNeighbors(embed, 1:6, 2, knn = 5, exclusionRadius = 2,
                      libOverlap = TRUE, backend = BR)
  fin <- is.finite(fn$distances[1, ])
  expect_equal(sum(fin), 2L)                           # only 2 valid
  expect_true(all(is.finite(fn$distances[1, 1:2])))    # finite prefix
  expect_true(all(is.infinite(fn$distances[1, 3:5])))  # Inf padding trails
  expect_identical(sort(fn$neighbors[1, fin]), c(5L, 6L))
  expect_false(any(fin & abs(2 - fn$neighbors[1, ]) <= 2))  # no leak
})

test_that("KNNQuery emits sentinels when k > nLib", {
  q <- KNNQuery(matrix(c(0, 1), ncol = 1), matrix(0, ncol = 1), k = 3, backend = BR)
  expect_identical(q$nnIdx[1, 3], 0L)
  expect_true(is.infinite(q$nnDist[1, 3]))
  expect_identical(q$nnIdx[1, 1:2], c(1L, 2L))
})

test_that("degenerate library caps kOut at min(knn, kQuery)", {
  fn <- FindNeighbors(embed, 1:2, 1:2, knn = 3, libOverlap = TRUE, backend = BR)
  expect_equal(ncol(fn$neighbors), 2L)
})

test_that("distances are unsquared Euclidean in 2-D", {
  embed2 <- matrix(c(0, 0, 1, 1, 3, 3), ncol = 2, byrow = TRUE)
  fn <- FindNeighbors(embed2, 1:3, 2, knn = 1, libOverlap = TRUE, backend = BR)
  expect_identical(fn$neighbors[1, 1], 1L)
  expect_equal(fn$distances[1, 1], sqrt(2))
})

test_that("RANN and brute backends agree when RANN is installed", {
  skip_if_not_installed("RANN")
  for (over in c(FALSE, TRUE)) {
    lib <- if (over) 1:10 else 1:5
    prd <- if (over) 1:10 else 6:10
    a <- FindNeighbors(embed, lib, prd, knn = 2, exclusionRadius = 1,
                       libOverlap = over, backend = "brute")
    b <- FindNeighbors(embed, lib, prd, knn = 2, exclusionRadius = 1,
                       libOverlap = over, backend = "RANN")
    expect_equal(a$distances, b$distances)   # neighbour identity may tie-differ
  }
})
