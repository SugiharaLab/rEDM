#------------------------------------------------------------------------
# Phase 3 : S-map (sequential locally weighted global linear map).
#
# Port of pyEDM SMap (SMap.py, EDM.py). See spec: rEDM_SMap_solver_spec.md.
#
# Solver replicates numpy.linalg.lstsq(rcond = None): SVD least squares with
# singular-value cutoff eps * max(M, N) * s_max (eps = .Machine$double.eps).
# Theta weighting uses the row MEAN distance (not Simplex 1e-6 floor).
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' S-map projection (top-level API).
#'
#' @inheritParams Simplex
#' @param theta local weighting parameter (0 = global linear map).
#' @return list(predictions, coefficients, singularValues) data.frames.
#' @keywords internal
#------------------------------------------------------------------------
SMap <- function(dataFrame = NULL, columns, target, lib, pred,
                 E = NULL, Tp = 1, knn = 0, tau = -1, theta = 0,
                 exclusionRadius = 0, embedded = FALSE, validLib = logical(0),
                 noTime = FALSE, ignoreNan = TRUE, backend = "RANN",
                 pathIn = "./", dataFile = "", pathOut = "./", predictFile = "",
                 parameterList = FALSE, includeState = FALSE,
                 verbose = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame", "includeState", "pathIn", "dataFile",
               "pathOut", "predictFile", "parameterList", "showPlot",
               "verbose", "backend")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  columns   <- SplitColumns(columns)
  target    <- SplitColumns(target)
  nRows     <- nrow(dataFrame)

  E      <- RequireE("SMap", E, embedded, columns) # required unless embedded
  idx    <- CreateIndices(lib, pred, E, tau, Tp, embedded, nRows, "SMap")
  lib_i  <- idx$lib_i
  pred_i <- idx$pred_i

  originalLibLen <- length(lib_i)
  if (knn < 1) { knn <- originalLibLen - 1L } # SMap default: full library
  originalKnn <- knn

  targetVec <- dataFrame[[target[1]]]
  timeVec   <- if (noTime) seq_len(nRows) else dataFrame[[1]]

  embedding <- if (embedded) {
                 as.matrix(dataFrame[, columns, drop = FALSE])
               }
               else {
                 MakeBlock(dataFrame, E, tau, columns)
               }

  # Drop rows with NaN in the embedding from the lib and pred index sets.
  rn <- RemoveNan(embedding, lib_i, pred_i, ignoreNan)

  # S-map's default knn is the full library (originalLibLen - 1). If NaN
  # removal shrank the library, that default is now too large, adjust
  # to the reduced library. A user-specified knn is left untouched.
  libraryShrank     <- length(rn$lib_i) < originalLibLen
  knnWasFullLibrary <- originalKnn == (originalLibLen - 1L)

  if (isTRUE(ignoreNan) && libraryShrank && knnWasFullLibrary) {
    # re-apply "whole library" to the smaller library
    knn <- length(rn$lib_i) - 1L
  }

  lib_i  <- rn$lib_i
  pred_i <- rn$pred_i

  nb <- FindNeighbors(embedding, lib_i, pred_i, knn,
                      exclusionRadius = exclusionRadius, validLib = validLib,
                      libOverlap = idx$libOverlap, backend = backend,
                      tieBreak = FALSE, verbose = verbose) # SMap: no tie-break

  # Flag if nan are present in targetVec for solver filtering
  targetVecNan <- any(is.na(targetVec))

  pr <- SMapProject(nb$neighbors, nb$distances, embedding, targetVec,
                    pred_i, E, Tp, theta, knn, targetVecNan)

  timeName = "Time"
  if (!noTime && !is.null(names(dataFrame))) {
    timeName <- names(dataFrame)[1]
  }

  predictions <- FormatProjection(idx$predList, pred_i, idx$pred_i_all,
                                  targetVec, timeVec, pr$projection,
                                  pr$variance, Tp, nRows, tau, timeName)

  # Coefficients & singular-value frames share the projection layout
  L <- ProjectionLayout(idx$predList, pred_i, idx$pred_i_all,
                        timeVec, Tp, nRows, tau)
  nDim     <- E + 1
  embNames <- EmbedColumnNames(if (embedded) columns
                               else columns, E, tau, embedded)
  coefNames <- paste0("\u2202", target[1], "/\u2202", embNames)

  coefOut <- matrix(NA_real_, L$outSize, nDim)
  svOut   <- matrix(NA_real_, L$outSize, nDim)
  coefOut[L$predOut_i + 1L, ] <- pr$coefficients
  svOut[L$predOut_i + 1L, ]   <- pr$singularValues

  coefDf <- data.frame(Time = L$timeOut, coefOut, check.names = FALSE,
                       stringsAsFactors = FALSE)
  names(coefDf) <- c("Time", "C0", coefNames)

  svDf <- data.frame(Time = L$timeOut, svOut, check.names = FALSE,
                     stringsAsFactors = FALSE)
  names(svDf) <- c("Time", paste0("C", 0:(nDim - 1)))

  if (showPlot) {
    PlotObsPred(predictions, "", E, Tp)
  }

  result <- list(predictions = predictions, coefficients = coefDf,
                 singularValues = svDf)

  if (showPlot) {
    PlotSmap(result, "", E, Tp)
  }

  internal <- if (isTRUE(includeState)) {
                  list(knn_neighbors = nb$neighbors,
                       knn_distances = nb$distances,
                       lib_i         = lib_i,
                       pred_i        = pred_i,
                       targetVec     = targetVec,
                       embedding     = embedding)
              }
              else { NULL }
    
    result = EdmFinalize(result, result $ predictions,
                         pathOut, predictFile, parameterList, parameters,
                         includeState = includeState, internal = internal)
}

#------------------------------------------------------------------------
#' Minimum-norm least squares via truncated SVD.
#'
#' @param A design matrix (m x n).
#' @param b response vector (length m).
#' @return list(C = coefficients length n, SV = singular values length min(m,n)).
#' @details matches numpy.linalg.lstsq( rcond=(eps * max(m, n)) )
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
#
# Reference is pyEDM's numpy.linalg.lstsq, a wrapper over LAPACK gelsd() 
# (SVD-based, singular-value truncation, minimum-norm solution). R has
# no base equivalent that both solves and returns singular values, so we
#   First: compute SVD with svd() -- a LAPACK dgesdd call -- and
#   Next:  reconstruct the least-squares solution by hand.
# Steps after svd() below are the extra work gelsd does internally.
#
# SMap is deliberately ill-conditioned at large theta, where SVD's
# minimum-norm truncation is what matches gelsd/pyEDM. The cutoff
# reproduces numpy.linalg.lstsq's default rcond, which is what makes the
# coefficient and singularValue fixtures agree bit-for-bit. Do not swap the
# solver or drop the cutoff without re-baselining those fixtures.
#
# Returns: C  = fit coefficients (pyEDM 'coefficients')
#          SV = singular values of A (pyEDM 'singularValues'), raw from LAPACK
#------------------------------------------------------------------------
SMapSolve <- function(A, b) {
  # (1) LAPACK dgesdd: A = u diag(d) t(v). Decomposition only, no solve.
  s <- svd(A)
  d <- s$d    # singular values, descending

  # (2) rcond-style threshold ~ numpy lstsq default
  cutoff <- .Machine$double.eps * max(nrow(A), ncol(A)) * d[1]

  # (3) truncated reciprocal: zero below cutoff -> minimum-norm, rank-safe
  dInv <- ifelse(d > cutoff, 1 / d, 0)

  # (4) reconstruct x = V diag(dInv) t(U) b  (the solve gelsd does internally)
  C <- s$v %*% (dInv * crossprod(s$u, b))

  list(C = as.numeric(C), SV = d)
}

#------------------------------------------------------------------------
#' S-map projection kernel (SMap.py::Project).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
SMapProject <- function(knnNeighbors, knnDistances, embedding, targetVec,
                        pred_i, E, Tp, theta, knn, targetVecNan) {
  nPred <- length(pred_i)
  nDim  <- E + 1

  projection     <- rep(NA_real_, nPred)
  variance       <- rep(NA_real_, nPred)
  coefficients   <- matrix(NA_real_, nPred, nDim)
  singularValues <- matrix(NA_real_, nPred, nDim)

  # Mean neighbour distance per row over FINITE (valid) neighbours only;
  # Inf padding slots are excluded. Rows with no finite neighbour give NaN,
  # which flags an empty (skipped) row below (issue #74).
  finiteD     <- is.finite(knnDistances)
  finiteCnt   <- rowSums(finiteD)
  distSum     <- rowSums(ifelse(finiteD, knnDistances, 0))
  distRowMean <- ifelse(finiteCnt > 0, distSum / finiteCnt, NaN)
  if (theta == 0) {
    W <- matrix(1, nrow(knnDistances), ncol(knnDistances))
  }
  else {
    W <- exp(-(theta / distRowMean) * knnDistances)   # row-scaled
  }

  neighborsTp <- knnNeighbors + Tp
  n <- length(targetVec)
  B <- matrix(NA_real_, nPred, ncol(knnNeighbors))
  inRange <- neighborsTp >= 1 & neighborsTp <= n
  B[inRange] <- targetVec[neighborsTp[inRange]]

  wB <- W * B

  for (row in seq_len(nPred)) {
    # Empty row : no finite-distance neighbour survived exclusion / validLib.
    # Leave projection / coefficients / variance at their NA defaults.
    if (is.nan(distRowMean[row])) { next }

    # Valid slots : finite distance (exclude Inf padding) and, when the
    # target has NA, finite target. Gating on distance finiteness (not the
    # weight) makes padding inert for all theta, so theta == 0 needs no
    # special case (issue #74).
    valid <- is.finite(knnDistances[row, ])
    if (targetVecNan) { valid <- valid & is.finite(B[row, ]) }
    if (!any(valid))  { next }

    # Build A over valid slots only : knnNeighbors padding carries the 0
    # sentinel index, which must never index the embedding.
    libRows <- knnNeighbors[row, valid]
    m       <- sum(valid)
    A       <- matrix(NA_real_, m, nDim)
    A[, 1]  <- W[row, valid]
      
    for (j in 2:nDim) {
      A[, j] <- W[row, valid] * embedding[libRows, j - 1]
    }
    wBrow <- wB[row, valid]

    sol <- SMapSolve(A, wBrow)
    C   <- sol$C
    coefficients[row, ]   <- C
    singularValues[row, ] <- c(sol$SV, rep(NA_real_, nDim - length(sol$SV)))

    proj <- if (is.na(C[1])) { 0 }
            else             { C[1] }

    for (e in 2:nDim) {
      proj <- proj + C[e] * embedding[pred_i[row], e - 1]
    }

    projection[row] <- proj

    wRow          <- W[row, valid]
    deltaSqr      <- (B[row, valid] - proj)^2
    variance[row] <- sum(wRow * deltaSqr) / sum(wRow)
  }

  list(projection = projection, variance = variance,
       coefficients = coefficients, singularValues = singularValues)
}
