#------------------------------------------------------------------------
# Phase 3 : S-map (sequential locally weighted global linear map).
#
# Port of pyEDM SMap (SMap.py, EDM.py). See spec: rEDM_SMap_solver_spec.md.
#
# Solver replicates numpy.linalg.lstsq(rcond = None): SVD least squares with
# singular-value cutoff eps * max(M, N) * s_max (eps = .Machine$double.eps).
# Theta weighting uses the row MEAN distance (no 1e-6 floor; that floor is
# Simplex-only).
#
# Naming: functions/constants uppercase-first, variables lowercase-first.
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
                 E = 0, Tp = 1, knn = 0, tau = -1, theta = 0,
                 exclusionRadius = 0, embedded = FALSE, validLib = logical(0),
                 noTime = FALSE, ignoreNan = TRUE, backend = "RANN",
                 pathIn = "./", dataFile = "", pathOut = "./", predictFile = "",
                 parameterList = FALSE, verbose = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  columns   <- SplitColumns(columns)
  target    <- SplitColumns(target)
  nRows     <- nrow(dataFrame)

  if (embedded) E <- length(columns)  # Any user supplied E ignored
  idx    <- CreateIndices(lib, pred, E, tau, Tp, embedded, nRows, "SMap")
  lib_i  <- idx$lib_i
  pred_i <- idx$pred_i

  originalLibLen   <- length(lib_i)
  if (knn < 1) knn <- originalLibLen - 1L          # SMap default: full library
  originalKnn      <- knn

  targetVec <- dataFrame[[target[1]]]
  timeVec   <- if (noTime) seq_len(nRows) else dataFrame[[1]]

  embedding <- if (embedded) as.matrix(dataFrame[, columns, drop = FALSE])
               else MakeBlock(dataFrame, E, tau, columns)

  rn <- RemoveNan(embedding, lib_i, pred_i, ignoreNan)
  if (ignoreNan && length(rn$lib_i) < originalLibLen &&
      originalKnn == originalLibLen - 1L)
  knn    <- length(rn$lib_i) - 1L                  # track reduced library
  lib_i  <- rn$lib_i
  pred_i <- rn$pred_i

  targetVecNan <- any(is.na(targetVec))

  nb <- FindNeighbors(embedding, lib_i, pred_i, knn,
                      exclusionRadius = exclusionRadius, validLib = validLib,
                      libOverlap = idx$libOverlap, backend = backend,
                      tieBreak = FALSE, verbose = verbose)   # SMap: no tie-break

  pr <- SMapProject(nb$neighbors, nb$distances, embedding, targetVec,
                    pred_i, E, Tp, theta, knn, targetVecNan)

  timeName = "Time"
  if (!noTime && !is.null(names(dataFrame))) timeName <- names(dataFrame)[1]

  predictions <- FormatProjection(idx$predList, pred_i, idx$pred_i_all,
                                  targetVec, timeVec, pr$projection,
                                  pr$variance, Tp, nRows, tau, timeName)

  # Coefficients & singular-value frames share the projection layout
  L <- ProjectionLayout(idx$predList, pred_i, idx$pred_i_all, timeVec,
                        Tp, nRows, tau)
  nDim     <- E + 1
  embNames <- EmbedColumnNames(if (embedded) columns else columns, E, tau, embedded)
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

  if (showPlot) PlotObsPred(predictions, "", E, Tp)
  result <- list(predictions = predictions, coefficients = coefDf,
                 singularValues = svDf)
  if (showPlot) PlotSmap(result, "", E, Tp)
    result = EdmFinalize(result, result $ predictions,
                         pathOut, predictFile, parameterList, parameters)
}

#------------------------------------------------------------------------
#' Minimum-norm least squares via SVD, matching numpy.linalg.lstsq(rcond=NULL).
#'
#' @param A design matrix (m x n).
#' @param b response vector (length m).
#' @return list(C = coefficients length n, SV = singular values length min(m,n)).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
SMapSolve <- function(A, b) {
  s <- svd(A)                                   # A = u diag(d) t(v)
  d <- s$d
  cutoff <- .Machine$double.eps * max(nrow(A), ncol(A)) * d[1]
  dInv <- ifelse(d > cutoff, 1 / d, 0)
  C <- s$v %*% (dInv * crossprod(s$u, b))       # v diag(dInv) t(u) b
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

  distRowMean <- rowMeans(knnDistances)
  if (theta == 0) {
    W <- matrix(1, nrow(knnDistances), ncol(knnDistances))
  } else {
    W <- exp(-(theta / distRowMean) * knnDistances)   # row-scaled
  }

  neighborsTp <- knnNeighbors + Tp
  n <- length(targetVec)
  B <- matrix(NA_real_, nPred, ncol(knnNeighbors))
  inRange <- neighborsTp >= 1 & neighborsTp <= n
  B[inRange] <- targetVec[neighborsTp[inRange]]

  wB <- W * B

  for (row in seq_len(nPred)) {
    libRows <- knnNeighbors[row, ]
    A <- matrix(NA_real_, knn, nDim)
    A[, 1] <- W[row, ]
    for (j in 2:nDim) A[, j] <- W[row, ] * embedding[libRows, j - 1]
    wBrow <- wB[row, ]

    valid <- rep(TRUE, knn)
    if (targetVecNan) {
      valid <- is.finite(B[row, ])
      A     <- A[valid, , drop = FALSE]
      wBrow <- wB[row, valid]
    }

    sol <- SMapSolve(A, wBrow)
    C   <- sol$C
    coefficients[row, ]   <- C
    singularValues[row, ] <- c(sol$SV, rep(NA_real_, nDim - length(sol$SV)))

    proj <- if (is.na(C[1])) 0 else C[1]
    for (e in 2:nDim) proj <- proj + C[e] * embedding[pred_i[row], e - 1]
    projection[row] <- proj

    wRow     <- W[row, valid]
    deltaSqr <- (B[row, valid] - proj)^2
    variance[row] <- sum(wRow * deltaSqr) / sum(wRow)
  }

  list(projection = projection, variance = variance,
       coefficients = coefficients, singularValues = singularValues)
}
