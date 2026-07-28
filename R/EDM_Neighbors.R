#------------------------------------------------------------------------
# Nearest-neighbour search and exclusion.
#
# Port of pyEDM/src/pyEDM/Neighbors.py::FindNeighbors (2.5.0).
# See spec: rEDM_Neighbors_spec.md.
#
# Naming follows rEDM convention: functions/constants uppercase-first,
# variables lowercase-first.
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Find k nearest library neighbours for each prediction point.
#'
#' @param embedding      numeric matrix (N x E), row-indexed by data row.
#' @param libRows        integer vector of 1-based library data-row indices.
#' @param predRows       integer vector of 1-based prediction data-row indices.
#' @param knn            neighbours per prediction (assumed already defaulted).
#' @param exclusionRadius Theiler window in rows (>= 0).
#' @param validLib       logical vector (length N) of admissible library rows,
#'                       or length-0 for "all".
#' @param libOverlap     TRUE if libRows and predRows intersect.
#' @param xRadKnnFactor  over-query multiplier for exclusionRadius (default 5).
#' @param backend        "RANN" or "brute".
#' @param verbose        emit progress / deficiency warnings.
#'
#' @return list(neighbors, distances): each an (nPred x knn) matrix of
#'   1-based library data-row indices and unsquared Euclidean distances,
#'   ascending by distance per row, after exclusion.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
FindNeighbors <- function(embedding, libRows, predRows, knn,
                          exclusionRadius = 0, validLib = logical(0),
                          libOverlap = FALSE, xRadKnnFactor = 5,
                          backend = c("RANN", "brute"), verbose = FALSE) {
  backend <- match.arg(backend)
  nPred   <- length(predRows)

  #-----------------------------------------------
  # (3.1) Determine if exclusionRadius filtering needed
  #-----------------------------------------------
  exclusionRadiusKnn <- FALSE
  if (exclusionRadius > 0) {
    if (libOverlap) {
      exclusionRadiusKnn <- TRUE
    } else {
      excludeRow <- 0
      nLibCur    <- length(libRows)
      if (predRows[1] > libRows[nLibCur]) {
        excludeRow <- predRows[1] - libRows[nLibCur]    # pred start beyond lib end
      } else if (libRows[1] > predRows[nPred]) {
        excludeRow <- libRows[1] - predRows[nPred]      # lib start beyond pred end
      }
      if (exclusionRadius >= excludeRow) exclusionRadiusKnn <- TRUE
    }
  }

  #-----------------------------------------------
  # (3.2) Filter library by validLib (order preserving)
  #-----------------------------------------------
  if (length(validLib)) {
    validLibIdx <- which(validLib)                       # 1-based data rows
    libRows     <- libRows[libRows %in% validLibIdx]
    if (length(libRows) == 0L) {
      stop("FindNeighbors: no valid library points; all excluded by validLib.")
    }
    if (length(libRows) < knn && verbose) {
      warning(sprintf("FindNeighbors: only %d valid library points, but knn = %d.",
                      length(libRows), knn))
    }
  }
  nLib <- length(libRows)

  #-----------------------------------------------
  # (3.3) Determine kQuery : neighbours to request
  #-----------------------------------------------
  kQuery <- knn
  if (exclusionRadiusKnn) {
    kQuery <- min(knn * xRadKnnFactor, nLib)
  } else if (libOverlap) {
    kQuery <- knn + 1L
  }
  if (length(validLib)) {
    kQuery <- nLib                                            # examine all
  }
  kQuery <- max(1L, min(kQuery, nLib))                        # never exceed nLib

  #-----------------------------------------------
  # (3.4) Query
  #-----------------------------------------------
  libData  <- embedding[libRows,  , drop = FALSE]
  predData <- embedding[predRows, , drop = FALSE]
  q        <- KNNQuery(libData, predData, kQuery, backend = backend)
  nnIdx    <- q$nnIdx       # (nPred x kQuery) 1-based into libRows, 0 = sentinel
  nnDist   <- q$nnDist

  #-----------------------------------------------
  # (3.5, 3.6) Sentinel guard + remap to data rows
  #-----------------------------------------------
  sentinel  <- nnIdx == 0L
  neighbors <- matrix(0L, nrow(nnIdx), ncol(nnIdx))
  if (any(!sentinel)) {
    neighbors[!sentinel] <- libRows[nnIdx[!sentinel]]        # 1-based, no offset
  }
  distances <- nnDist
  distances[sentinel] <- Inf

  #-----------------------------------------------
  # (3.7) Exclusion mask + compaction
  #-----------------------------------------------
  needsFiltering <- libOverlap || exclusionRadiusKnn ||
                    kQuery > knn || any(sentinel)

  if (!needsFiltering) {
    return(list(neighbors = neighbors, distances = distances))
  }

  predCol <- matrix(predRows, nPred, kQuery)               # broadcast down columns

  if (exclusionRadiusKnn) {
    mask <- abs(predCol - neighbors) <= exclusionRadius    # subsumes self-match
  } else if (libOverlap) {
    mask <- predCol == neighbors                           # self-match only
  } else {
    mask <- matrix(FALSE, nPred, kQuery)                   # validLib trim only
  }
  mask <- mask | sentinel                                  # always exclude sentinels

  valid   <- !mask
  cs      <- RowCumsum(valid)
  firstK  <- valid & (cs <= knn)

  # Deficiency: fewer than knn valid neighbours in a row
  validCounts <- cs[, kQuery]
  deficient   <- validCounts < knn
  if (any(deficient)) {
    if (verbose) {
      warning(sprintf(
        "FindNeighbors: failed to find knn=%d outside exclusionRadius=%d for %d prediction(s); consider reducing knn.",
        knn, exclusionRadius, sum(deficient)))
    }
    fallbackCols <- seq_len(min(knn, kQuery))
    for (i in which(deficient)) {
      firstK[i, ]            <- FALSE
      firstK[i, fallbackCols] <- TRUE                        # raw nearest knn
    }
  }

  # Compact: leftmost selected columns per row (stable), take first kOut
  kOut         <- min(knn, kQuery)
  outNeighbors <- matrix(0L,  nPred, kOut)
  outDist      <- matrix(Inf, nPred, kOut)
  for (i in seq_len(nPred)) {
    sel <- which(firstK[i, ])
    if (length(sel) > kOut) sel <- sel[seq_len(kOut)]
    m <- length(sel)
    if (m) {
      outNeighbors[i, seq_len(m)] <- neighbors[i, sel]
      outDist[i,      seq_len(m)] <- distances[i, sel]
    }
  }

  list(neighbors = outNeighbors, distances = outDist)
}

#------------------------------------------------------------------------
#' Exact k nearest-neighbour query (backend abstraction).
#'
#' Returns 1-based indices into \code{libData} and unsquared Euclidean
#' distances, ascending per row. When \code{k} exceeds the number of
#' library rows the surplus columns are filled with the "no neighbour"
#' sentinel (index 0, distance Inf), matching RANN's convention.
#'
#' @param libData  numeric matrix (nLib x E), the library embedding.
#' @param predData numeric matrix (nPred x E), the query embedding.
#' @param k        neighbours to request.
#' @param backend  "RANN" (production) or "brute" (exact base-R fallback).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
KNNQuery <- function(libData, predData, k, backend = c("RANN", "brute")) {
  backend <- match.arg(backend)
  nLib    <- nrow(libData)
  nPred   <- nrow(predData)
  kEff    <- min(k, nLib)

  nnIdx  <- matrix(0L,  nPred, k)   # sentinel index 0
  nnDist <- matrix(Inf, nPred, k)   # sentinel distance Inf

  if (backend == "RANN") {
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop("KNNQuery: package 'RANN' not available; ",
           "install RANN or call with backend = 'brute'.")
    }
    res <- RANN::nn2(data = libData, query = predData, k = kEff,
                     treetype = "kd", searchtype = "standard", eps = 0)
    nnIdx[,  seq_len(kEff)] <- res$nn.idx
    nnDist[, seq_len(kEff)] <- res$nn.dists
  } else {
    # Vectorised exact brute force: d^2 = |p|^2 + |l|^2 - 2 p.l
    libSq  <- rowSums(libData  * libData)          # nLib
    predSq <- rowSums(predData * predData)          # nPred
    cross  <- predData %*% t(libData)               # nPred x nLib
    d2     <- outer(predSq, libSq, "+") - 2 * cross # nPred x nLib
    d2[d2 < 0] <- 0                                 # guard FP negatives
    d  <- sqrt(d2)
    for (i in seq_len(nPred)) {
      o <- order(d[i, ])[seq_len(kEff)]            # stable; ties by index
      nnIdx[i,  seq_len(kEff)] <- o
      nnDist[i, seq_len(kEff)] <- d[i, o]
    }
  }

  list(nnIdx = nnIdx, nnDist = nnDist)
}

#------------------------------------------------------------------------
#' Row-wise cumulative sum of a logical / numeric matrix (vectorised).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
RowCumsum <- function(m) {
  k <- ncol(m)
  upper <- upper.tri(matrix(0, k, k), diag = TRUE) * 1   # k x k, U[a,b]=1 if a<=b
  (m + 0) %*% upper
}
