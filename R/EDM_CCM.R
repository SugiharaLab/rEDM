#------------------------------------------------------------------------
# Convergent Cross Mapping (CCM), pyEDM 2.5.0 model.
#
# Port of pyEDM CCM.py. Strategy: precompute delay embeddings and validity
# filtering once; for each (direction x libSize) task, run `sample` random
# subsample trials. Each trial builds a KDTree from L library vectors,
# queries all M valid points, computes CCM Simplex-weighted cross-map
# predictions, and takes the Pearson correlation.
#
# CCM has its OWN Simplex weighting (NOT the Simplex kernel): no 1e-6 floor,
# d_min_nz = where(d>0, d, 1), explicit zero-distance branch, explicit
# normalization.
#
# RNG: NumPy's PCG64 / SeedSequence subsampling cannot be bit-reproduced in
# R, so results cannot match pyEDM's fixtures exactly (see the tiered
# validation in the tests). For reproducibility, each (direction x libSize)
# task draws from its own L'Ecuyer-CMRG stream, making results independent of
# execution order and of sequential-vs-parallel dispatch.
#
#------------------------------------------------------------------------


#------------------------------------------------------------------------
#' Convergent Cross Mapping (top-level API).
#'
#' @param dataFrame data.frame of source data.
#' @param columns,target variable name(s); univariate targets use the first.
#' @param E,Tp,knn,tau,exclusionRadius standard EDM parameters (CCM Tp default 0).
#' @param libSizes explicit vector/string or (start stop increment).
#' @param sample subsamples per library size.
#' @param seed RNG seed (reproducible; not bit-comparable to pyEDM/NumPy).
#' @param embedded TRUE if columns/target already form embeddings.
#' @param validLib logical mask of admissible rows (or length 0).
#' @param includeData add per-libSize sample variance columns.
#' @param numProcess worker processes for the (direction x libSize) grid.
#' @param backend "RANN" or "brute".
#' @param pathIn file path for input dataFile.
#' @param dataFile input dataFile, .csv format.
#' @param pathOut file path for predictFile.
#' @param predictFile output file name.
#' @return data.frame(LibSize, <col:tgt>, <tgt:col> [, _var columns]).
#' @keywords internal
#------------------------------------------------------------------------
CCM <- function(dataFrame = NULL, columns, target, E, Tp = 0, knn = 0, tau = -1,
                exclusionRadius = 0, libSizes, sample = 30, seed = NULL,
                embedded = FALSE, validLib = logical(0), includeData = FALSE,
                numProcess = 4, backend = "RANN",
                pathIn = "./", dataFile = "", pathOut = "./", predictFile = "",
                parameterList = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  columns   <- SplitColumns(columns)
  target    <- SplitColumns(target)
  RequireValidData("CCM", dataFrame, columns, target)
  nRows     <- nrow(dataFrame)

  if (embedded) { E   <- length(columns) }
  if (knn < 1)  { knn <- E + 1L }

  libSizes <- CcmLibSizes(libSizes, E, nRows)
  shifts   <- (0:(E - 1)) * tau

  tgtVec     <- as.numeric(dataFrame[[target[1]]])
  colPredVec <- as.numeric(dataFrame[[columns[1]]])
  N          <- nRows

  if (embedded) {
    embCol   <- as.matrix(dataFrame[, columns, drop = FALSE])
    embTgt   <- as.matrix(dataFrame[, target,  drop = FALSE])
    validCol <- !apply(embCol, 1, function(r) any(is.na(r)))
    validTgt <- !apply(embTgt, 1, function(r) any(is.na(r)))
  }
  else {
    colVecs <- lapply(columns, function(c) as.numeric(dataFrame[[c]]))
    tgtVecs <- lapply(target,  function(c) as.numeric(dataFrame[[c]]))
    bc      <- CcmBuildEmbedding(colVecs, shifts)
    embCol  <- bc$emb; validCol <- bc$valid
    bt      <- CcmBuildEmbedding(tgtVecs, shifts)
    embTgt  <- bt$emb; validTgt <- bt$valid
  }

  if (length(validLib)) {
    validCol <- validCol & validLib
    validTgt <- validTgt & validLib
  }
  validCol <- validCol & CcmTpValidMask(N, tgtVec, Tp)     & !is.na(colPredVec)
  validTgt <- validTgt & CcmTpValidMask(N, colPredVec, Tp) & !is.na(tgtVec)

  idxCol <- which(validCol)
  idxTgt <- which(validTgt)

  predValsCol   <- tgtVec[idxCol + Tp]
  predValsTgt   <- colPredVec[idxTgt + Tp]
  embedValidCol <- embCol[idxCol, , drop = FALSE]
  embedValidTgt <- embTgt[idxTgt, , drop = FALSE]

  # ---- Task grid: 2 x nLib, each with its own L'Ecuyer stream ----
  nLib    <- length(libSizes)
  nTasks  <- 2L * nLib
  oldKind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(oldKind[1]), add = TRUE)
  if (is.null(seed)) seed <- as.integer(stats::runif(1, 1, .Machine$integer.max))
  set.seed(seed)
  streams <- vector("list", nTasks)
  cur <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  for (i in seq_len(nTasks)) { streams[[i]] <- cur;
                               cur <- parallel::nextRNGStream(cur) }

  tasks <- vector("list", nTasks)
  for (i in seq_len(nLib)) {
      tasks[[2L * i - 1L]] <- list(dir = "fwd", L = libSizes[i],
                                   stream = streams[[2L * i - 1L]])
      tasks[[2L * i]]      <- list(dir = "rev",
                                   L = libSizes[i], stream = streams[[2L * i]])
  }

  worker <- function(task) {
    assign(".Random.seed", task$stream, envir = .GlobalEnv)
    if (task$dir == "fwd") {
      rhos <- CcmForLibSize(embedValidCol, predValsCol, idxCol, task$L,
                            sample, knn, exclusionRadius, backend)
    }
    else {
      rhos <- CcmForLibSize(embedValidTgt, predValsTgt, idxTgt, task$L,
                            sample, knn, exclusionRadius, backend)
    }
    
    c(mean = mean(rhos, na.rm = TRUE),
      var  = if (includeData) {
                stats::var(rhos[!is.na(rhos)]) * (sum(!is.na(rhos))-1) /
                sum(!is.na(rhos)) # population var
             }
             else {
               NA_real_
             }
      ) 
  }

  results <- SweepApply(tasks, worker, numProcess)

  fwdMean <- vapply(seq_len(nLib),
                    function(i) results[[2L * i - 1L]]["mean"], numeric(1))
  revMean <- vapply(seq_len(nLib),
                    function(i) results[[2L * i]]["mean"],      numeric(1))

  cn <- paste0(columns[1], ":", target[1])
  tn <- paste0(target[1],  ":", columns[1])
  out <- data.frame(LibSize = libSizes, a = fwdMean, b = revMean,
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c("LibSize", cn, tn)

  if (includeData) {
    out[[paste0(cn, "_var")]] <- vapply(seq_len(nLib),
                                        function(i) results[[2L*i - 1L]]["var"],
                                        numeric(1))
    out[[paste0(tn, "_var")]] <- vapply(seq_len(nLib),
                                        function(i) results[[2L * i]]["var"],
                                        numeric(1))
  }
  if (showPlot) { PlotCCM(out, E) }
  out = EdmFinalize(out, out, pathOut, predictFile, parameterList, parameters)
}

#------------------------------------------------------------------------
#' Single-variable delay embedding (CCM convention; NA in invalid rows).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CcmSingleEmbed <- function(vec, shifts) {
  N <- length(vec); E <- length(shifts)
  emb <- matrix(NA_real_, N, E)
  for (dim in seq_len(E)) {
    s <- shifts[dim]
    if (s <= 0) { emb[(1 - s):N, dim] <- vec[1:(N + s)] }
    else        { emb[1:(N - s), dim] <- vec[(s + 1):N] }
  }
  emb
}

#------------------------------------------------------------------------
#' Mixed-multivariate delay embedding: horizontal stack of per-variable
#' embeddings. Returns emb and the AND of per-variable validity masks.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CcmBuildEmbedding <- function(vectors, shifts) {
  blocks <- lapply(vectors, CcmSingleEmbed, shifts = shifts)
  emb    <- do.call(cbind, blocks)
  valid  <- !apply(emb, 1, function(r) any(is.na(r)))
  list(emb = emb, valid = valid)
}

#------------------------------------------------------------------------
#' True where t+Tp is in bounds and vec[t+Tp] is finite.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CcmTpValidMask <- function(N, vec, Tp) {
  shifted  <- seq_len(N) + Tp
  inBounds <- shifted >= 1 & shifted <= N
  clipped  <- pmin(pmax(shifted, 1), N)
  inBounds & !is.na(vec[clipped])
}

#------------------------------------------------------------------------
#' Pearson r over finite pairs; NA if fewer than 3 valid pairs.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
NanSafePearson <- function(pred, act) {
  ok <- !(is.na(pred) | is.na(act))

  if (sum(ok) < 3) { return(NA_real_) }

  p <- pred[ok]; a <- act[ok]
  pm <- p - mean(p); am <- a - mean(a)
  denom <- sqrt(sum(pm^2) * sum(am^2))

  if (denom == 0) { return(0) }

  sum(pm * am) / denom
}

#------------------------------------------------------------------------
#' KNN query of all points against an L-point library (unsquared Euclidean).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CcmKnnQuery <- function(libEmbed, queryEmbed, kQuery, backend) {
  if (backend == "RANN") {
    r <- RANN::nn2(libEmbed, queryEmbed, k = kQuery,
                   treetype = "kd", searchtype = "standard", eps = 0)
    return(list(idx = r$nn.idx, dist = r$nn.dists))
  }
  M    <- nrow(queryEmbed)
  idx  <- matrix(0L, M, kQuery)
  dist <- matrix(0, M, kQuery)
  libT <- t(libEmbed)
  for (i in seq_len(M)) {
    d <- sqrt(colSums((libT - queryEmbed[i, ])^2))
    o <- order(d)[seq_len(kQuery)]
    idx[i, ]  <- o
    dist[i, ] <- d[o]
  }
  list(idx = idx, dist = dist)
}

#------------------------------------------------------------------------
#' Run S subsample CCM trials at one library size (CCM.py::_ccm_for_libsize).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CcmForLibSize <- function(embed, predVals, timeIndices, L, S, k,
                          exclusionRadius, backend) {
  M <- nrow(embed)
  L <- min(L, M)
  k <- min(k, L - 1L)
  if (k < 1L) return(rep(NA_real_, S))

  kQuery <- k + 1L
  if (exclusionRadius > 0) { kQuery <- kQuery + 2L * exclusionRadius }
  kQuery <- min(kQuery, L)

  rhos <- numeric(S)
  rowM <- seq_len(M)

  for (s in seq_len(S)) {
    libIdx   <- sort(sample.int(M, L))               # replace = FALSE
    libEmbed <- embed[libIdx, , drop = FALSE]
    q        <- CcmKnnQuery(libEmbed, embed, kQuery, backend)
    nnGlobal <- matrix(libIdx[q$idx], nrow = M)      # local -> M-space
    nnDist   <- q$dist

    isSelf <- nnGlobal == rowM                       # column-recycled compare
    if (exclusionRadius > 0) {
      nnT    <- matrix(timeIndices[nnGlobal], nrow = M)
      isExcl <- abs(timeIndices[rowM] - nnT) <= exclusionRadius
      mask   <- isSelf | isExcl
    }
    else {
      mask <- isSelf
    }

    valid <- !mask
    cs    <- RowCumsum(valid)
    firstK       <- valid & (cs <= k)
    totalFound   <- cs[, ncol(cs)]
    insufficient <- totalFound < k

    if (all(insufficient)) { rhos[s] <- NA_real_; next }

    firstK[insufficient, ] <- FALSE

    suf  <- which(!insufficient)
    sel  <- t(vapply(suf, function(r) which(firstK[r, ])[seq_len(k)],
                     integer(k)))                    # nSuf x k column positions
    if (k == 1L) sel <- matrix(sel, ncol = 1L)

    rr   <- matrix(suf, length(suf), k)
    Dsel <- matrix(nnDist[cbind(as.vector(rr), as.vector(sel))],
                   nrow = length(suf))
    Gsel <- matrix(nnGlobal[cbind(as.vector(rr), as.vector(sel))],
                   nrow = length(suf))

    # CCM Simplex weighting (no 1e-6 floor)
    dMin   <- Dsel[, 1]
    dMinNz <- ifelse(dMin > 0, dMin, 1)
    W      <- exp(-Dsel / dMinNz)
    zero   <- dMin == 0
    if (any(zero)) {
      W[zero, ] <- (Dsel[zero, , drop = FALSE] == 0) * 1
    }
    wSum <- rowSums(W); wSum <- ifelse(wSum > 0, wSum, 1)
    W    <- W / wSum

    predictions      <- rep(NA_real_, M)
    predictions[suf] <- rowSums(W * matrix(predVals[Gsel], nrow = length(suf)))
    rhos[s] <- NanSafePearson(predictions, predVals)
  }
  rhos
}

#------------------------------------------------------------------------
#' Expand a libSizes argument (explicit list or start/stop/increment).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CcmLibSizes <- function(libSizes, E, nRows) {
  if (is.character(libSizes)) {
    libSizes <- as.integer(strsplit(trimws(libSizes), "\\s+")[[1]])
  }
  libSizes <- as.integer(libSizes)
  
  if (length(libSizes) == 3L) {
    start <- libSizes[1]; stop <- libSizes[2]; incr <- libSizes[3]
    if (incr < stop) libSizes <- seq.int(start, stop, by = incr)
  }

  if (libSizes[length(libSizes)] > nRows) {
    stop("CCM: maximum libSize exceeds data size.")
  }

  libSizes
}
