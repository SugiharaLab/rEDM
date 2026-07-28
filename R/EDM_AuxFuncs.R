#------------------------------------------------------------------------
# Shared marshalling utilities for the pure-R EDM kernels.
# Ported from rEDM/R/EDM_AuxFuncs.R; behaviour preserved.
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Test for a length-one atomic value.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

#------------------------------------------------------------------------
#' Assemble the Simplex/SMap prediction data.frame.
#'
#' Port of Formatting.py::FormatProjection for the projection/variance path
#' (Time, Observations, Predictions, Pred_Variance). Internal position
#' arithmetic mirrors pyEDM (0-based), converted to 1-based at R indexing.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
ProjectionLayout <- function(predList, pred_i, pred_i_all, timeVec,
                             Tp, nRows, tau) {
  TpMag <- abs(Tp)

  # 0-based locals to mirror pyEDM arithmetic
  predList0   <- lapply(predList, function(s) s - 1L)
  pred_i0     <- pred_i - 1L
  pred_i_all0 <- pred_i_all - 1L

  # --- Observations index vector obs_i (0-based) and total outSize ---
  obs_i   <- integer(0)
  outSize <- 0L
  for (seg in predList0) {
    nPred <- length(seg)
    outSize <- outSize + nPred + TpMag
    if (nPred == 0L) next
    if (Tp == 0) {
      appendI <- seg
    } else if (Tp > 0) {
      if (seg[nPred] + Tp < nRows) {
        appendI <- c(seg, (seg[nPred] + 1):(seg[nPred] + Tp))
      } else appendI <- seg
    } else {
      if (seg[1] + Tp > -1) {
        appendI <- c((seg[1] + Tp):(seg[1] - 1), seg)
      } else appendI <- seg
    }
    obs_i <- c(obs_i, appendI)
  }

  # --- Projection output positions predOut_i (0-based) ---
  predOut_i <- integer(0)
  predOut0  <- 0L
  for (seg in predList0) {
    nPred <- length(seg)
    if (nPred == 0L) next
    if (Tp == 0) {
      idx <- predOut0:(predOut0 + nPred - 1L)
      predOut_i <- c(predOut_i, idx); predOut0 <- idx[length(idx)] + 1L
    } else if (Tp > 0) {
      idx <- (predOut0 + Tp):(predOut0 + Tp + nPred - 1L)
      predOut_i <- c(predOut_i, idx); predOut0 <- idx[length(idx)] + 1L
    } else {
      idx <- predOut0:(predOut0 + nPred - 1L)
      predOut_i <- c(predOut_i, idx); predOut0 <- idx[length(idx)] + TpMag + 1L
    }
  }

  # NaN removal may have shortened pred_i vs predOut_i: remap
  if (length(pred_i0) != length(predOut_i)) {
    nMap <- length(predOut_i)
    if (tau < 0) keys <- pred_i_all0[(length(pred_i_all0) -
                                      nMap + 1L):length(pred_i_all0)]
    else         keys <- pred_i_all0[seq_len(nMap)]
    mp <- predOut_i
    names(mp) <- as.character(keys)
    predOut_i <- unname(mp[as.character(pred_i0)])
  }

  # --- Observation output positions obsOut_i (0-based) ---
  if (Tp > 0) {
    if (obs_i[length(obs_i)] + Tp > nRows - 1L) obsOut_i <- seq_len(length(obs_i)) - 1L
    else                                        obsOut_i <- seq_len(outSize) - 1L
  } else {
    obsOut_i <- seq_len(length(obs_i)) - 1L
    if (pred_i0[1] + Tp < 0) obsOut_i <- obsOut_i + (TpMag - pred_i0[1])
  }

  # --- Time ---
  timeOut <- rep(NA, outSize)
  if (Tp == 0 ||
      (Tp > 0 && (pred_i_all[length(pred_i_all)] + Tp) <= length(timeVec)) ||
      (Tp < 0 && (pred_i0[1] + Tp >= 0))) {
    timeOut[obsOut_i + 1L] <- timeVec[obs_i + 1L]
  } else {
    timeOut <- AddTime(timeVec, TpMag, outSize, obs_i, obsOut_i, Tp)
  }

  list(obs_i = obs_i, outSize = outSize, predOut_i = predOut_i,
       obsOut_i = obsOut_i, timeOut = timeOut)
}

#------------------------------------------------------------------------
#' Assemble the Simplex/SMap prediction data.frame.
#'
#' Port of Formatting.py::FormatProjection for the projection/variance path
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
FormatProjection <- function(predList, pred_i, pred_i_all, targetVec, timeVec,
                             projection, variance, Tp, nRows, tau, timeName) {

  L <- ProjectionLayout(predList, pred_i, pred_i_all, timeVec, Tp, nRows, tau)

  observationOut <- rep(NA_real_, L$outSize)
  projectionOut  <- rep(NA_real_, L$outSize)
  varianceOut    <- rep(NA_real_, L$outSize)
  observationOut[L$obsOut_i  + 1L] <- targetVec[L$obs_i + 1L]
  projectionOut [L$predOut_i + 1L] <- projection
  varianceOut   [L$predOut_i + 1L] <- variance

  dfNames = c( timeName, "Observations", "Predictions", "Pred_Variance" )
  df = data.frame(L$timeOut, observationOut, projectionOut, varianceOut,
                  stringsAsFactors = FALSE, check.names = FALSE)
  names(df) = dfNames
  df
}

#------------------------------------------------------------------------
#' Extend the time vector for forecast rows beyond the data (Tp shift).
#' Numeric time only; assumes a constant step from the last two values.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
AddTime <- function(timeVec, TpMag, outSize, obs_i, obsOut_i, Tp) {
  timeOut <- rep(NA, outSize)
  timeOut[obsOut_i + 1L] <- timeVec[obs_i + 1L]
  # Fill any trailing/leading NA by extrapolating a constant step
  known <- which(!is.na(timeOut))
  if (length(known) >= 2 && is.numeric(timeVec)) {
    step <- timeVec[2] - timeVec[1]
    for (i in seq_len(outSize)) {
      if (is.na(timeOut[i])) {
        earlier <- known[known < i]
        if (length(earlier)) {
          ref <- max(earlier)
          timeOut[i] <- timeOut[ref] + step * (i - ref)
        }
      }
    }
    # leading NA
    for (i in rev(seq_len(outSize))) {
      if (is.na(timeOut[i])) {
        later <- known[known > i]
        if (length(later)) {
          ref <- min(later)
          timeOut[i] <- timeOut[ref] - step * (ref - i)
        }
      }
    }
  }
  timeOut
}

#------------------------------------------------------------------------
#' Split a columns/target argument (space string or vector) into names.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
SplitColumns <- function(x) {
  if (length(x) == 1L && is.character(x) && grepl("\\s", x))
    strsplit(trimws(x), "\\s+")[[1]]
  else
    x
}

#------------------------------------------------------------------------
#' Build 1-based library and prediction indices from (start end) pairs.
#'
#' Port of EDM.py::CreateIndices with the 0-offset removed.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
CreateIndices <- function(lib, pred, E, tau, Tp, embedded, nRows,
                          name = "Simplex") {
  embedShift <- abs(tau) * (E - 1)

  libPairs <- matrix(ParsePairs(lib), ncol = 2, byrow = TRUE)
  libList  <- vector("list", nrow(libPairs))
  for (r in seq_len(nrow(libPairs))) {
    start <- libPairs[r, 1]; stop <- libPairs[r, 2]
    if (name %in% c("Simplex", "SMap", "Multiview") && start >= stop)
      stop(sprintf("CreateIndices: lib start %d exceeds lib end %d.", start, stop))
    if (start < 1 || stop < 1) stop("CreateIndices: lib indices < 1 not allowed.")
    if (!embedded) {
      if (tau < 0) start <- start + embedShift else stop <- stop - embedShift
    }
    if (Tp < 0) {
      if (!embedded) start <- max(start, start + abs(Tp) - 1)
    } else if (r == nrow(libPairs)) {
      stop <- stop - Tp
    }
    libList[[r]] <- start:stop
  }
  lib_i       <- unlist(libList)
  disjointLib <- length(libList) > 1

  predPairs <- matrix(ParsePairs(pred), ncol = 2, byrow = TRUE)
  predList  <- vector("list", nrow(predPairs))
  for (r in seq_len(nrow(predPairs))) {
    pStart <- predPairs[r, 1]; pStop <- predPairs[r, 2]
    if (name %in% c("Simplex", "SMap", "Multiview") && pStart >= pStop)
      stop(sprintf("CreateIndices: pred start %d exceeds pred end %d.", pStart, pStop))
    if (pStart < 1 || pStop < 1) stop("CreateIndices: pred indices < 1 not allowed.")
    predList[[r]] <- pStart:pStop
  }
  pred_i_all <- unlist(predList)

  # Remove leading/trailing embedding-NaN rows from predList (not embedded)
  if (!embedded && embedShift > 0) {
    nanStart <- seq_len(embedShift)                       # rows 1..embedShift
    nanEnd   <- (nRows - embedShift + 1):nRows            # last embedShift rows
    for (i in seq_along(predList)) {
      seg <- predList[[i]]
      if (tau > 0) predList[[i]] <- seg[!(seg %in% nanEnd)]
      else         predList[[i]] <- seg[!(seg %in% nanStart)]
    }
  }
  pred_i <- unlist(predList)

  if (lib_i[length(lib_i)] > nRows)
    stop("CreateIndices: library index exceeds number of data rows.")
  if (pred_i[length(pred_i)] > nRows)
    stop("CreateIndices: prediction index exceeds number of data rows.")

  libOverlap <- length(intersect(lib_i, pred_i)) > 0

  list(lib_i = lib_i, pred_i = pred_i, pred_i_all = pred_i_all,
       predList = predList, disjointLib = disjointLib,
       libOverlap = libOverlap, embedShift = embedShift)
}

#------------------------------------------------------------------------
#' Parse (start end) pairs from a space string or numeric vector.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
ParsePairs <- function(spec) {
  if (is.character(spec)) {
    v <- suppressWarnings(as.integer(strsplit(trimws(spec), "\\s+")[[1]]))
  } else v <- as.integer(spec)
  if (length(v) == 0L || anyNA(v) || length(v) %% 2L != 0L)
    stop("CreateIndices: lib/pred must be an even set of (start end) pairs.")
  v
}

#------------------------------------------------------------------------
#' Remove embedding-NaN rows from lib_i / pred_i (ignoreNan path).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
RemoveNan <- function(embedding, lib_i, pred_i, ignoreNan) {
  if (!ignoreNan) return(list(lib_i = lib_i, pred_i = pred_i))
  naRow <- function(idx) apply(embedding[idx, , drop = FALSE], 1,
                               function(r) any(is.na(r)))
  naLib  <- naRow(lib_i)
  naPred <- naRow(pred_i)
  list(lib_i  = lib_i[!naLib],
       pred_i = pred_i[!naPred])
}

#------------------------------------------------------------------------
#' Flatten a vector / list / data.frame / matrix to a delimited string.
#'
#' Preserves the rEDM convention used to pass \code{columns}, \code{target},
#' \code{lib} and \code{pred} to the kernels as space-delimited strings.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
FlattenToString <- function(x, delimiter = " ") {
  if (is.data.frame(x) || is.matrix(x)) {
    s <- ""
    for (row in seq_len(nrow(x))) {
      s <- paste(s, paste(x[row, ], collapse = delimiter),
                 collapse = delimiter)
    }
  } else if (is.list(x)) {
    s <- paste(unlist(x), collapse = delimiter)
  } else if (is.vector(x)) {
    s <- paste(x, collapse = delimiter)
  } else {
    s <- x
  }
  s
}

#------------------------------------------------------------------------
#' Parse a lib/pred specification into 1-based row indices.
#'
#' Accepts either a whitespace string ("1 100" or "1 100 201 300") or a
#' numeric vector (c(1, 100) or c(1, 100, 201, 300)). Values are interpreted
#' as inclusive (start, end) pairs and expanded to a 1-based integer index
#' vector. This is the R-native form: no 0-offset conversion is applied
#' (contrast pyEDM, which subtracts 1 for NumPy).
#'
#' @param spec character scalar or numeric vector of (start end) pairs.
#' @return integer vector of 1-based indices (sorted, de-duplicated within
#'   overlapping pairs preserved in order of pairs).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
ParseIndexPairs <- function(spec) {
  if (is.character(spec)) {
    toks <- strsplit(trimws(spec), "\\s+")[[1]]
    toks <- toks[nzchar(toks)]
    vals <- suppressWarnings(as.integer(toks))
  } else {
    vals <- as.integer(spec)
  }

  if (length(vals) == 0L || anyNA(vals)) {
    stop("ParseIndexPairs: spec must be integers as (start end) pairs, got: ",
         FlattenToString(spec))
  }
  if (length(vals) %% 2L != 0L) {
    stop("ParseIndexPairs: expected an even number of values (start end pairs), got ",
         length(vals), ": ", FlattenToString(vals))
  }

  idx <- integer(0)
  for (p in seq(1L, length(vals), by = 2L)) {
    start <- vals[p]
    end   <- vals[p + 1L]
    if (start < 1L) {
      stop("ParseIndexPairs: 1-based indices must be >= 1; got start = ", start)
    }
    if (end < start) {
      stop("ParseIndexPairs: pair end (", end, ") < start (", start, ")")
    }
    idx <- c(idx, start:end)
  }
  idx
}

#------------------------------------------------------------------------
#' Is the object a non-empty data.frame with numeric value columns?
#'
#' Port of rEDM isValidDataFrame: columns 2..ncol must be numeric (the first
#' column may be time / date / character).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
isValidDataFrame <- function(dataFrame) {
  if (!inherits(dataFrame, "data.frame")) {
    message("isValidDataFrame(): object is not a data.frame."); return(FALSE)
  }
  if (nrow(dataFrame) == 0 || ncol(dataFrame) == 0) {
    message("isValidDataFrame(): dataFrame is empty."); return(FALSE)
  }
  if (ncol(dataFrame) >= 2) {
    cls <- vapply(dataFrame[, 2:ncol(dataFrame), drop = FALSE],
                  function(col) class(col)[1], character(1))
    if (any(cls %in% c("character", "factor"))) {
      message("isValidDataFrame(): non-numeric value column detected."); return(FALSE)
    }
  }
  TRUE
}

#------------------------------------------------------------------------
#' Validate a data.frame and that columns/target exist in it.
#'
#' Port of rEDM ValidateDataFrame (dataFrame path). Returns TRUE/FALSE.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
ValidateDataFrame <- function(dataFrame, columns, target, verbose = FALSE) {
  if (!isValidDataFrame(dataFrame)) return(FALSE)
  nm <- names(dataFrame)
  for (c in SplitColumns(columns))
    if (!c %in% nm) { message("ValidateDataFrame(): column '", c, "' not found."); return(FALSE) }
  for (t in SplitColumns(target))
    if (!t %in% nm) { message("ValidateDataFrame(): target '", t, "' not found."); return(FALSE) }
  if (verbose) message("ValidateDataFrame(): dataFrame validated.")
  TRUE
}

#------------------------------------------------------------------------
#' Stop with a clear message when input is not usable (API entry guard).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
RequireValidData <- function(name, dataFrame, columns, target) {
  if (!ValidateDataFrame(dataFrame, columns, target))
    stop(sprintf("%s(): invalid dataFrame, or columns/target not found.", name),
         call. = FALSE)
  invisible(TRUE)
}
#------------------------------------------------------------------------
#' Generate surrogate time series for significance testing.
#'
#' Three methods:
#' \code{random_shuffle} (permutation), \code{ebisuzaki} (Fourier
#' phase-randomisation preserving the power spectrum / autocorrelation), and
#' \code{seasonal} (period-\code{T_period} spline cycle plus shuffled
#' residuals plus Gaussian noise).
#'
#' @param ts numeric time-series vector.
#' @param method one of "random_shuffle", "ebisuzaki", "seasonal".
#' @param num_surr number of surrogate series to generate.
#' @param T_period seasonal period (seasonal method).
#' @param alpha Gaussian-noise standard deviation (seasonal method).
#' @return numeric matrix with \code{length(ts)} rows and \code{num_surr} columns.
#' @keywords internal
#------------------------------------------------------------------------
SurrogateData <- function(ts,
                          method   = c("random_shuffle", "ebisuzaki", "seasonal"),
                          num_surr = 100, T_period = 1, alpha = 0) {
  method <- match.arg(method)

  if (method == "random_shuffle") {
    return(sapply(seq_len(num_surr), function(i) sample(ts, size = length(ts))))
  }

  if (method == "ebisuzaki") {
    if (any(!is.finite(ts)))
      stop("SurrogateData(): input time series contained invalid values")

    n  <- length(ts)
    n2 <- floor(n / 2)
    sigma <- stats::sd(ts)
    amplitudes    <- abs(stats::fft(ts))
    amplitudes[1] <- 0

    return(sapply(seq_len(num_surr), function(i) {
      if (n %% 2 == 0) {                                  # even length
        thetas    <- 2 * pi * stats::runif(n2 - 1)
        angles    <- c(0, thetas, 0, -rev(thetas))
        recf      <- amplitudes * exp(complex(imaginary = angles))
        recf[n2]  <- complex(real = sqrt(2) * amplitudes[n2] *
                             cos(stats::runif(1) * 2 * pi))
      } else {                                            # odd length
        thetas <- 2 * pi * stats::runif(n2)
        angles <- c(0, thetas, -rev(thetas))
        recf   <- amplitudes * exp(complex(imaginary = angles))
      }
      temp <- Re(stats::fft(recf, inverse = TRUE) / n)
      temp / stats::sd(temp) * sigma                      # match original variance
    }))
  }

  # seasonal
  if (any(!is.finite(ts)))
    stop("SurrogateData(): input time series contained invalid values")

  n        <- length(ts)
  I_season <- suppressWarnings(matrix(1:T_period, nrow = n, ncol = 1))
  seasonalF <- stats::smooth.spline(
    c(I_season - T_period, I_season, I_season + T_period), c(ts, ts, ts))
  seasonalCyc   <- stats::predict(seasonalF, I_season)$y
  seasonalResid <- ts - seasonalCyc

  sapply(seq_len(num_surr), function(i)
    seasonalCyc + sample(seasonalResid, n) + stats::rnorm(n, 0, alpha))
}
