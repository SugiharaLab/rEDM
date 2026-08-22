#------------------------------------------------------------------------
# Diagnostic parameter sweeps.
#
# Port of pyEDM EmbedDimension / PredictInterval / PredictExclusionRadius /
# PredictNonlinear (API.py + PoolFunc.py). Each sweeps one parameter over a
# list, runs Simplex (or SMap for PredictNonlinear) at each value, and
# collects ComputeError()$rho.
#
# The sweeps are RNG-free and order-preserving, so parallel and sequential
# evaluation give identical results. Column naming follows the rEDM
# convention (E, Tp, ExclusionRadius, Theta); pyEDM uses lower-case theta /
# 'Exclusion radius', which the validation harness reconciles.
#
# Naming: functions/constants uppercase-first, variables lowercase-first.
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' EmbedDimension : prediction skill vs embedding dimension E in 1..maxE.
#'
#' @inheritParams Simplex
#' @param maxE largest embedding dimension to evaluate.
#' @param numProcess worker processes for the sweep.
#' @return data.frame(E, rho).
#' @keywords internal
#------------------------------------------------------------------------
EmbedDimension <- function(dataFrame = NULL, columns, target, lib, pred,
                           maxE = 10, Tp = 1, tau = -1, exclusionRadius = 0,
                           embedded = FALSE, validLib = logical(0),
                           noTime = FALSE, ignoreNan = TRUE, numProcess = 4,
                           pathIn = "./", dataFile = "", pathOut = "./",
                           predictFile = "", backend = "RANN",
                           parameterList = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  Evals     <- seq_len(maxE)
  worker <- function(E) SweepRho(
    Simplex(dataFrame, columns, target, lib, pred, E = E, Tp = Tp, tau = tau,
            exclusionRadius = exclusionRadius, embedded = embedded,
            validLib = validLib, noTime = noTime, ignoreNan = ignoreNan,
            backend = backend))

  rho <- unlist(SweepApply(Evals, worker, numProcess))
  out <- data.frame(E = Evals, rho = rho)

  if (showPlot) {
    title = paste( columns, ":", target, "\nTp=", Tp)
    PlotSweep(out, "E", title)
  }
  out = EdmFinalize( out, out, pathOut, predictFile, parameterList, parameters )
}

#------------------------------------------------------------------------
#' PredictInterval : prediction skill vs forecast interval Tp in 1..maxTp.
#'
#' @inheritParams Simplex
#' @param maxTp largest forecast interval to evaluate.
#' @return data.frame(Tp, rho).
#' @keywords internal
#------------------------------------------------------------------------
PredictInterval <- function(dataFrame = NULL, columns, target, lib, pred,
                            maxTp = 10, E = NULL, tau = -1, exclusionRadius = 0,
                            embedded = FALSE, validLib = logical(0),
                            noTime = FALSE, ignoreNan = TRUE, numProcess = 4,
                            backend = "RANN", pathIn = "./", dataFile = "",
                            pathOut = "./", predictFile = "", 
                            parameterList = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  E         <- RequireE("PredictInterval", E, embedded, SplitColumns(columns))
  Evals     <- seq_len(maxTp)

  worker <- function(Tp) SweepRho(
    Simplex(dataFrame, columns, target, lib, pred, E = E, Tp = Tp, tau = tau,
            exclusionRadius = exclusionRadius, embedded = embedded,
            validLib = validLib, noTime = noTime, ignoreNan = ignoreNan,
            backend = backend))

  rho <- unlist(SweepApply(Evals, worker, numProcess))
  out <- data.frame(Tp = Evals, rho = rho)

  if (showPlot) {
    title = paste( columns, ":", target, "\nE=", E) 
    PlotSweep(out, "Tp", title)
  }

  out = EdmFinalize( out, out, pathOut, predictFile, parameterList, parameters )
}

#------------------------------------------------------------------------
#' PredictExclusionRadius : prediction skill vs Theiler exclusion radius.
#'
#' @inheritParams Simplex
#' @param exclusionRadius numeric vector of radii (default pyEDM sweep list).
#' @return data.frame(ExclusionRadius, rho).
#' @keywords internal
#------------------------------------------------------------------------
PredictExclusionRadius <- function(dataFrame = NULL, columns, target, lib, pred,
                                   exclusionRadius = NULL, E = 1, Tp = 1,
                                   tau = -1, embedded = FALSE,
                                   validLib = logical(0), noTime = FALSE,
                                   ignoreNan = TRUE, numProcess = 4,
                                   backend = "RANN", pathIn = "./", 
                                   dataFile = "", pathOut = "./",
                                   predictFile = "", 
                                   parameterList = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  if (is.null(exclusionRadius)) {
    exclusionRadius <- c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30)
  }
  worker <- function(er) SweepRho(
    Simplex(dataFrame, columns, target, lib, pred, E = E, Tp = Tp, tau = tau,
            exclusionRadius = er, embedded = embedded, validLib = validLib,
            noTime = noTime, ignoreNan = ignoreNan, backend = backend))

  rho <- unlist(SweepApply(exclusionRadius, worker, numProcess))
  out <- data.frame(ExclusionRadius = exclusionRadius, rho = rho)
    
  if (showPlot) {
      title = paste( columns, ":", target, "\nE=", E, ' Tp=', Tp)
      PlotSweep(out, "Exclusion Radius", title)
  }
  out = EdmFinalize( out, out, pathOut, predictFile, parameterList, parameters )
}

#------------------------------------------------------------------------
#' PredictNonlinear : S-map prediction skill vs localisation theta.
#'
#' @inheritParams SMap
#' @param theta numeric vector of theta values (default pyEDM sweep list).
#' @return data.frame(Theta, rho).
#' @keywords internal
#------------------------------------------------------------------------
PredictNonlinear <- function(dataFrame = NULL, columns, target, lib, pred,
                             theta = NULL, E = NULL, Tp = 1, knn = 0, tau = -1,
                             exclusionRadius = 0, embedded = FALSE,
                             validLib = logical(0), noTime = FALSE,
                             ignoreNan = TRUE, numProcess = 4,
                             backend = "RANN", pathIn = "./", dataFile = "",
                             pathOut = "./", predictFile = "",
                             parameterList = FALSE, showPlot = FALSE,
                             verbose = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  E         <- RequireE("PredictNonlinear", E, embedded, SplitColumns(columns))
    
  if (is.null(theta)) {
    theta <- c(0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9)
  }

  worker <- function(th) SweepRho(
    SMap(dataFrame, columns, target, lib, pred, E = E, Tp = Tp, knn = knn,
         tau = tau, theta = th, exclusionRadius = exclusionRadius,
         embedded = embedded, validLib = validLib, noTime = noTime,
         ignoreNan = ignoreNan, backend = backend)$predictions)

  rho <- unlist(SweepApply(theta, worker, numProcess))
  out <- data.frame(Theta = theta, rho = rho)

  if (showPlot) {
    title = paste( columns, ":", target, "\nE=", E, ' Tp=', Tp) 
    PlotSweep(out, expression(paste("S-map Localization",~(theta))), title)
  }
  out = EdmFinalize( out, out, pathOut, predictFile, parameterList, parameters )
}

#------------------------------------------------------------------------
#' Map a worker over a sweep list, in parallel where supported.
#'
#' Uses forked parallel::mclapply on unix when numProcess > 1; falls back to
#' sequential lapply otherwise. Results are deterministic either way.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
SweepApply <- function(values, worker, numProcess = 4) {
  useParallel <- numProcess>1 && .Platform$OS.type == "unix" && length(values)>1
  if (useParallel) {
    # Respect the R CMD check core limit (_R_CHECK_LIMIT_CORES_ caps at 2);
    # otherwise use the machine core count.
    limited  <- isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "false")))
    maxCores <- if (limited) 2L else parallel::detectCores()
      
    if (is.na(maxCores) || maxCores < 1L) { maxCores <- 1L }

    cores <- max(1L, min(numProcess, length(values), maxCores))

    if (cores > 1L) {
      out <- tryCatch(
        parallel::mclapply(values, worker, mc.cores = cores,
                           mc.preschedule = TRUE),
        error = function(e) NULL)
      
      if (!is.null(out) &&
          !any(vapply(out, function(x) inherits(x, "try-error"), logical(1)))) {
        return(out)
      }
    }
  }
  lapply(values, worker) # sequential fallback
}

#------------------------------------------------------------------------
#' rho from a Simplex/SMap prediction frame.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
SweepRho <- function(predictions) {
  ComputeError(predictions$Observations, predictions$Predictions)$rho
}
