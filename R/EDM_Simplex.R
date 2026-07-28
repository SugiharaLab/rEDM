#------------------------------------------------------------------------
# Simplex projection
#
# Port of pyEDM Simplex (Simplex.py, EDM.py, Formatting.py, API.Embed).
# See spec: rEDM_Simplex_spec.md.
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Simplex projection (top-level API).
#'
#' @param dataFrame data.frame; first column is time unless noTime = TRUE.
#' @param columns,target column name(s) (vector or space string).
#' @param lib,pred (start end) index pairs (vector or space string).
#' @param E,Tp,knn,tau,exclusionRadius standard EDM parameters.
#' @param embedded TRUE if columns already form the embedding.
#' @param validLib logical mask of admissible library rows (or length 0).
#' @param noTime TRUE to synthesise a 1..N time index.
#' @param ignoreNan remove embedding-NaN rows from lib/pred.
#' @param backend "RANN" or "brute".
#' @param pathIn file path for input dataFile.
#' @param dataFile input dataFile, .csv format.
#' @param pathOut file path for predictFile.
#' @param predictFile output file name.
#' @param parameterList append list of parameters to return.
#' @return data.frame(Time, Observations, Predictions, Pred_Variance).
#' @keywords internal
#------------------------------------------------------------------------
Simplex <- function(dataFrame = NULL, columns, target, lib, pred,
                    E = 0, Tp = 1, knn = 0, tau = -1, exclusionRadius = 0,
                    embedded = FALSE, validLib = logical(0), noTime = FALSE,
                    ignoreNan = TRUE, backend = "RANN", pathIn = "./",
                    dataFile = "", pathOut = "./", predictFile = "",
                    parameterList = FALSE, verbose = FALSE, showPlot = FALSE) {

  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  columns   <- SplitColumns(columns)
  target    <- SplitColumns(target)
  RequireValidData("Simplex", dataFrame, columns, target)
  nRows     <- nrow(dataFrame)


  if (embedded) E  <- length(columns)  # Any user supplied E ignored
  if (knn < 1) knn <- E + 1            # Simplex default

  idx <- CreateIndices(lib, pred, E, tau, Tp, embedded, nRows, "Simplex")

  targetVec <- dataFrame[[target[1]]]
  timeVec   <- if (noTime) seq_len(nRows) else dataFrame[[1]]

  embedding <- if (embedded) as.matrix(dataFrame[, columns, drop = FALSE])
               else Embed(dataFrame, E, tau, columns)

  rn        <- RemoveNan(embedding, idx$lib_i, idx$pred_i, ignoreNan)
  lib_i     <- rn$lib_i
  pred_i    <- rn$pred_i

  nb <- FindNeighbors(embedding, lib_i, pred_i, knn,
                      exclusionRadius = exclusionRadius, validLib = validLib,
                      libOverlap = idx$libOverlap, backend = backend,
                      verbose = verbose)

  pr <- SimplexProject(nb$neighbors, nb$distances, targetVec, Tp)

  timeName = "Time"
  if (!noTime && !is.null(names(dataFrame))) timeName <- names(dataFrame)[1]

  out <- FormatProjection(idx$predList, pred_i, idx$pred_i_all, targetVec,
                          timeVec, pr$projection, pr$variance, Tp, nRows,
                          tau, timeName)
  if (showPlot) PlotObsPred(out, "", E, Tp)

  out = EdmFinalize(out, out, pathOut, predictFile, parameterList, parameters)
}

#------------------------------------------------------------------------
#' Simplex projection kernel (Simplex.py::Project).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
SimplexProject <- function(knnNeighbors, knnDistances, targetVec, Tp) {
  minDistances    <- pmax(knnDistances[, 1], 1e-6)
  scaledDistances <- knnDistances / minDistances
  weights         <- exp(-scaledDistances)
  weightRowSum    <- rowSums(weights)

  neighborsTp <- knnNeighbors + Tp
  n <- length(targetVec)
  libTargets  <- matrix(NA_real_, nrow(knnNeighbors), ncol(knnNeighbors))
  inRange     <- neighborsTp >= 1 & neighborsTp <= n & knnNeighbors >= 1
  libTargets[inRange] <- targetVec[neighborsTp[inRange]]

  projection <- rowSums(weights * libTargets) / weightRowSum
  delta      <- libTargets - projection
  variance   <- rowSums(weights * delta^2) / weightRowSum
  list(projection = projection, variance = variance)
}
