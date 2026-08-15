#------------------------------------------------------------------------
# Multiview embedding (Ye & Sugihara 2016).
#
# Port of pyEDM Multiview (Multiview.py, API.Multiview, PoolFunc). Build the
# multivariate delay embedding, form all D-column combinations ("views"),
# rank each view by in-/out-of-sample Simplex rho, average the predictions of
# the top sqrt(nCombo) views, and report the ranked view table.
#
# Built entirely on the validated Simplex path (embedded, noTime), so results
# follow Simplex exactly; only view selection and averaging are new.
#
# Naming: functions/constants uppercase-first, variables lowercase-first.
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Multiview embedding forecast.
#'
#' @param dataFrame data.frame of source data.
#' @param columns library variable name(s) (vector or space string).
#' @param target target variable name.
#' @param lib,pred (start end) index pairs.
#' @param D view dimension (number of embedding columns combined; default
#'   = number of input columns).
#' @param E per-variable embedding dimension.
#' @param Tp,knn,tau,exclusionRadius standard EDM parameters.
#' @param multiview number of top views to average (default floor(sqrt(nCombo))).
#' @param trainLib rank on the library (in-sample) rather than pred.
#' @param excludeTarget drop the target variable from candidate views.
#' @param numProcess worker processes for the ranking grid.
#' @param backend "RANN" or "brute".
#' @return list(Predictions = data.frame(Time, Observations, Predictions),
#'   View = data.frame(Columns, rho, MAE, CAE, RMSE)).
#' @keywords internal
#------------------------------------------------------------------------
Multiview <- function(dataFrame = NULL, columns, target, lib = numeric(0),
                      pred = numeric(0), D = 0, E = 1, Tp = 1, knn = 0,
                      tau = -1, multiview = 0, exclusionRadius = 0,
                      trainLib = TRUE, excludeTarget = FALSE, ignoreNan = TRUE,
                      numProcess = 4, backend = "RANN",
                      pathIn = "./", dataFile = "", pathOut = "./", predictFile = "",
                      parameterList = FALSE, showPlot = FALSE) {
    
  parameters <- as.list(environment())
  parameters[c("dataFrame")] <- NULL

  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  columns   <- SplitColumns(columns)
  target    <- SplitColumns(target)
  target1   <- paste0(target[1], "(t-0)")             # embedding-column label
  N         <- nrow(dataFrame)

  if (trainLib && !length(lib) && !length(pred)) {
    lib  <- c(1, floor(N / 2))
    pred <- c(floor(N / 2) + 1, N)
  }
  if (D == 0) D <- length(columns)
  comboCols <- if (excludeTarget) setdiff(columns, target) else columns

  # Multivariate delay embedding as a named, time-less data.frame
  emb      <- MakeBlock(dataFrame, E, tau, comboCols)
  embNames <- EmbedColumnNames(comboCols, E, tau, embedded = FALSE)
  embDf    <- as.data.frame(emb, stringsAsFactors = FALSE)
  names(embDf) <- embNames

  combos <- utils::combn(embNames, D, simplify = FALSE)
  if (multiview < 1)              multiview <- floor(sqrt(length(combos)))
  if (multiview > length(combos)) multiview <- length(combos)

  RunView <- function(combo, predRange)
    Simplex(embDf, columns = combo, target = target1, lib = lib,
            pred = predRange, E = D, Tp = Tp, tau = tau,
            exclusionRadius = exclusionRadius, embedded = TRUE,
            noTime = TRUE, ignoreNan = ignoreNan, backend = backend,
            .tieBreak = FALSE)   # Multiview ranking keeps 2.0.2 neighbour order

  # ---- Rank views by Simplex rho ----
  rankPred <- if (trainLib) lib else pred
  rhoVec <- unlist(SweepApply(combos, function(combo) {
    df <- RunView(combo, rankPred)
    ComputeError(df$Observations, df$Predictions)$rho
  }, numProcess))

  # rev(order()) matches numpy argsort()[::-1]: NA last ascending -> first here
  topRank   <- rev(order(rhoVec))[seq_len(multiview)]
  topCombos <- combos[topRank]

  # ---- Project top views and average ----
  projList <- SweepApply(topCombos, function(combo) RunView(combo, pred),
                         numProcess)

  predMat <- vapply(projList, function(d) d$Predictions,
                    numeric(nrow(projList[[1]])))
  if (is.null(dim(predMat))) predMat <- matrix(predMat, ncol = length(projList))
  multiviewPredict <- rowMeans(predMat)

  first <- projList[[1]]
  predictions <- data.frame(Time = first$Time,
                            Observations = first$Observations,
                            Predictions = multiviewPredict,
                            stringsAsFactors = FALSE, check.names = FALSE)

  # ---- View ranking table ----
  stats <- lapply(projList, function(d) ComputeError(d$Observations, d$Predictions))
  view  <- data.frame(
    Columns = vapply(topCombos,
                     function(c) paste0("(",paste0("'", c, "'", collapse = ", "), ")"),
                     character(1)),
    rho  = vapply(stats, `[[`, numeric(1), "rho"),
    MAE  = vapply(stats, `[[`, numeric(1), "MAE"),
    CAE  = vapply(stats, `[[`, numeric(1), "CAE"),
    RMSE = vapply(stats, `[[`, numeric(1), "RMSE"),
    stringsAsFactors = FALSE, check.names = FALSE)

  if (showPlot) PlotObsPred(predictions, "Multiview  ", D, Tp)
  if (showPlot) PlotObsPred(predictions, "", D, Tp)
  result = list(Predictions = predictions, View = view)
  result = EdmFinalize( result, predictions, pathOut, predictFile,
                        parameterList, parameters )
}
