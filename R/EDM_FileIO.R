#------------------------------------------------------------------------
# Legacy rEDM v1 compatibility: file-based I/O and parameterList.
#
# Self-contained helpers used by the public API functions. Adding the file
# parameters (pathIn, dataFile, pathOut, predictFile) and parameterList at
# their defaults leaves every function's numeric behaviour unchanged.
#
# Naming: functions/constants uppercase-first, variables lowercase-first.
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#' Resolve the input data.frame from a file or an in-memory frame.
#'
#' If \code{dataFile} is non-empty, read \code{file.path(pathIn, dataFile)}
#' (with \code{check.names = FALSE} so headers such as \code{x_t-1} or
#' \code{Var 5 1} survive); this takes precedence over \code{dataFrame}.
#' Otherwise return \code{dataFrame} unchanged.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
ResolveInput <- function(pathIn, dataFile, dataFrame) {
  if (nzchar(dataFile))
    return(utils::read.csv(file.path(pathIn, dataFile),
                           check.names = FALSE, stringsAsFactors = FALSE))
  dataFrame
}

#------------------------------------------------------------------------
#' Write the primary output frame and shape the return value.
#'
#' Optionally writes \code{primaryFrame} to \code{file.path(pathOut,
#' predictFile)}, then applies the \code{parameterList} contract:
#'   - frame result  -> list(predictions = result, parameters = parameters)
#'   - list result   -> result with a $parameters element appended
#' When \code{parameterList} is FALSE the result is returned unchanged.
#'
#' @param result the function's normal return object (data.frame or list).
#' @param primaryFrame the frame written to disk (result itself, or the
#'   predictions element of a list result).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
EdmFinalize <- function(result, primaryFrame, pathOut, predictFile,
                        parameterList, parameters,
                        includeState = FALSE, internal = NULL) {
  if (nzchar(predictFile))
    utils::write.csv(primaryFrame, file.path(pathOut, predictFile),
                     row.names = FALSE)
  if (isTRUE(parameterList) || isTRUE(includeState)) {
    if (is.data.frame(result))
      result <- list(predictions = result)
    if (isTRUE(parameterList)) result$parameters <- parameters
    if (isTRUE(includeState))  result$internal   <- internal
    return(result)
  }
  result
}

#------------------------------------------------------------------------
#' Resolve / require the embedding dimension E.
#'
#' pyEDM makes E a required argument for Simplex, SMap, PredictInterval and
#' PredictNonlinear. When embedded, E is inferred as length(columns) and any
#' user value is ignored; otherwise E must be supplied (non-NULL, >= 1).
#' @keywords internal
#------------------------------------------------------------------------
RequireE <- function(fn, E, embedded, columns) {
  if (isTRUE(embedded)) return(length(columns))
  if (is.null(E) || length(E) == 0L || E < 1L)
    stop(sprintf("%s(): E is required (integer >= 1).", fn), call. = FALSE)
  E
}
