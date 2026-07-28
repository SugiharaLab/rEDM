#------------------------------------------------------------------------
#' Takens time-delay embedding of one or more columns.
#'
#' @param dataFrame data.frame of source data.
#' @param E embedding dimension.
#' @param tau delay (negative = past lags, pyEDM default -1).
#' @param columns character vector (or space string) of column names.
#' @return numeric matrix (N x (E * length(columns))); NA in |tau|*(E-1) rows.
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
Embed <- function(dataFrame, E, tau, columns) {
  if (E < 1) stop("Embed: E must be positive.")
  if (tau == 0) stop("Embed: tau must be non-zero.")
  columns <- SplitColumns(columns)
  shiftVec <- seq(0, by = -tau, length.out = E)     # tau<0 -> 0,|tau|,2|tau|,...
  blocks <- vector("list", length(columns))
  for (ci in seq_along(columns)) {
    col <- dataFrame[[columns[ci]]]
    m <- matrix(NA_real_, length(col), E)
    for (j in seq_len(E)) m[, j] <- ShiftVec(col, shiftVec[j])
    blocks[[ci]] <- m
  }
  do.call(cbind, blocks)
}

#------------------------------------------------------------------------
#' Shift a vector by p (pandas .shift convention): p>0 moves values down
#' (NA at head), p<0 moves up (NA at tail).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
ShiftVec <- function(x, p) {
  n <- length(x)
  if (p == 0) return(x)
  if (p > 0)  return(c(rep(NA_real_, p), x[seq_len(n - p)]))
  p <- -p
  c(x[(p + 1):n], rep(NA_real_, p))
}
