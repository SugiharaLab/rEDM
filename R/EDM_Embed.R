#------------------------------------------------------------------------
#' Takens time-delay embedding block (numeric matrix).
#'
#' Internal engine builder: returns the unnamed numeric embedding matrix
#' consumed by the neighbour search and projection kernels. All rows are
#' retained; the first |tau|*(E-1) partial rows contain NA.
#'
#' @param dataFrame data.frame of source data.
#' @param E embedding dimension.
#' @param tau delay (negative = past lags, pyEDM default -1).
#' @param columns character vector (or space string) of column names.
#' @return numeric matrix (N x (E * length(columns))).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
MakeBlock <- function(dataFrame, E, tau, columns) {
  if (E < 1)    stop("MakeBlock(): E must be positive.")
  if (tau == 0) stop("MakeBlock(): tau must be non-zero.")
  columns <- SplitColumns(columns)
  if (length(columns) == 0L) stop("MakeBlock(): columns required.")
  absent <- setdiff(columns, names(dataFrame))
  if (length(absent))
    stop(sprintf("MakeBlock(): column(s) not found: %s",
                 paste(absent, collapse = ", ")))

  shiftVec <- seq(0, by = -tau, length.out = E)   # tau<0 -> 0,|tau|,2|tau|,...
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
#' Takens time-delay embedding of one or more columns.
#'
#' Wrapper over \code{MakeBlock} that returns a named data.frame matching
#' pyEDM's \code{Embed()}: each embedding column is named
#' \code{column(t-lag)} (or \code{(t+lag)} for positive \code{tau}), where
#' \code{lag = i * |tau|} for \code{i = 0..E-1}. All rows are kept, with NA
#' in the first \code{|tau|*(E-1)} partial rows.
#'
#' @param dataFrame data.frame of source data; column 1 is time.
#' @param E embedding dimension.
#' @param tau delay (negative = past lags, pyEDM default -1).
#' @param columns column name(s): a character vector or a space-separated string.
#' @param includeTime if TRUE, prepend column 1 of \code{dataFrame} (time),
#'   keeping its original name.
#' @param pathIn directory for \code{dataFile} (used when reading from disk).
#' @param dataFile CSV file name to read as input; when set it takes
#'   precedence over \code{dataFrame}.
#' @return a data.frame of the embedding (with the time column first when
#'   \code{includeTime = TRUE}).
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
Embed <- function(dataFrame = NULL, E = 0, tau = -1, columns = "",
                  includeTime = FALSE, pathIn = "./", dataFile = "") {
  dataFrame <- ResolveInput(pathIn, dataFile, dataFrame)
  columns <- SplitColumns(columns)
  block   <- MakeBlock(dataFrame, E, tau, columns)   # validates E / tau / columns
  df      <- as.data.frame(block, stringsAsFactors = FALSE)
  names(df) <- EmbedColumnNames(columns, E, tau, embedded = FALSE)
  if (isTRUE(includeTime))
    df <- cbind(dataFrame[, 1, drop = FALSE], df)
  df
}

#------------------------------------------------------------------------
#' Embedding column names, pyEDM convention.
#'
#' For column \code{c} and step \code{i = 0..E-1} the lag is \code{i*|tau|};
#' the name is \code{c(t-lag)} when \code{tau < 0} and \code{c(t+lag)} when
#' \code{tau > 0}. When \code{embedded}, \code{columns} already name the
#' embedding and are returned unchanged. Shared by \code{Embed} and the
#' S-map coefficient labelling so the two never disagree.
#'
#' @keywords internal
#' @noRd
#------------------------------------------------------------------------
EmbedColumnNames <- function(columns, E, tau, embedded = FALSE) {
  if (embedded) return(columns)
  sgn  <- if (tau < 0) "t-" else "t+"
  lags <- (0:(E - 1)) * abs(tau)
  unlist(lapply(columns, function(c) sprintf("%s(%s%d)", c, sgn, lags)))
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
