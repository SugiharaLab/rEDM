\name{Embed}
\alias{Embed}
\title{Time-delay embedding}
\description{Takens time-delay embedding of one or more columns.}
\usage{
Embed(dataFrame, E, tau, columns)
}
\arguments{
\item{dataFrame}{A data.frame of input data. The first column is time unless \code{noTime = TRUE}.}
\item{E}{Embedding dimension.}
\item{tau}{Embedding delay specified as number of
time column rows. (negative selects past lags).}
\item{columns}{Column name(s) to build the embedding: string of
whitespace separated column name(s), or vector
of column names used to create the library. If individual column names
contain whitespace place names in a vector, or, append ',' to the name.}
}
\value{A numeric matrix of the embedding (NA in the first/last |tau|*(E-1) rows).}
\details{
  Each \code{columns} item will have E-1 time-lagged vectors created.
  The column name is appended with \code{(t-n)}.  For example, data
  columns X, Y, with E = 2 will have columns named
  \code{X(t-0) X(t-1) Y(t-0) Y(t-1)}.

  The returned data.frame does not have a time column.  The returned
  data.frame is truncated by tau * (E-1) rows to remove state vectors
  with partial data (NaN elements).
}
\examples{
Embed(data.frame(x = 1:10), E = 3, tau = -1, columns = "x")
}
\keyword{nonlinear,ts}
