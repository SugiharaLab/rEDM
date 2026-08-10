\name{Embed}
\alias{Embed}
\title{Time-delay embedding}
\description{Takens time-delay embedding of one or more columns, returned as a
named data.frame matching pyEDM's \code{Embed()}.}
\usage{
Embed(dataFrame = NULL, E = 0, tau = -1, columns = "",
      includeTime = FALSE, pathIn = "./", dataFile = "")
}
\arguments{
\item{dataFrame}{A data.frame of input data. Column 1 is the time vector.}
\item{E}{Embedding dimension.}
\item{tau}{Embedding delay in time-column rows (negative selects past lags).}
\item{columns}{Column name(s) to embed: a whitespace-separated string, or a
character vector of column names.}
\item{includeTime}{If \code{TRUE}, prepend column 1 of \code{dataFrame} (the
time vector), keeping its original name, as the first output column.}
\item{pathIn}{Directory for \code{dataFile} (used when reading input from disk).}
\item{dataFile}{CSV file name to read as input; when set it takes precedence
over \code{dataFrame}.}
}
\value{A data.frame of the embedding (with the time column first when
\code{includeTime = TRUE}).}
\details{
  Each \code{columns} item produces E lagged vectors. For column \code{c} and
  step \code{i = 0..E-1} the lag is \code{i * |tau|}, and the column is named
  \code{c(t-lag)} when \code{tau < 0} or \code{c(t+lag)} when \code{tau > 0}.
  For example columns X, Y with E = 2, tau = -1 yield
  \code{X(t-0) X(t-1) Y(t-0) Y(t-1)}; with tau = -2 they yield
  \code{X(t-0) X(t-2) Y(t-0) Y(t-2)}.

  If \code{dataFile} is supplied the CSV at \code{file.path(pathIn, dataFile)}
  is read (preserving column names) and used as the input, taking precedence
  over \code{dataFrame}.

  All rows are retained; the first \code{|tau| * (E-1)} partial rows contain
  \code{NA}.
}
\examples{
Embed(data.frame(x = 1:10), E = 3, tau = -1, columns = "x")
}
\keyword{nonlinear}
