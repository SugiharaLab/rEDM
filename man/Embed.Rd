\name{Embed}
\alias{Embed}
\title{Embed data with time lags}

\usage{
Embed(path = "./", dataFile = "", dataFrame = NULL, E = 0, tau = -1, 
columns = "", verbose = FALSE)  
}

\arguments{
\item{path}{path to \code{dataFile}.}

\item{dataFile}{.csv format data file name. The first column must be a time
index or time values. The first row must be column names. One of
\code{dataFile} or \code{dataFrame} are required.}

\item{dataFrame}{input data.frame. The first column must be a time
index or time values. The columns must be named. One of
\code{dataFile} or \code{dataFrame} are required.}

\item{E}{embedding dimension.}

\item{tau}{integer time delay embedding lag specified as number of
time column rows.}

\item{columns}{string of whitespace separated column name(s), or vector
of column names used to create the library. If individual column names
contain whitespace place names in a vector, or, append ',' to the name.}

\item{verbose}{logical to produce additional console reporting.}
}

\description{
\code{\link{Embed}} performs Takens time-delay embedding on \code{columns}.
}

\value{
A data.frame with lagged columns. E columns for each variable specified
in \code{columns}.
}

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
data(circle)
embed <- Embed( dataFrame = circle, E = 2, tau = -1, columns = "x y" ) 
}
