\name{MakeBlock}
\alias{MakeBlock}
\title{Make embedded data block}

\usage{
MakeBlock(dataFrame, E = 0, tau = 1, columns = "", verbose = FALSE)  
}

\arguments{
\item{dataFrame}{input data.frame. The first column must be a time
index or time values. The columns must be named.}

\item{E}{embedding dimension.}

\item{tau}{positive integer time delay embedding lag specified as number of
time column rows.}

\item{columns}{string of whitespace separated column name(s) in the
input data to be embedded.}

\item{verbose}{logical to produce additional console reporting.}
}

\description{
\code{\link{MakeBlock}} performs Takens time-delay embedding on
\code{columns}. It is an internal function called by \code{\link{Embed}}
that does not perform input error checking or validation. 
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
data(TentMap)
embed <- MakeBlock(TentMap, 3, 1, "TentMap")
}
