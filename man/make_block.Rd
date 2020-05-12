\name{make_block}
\alias{make_block}
\title{Make a time delay offset block}
\usage{
make_block(block, columns = NULL, t = NULL, max_lag = 3, tau = -1, lib =
NULL, restrict_to_lib = TRUE)
}
\arguments{
\item{block}{a data.frame or matrix where each column is a time series}

\item{columns}{list of column names to time delay.}

\item{t}{Not used}

\item{max_lag}{the total number of lags to include for each variable. So if 
max_lag == 3, a variable X is offset with lags X[t], X[t + tau], 
X[t + 2*tau]}

\item{tau}{the time delay offset for embedding}

\item{lib}{not used}

\item{restrict_to_lib}{not used}
}
\value{
A data.frame with time offset columns. If the 
  original block had columns X, Y and max_lag = 3, then the returned 
  data.frame will have columns X(t-0) X(t-1) X(t-2) Y(t-0) Y(t-1) Y(t-2).
}
\description{
\code{\link{make_block}} generates a time offset block with the
appropriate max_lag and tau. The first column is presumed to be a time
or index vector, and is not included in the embedding.
}
\examples{
data("block_3sp")
make_block(block_3sp[, c(1, 2, 5)])
}
