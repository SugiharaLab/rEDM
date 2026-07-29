\name{Simplex}
\alias{Simplex}
\title{Simplex projection}

\description{
\code{\link{Simplex}} performs time series forecasting based on
weighted nearest neighbors projection in the time series phase space as
described in \cite{Sugihara and May}.
}

\usage{
Simplex(dataFrame = NULL, columns, target, lib, pred, E = 0, Tp = 1, knn = 0,
        tau = -1, exclusionRadius = 0, embedded = FALSE,
        validLib = logical(0), noTime = FALSE, ignoreNan = TRUE,
        backend = "RANN", pathIn = "./", dataFile = "", pathOut = "./",
        predictFile = "", parameterList = FALSE, verbose = FALSE, showPlot = FALSE)
}

\arguments{
\item{dataFrame}{A data.frame of input data. The first column is time
  unless \code{noTime = TRUE}. The columns must be named.}
\item{columns}{Column name(s) to build the embedding: character vector string of whitespace separated column name(s), or vector
of column names used to create the library. If individual column names
contain whitespace place names in a vector, or, append ',' to the name.}
\item{target}{Target column name to predict.}
\item{lib}{Library (training) index range as (start end) pairs.  Mulitple row index
pairs can be specified with each pair defining the first and last
rows of time series observation segments used to create the library.}
\item{pred}{Prediction index range as (start end) pairs. A single contiguous range is supported.}
\item{E}{Embedding dimension.}
\item{Tp}{Forecast interval (prediction horizon, number of time column rows).}
\item{knn}{Number of nearest neighbours. If \code{knn=0} knn is set to E+1..}
\item{tau}{Embedding delay (negative selects past lags).}
\item{exclusionRadius}{Temporal (Theiler) exclusion radius around each
  prediction point excludes vectors from the search space of nearest 
neighbors if their relative time index is within \code{exclusionRadius}.}
\item{embedded}{If \code{TRUE}, \code{columns} already form the embedding.}
\item{validLib}{Logical vector marking admissible library rows (or length 0 for all).}
\item{noTime}{If \code{TRUE}, synthesise a 1..N time index instead of using column 1.}
\item{ignoreNan}{Remove rows with NaN in the embedding from the library
  and prediction sets.}
\item{backend}{Nearest-neighbour backend: \code{"RANN"} (default) or \code{"brute"}.}
\item{pathIn}{File path for input \code{dataFile}.}
\item{dataFile}{Input dataFile, .csv format. The first column must be a time
index or time values unless noTime is TRUE. The first row must be column names.}
\item{pathOut}{Output file path for \code{predictFile}}
\item{predictFile}{Output file name, .csv format.}
\item{parameterList}{Append named list of parameters/values to return.}
\item{verbose}{Emit diagnostic messages.}
\item{showPlot}{If \code{TRUE}, draw a base-graphics plot of the result.}
}

\value{A data.frame with columns Time, Observations, Predictions and
  Pred_Variance.}

\examples{
data( block_3sp )
smplx = Simplex(block_3sp, "x_t", "x_t", lib = c(1, 100), pred = c(101, 195), E = 3)
ComputeError( smplx $ Predictions, smplx $ Observations )
}

\details{
  If embedded is \code{FALSE}, the data \code{column(s)} are embedded to
  dimension \code{E} with time lag \code{tau}. This embedding forms an
  E-dimensional phase space for the \code{\link{Simplex}} projection.
  If embedded is \code{TRUE}, the data are assumed to contain an
  E-dimensional embedding with E equal to the number of \code{columns}.
  Predictions are made using leave-one-out cross-validation, i.e.
  observation vectors are excluded from the prediction simplex. 

  To assess an optimal embedding dimension \code{\link{EmbedDimension}}
  can be applied. Accuracy statistics can be estimated by
  \code{\link{ComputeError}}.
}

\references{Sugihara G. and May R. 1990. Nonlinear forecasting as a way
of distinguishing chaos from measurement error in time series.
Nature, 344:734-741.
}

\keyword{nonlinear}
