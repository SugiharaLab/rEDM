\name{SMap}
\alias{SMap}
\title{Sequential locally weighted global linear map (S-map)}
\description{
  \code{\link{SMap}} performs time series forecasting based on localised
  (or global) nearest neighbor projection in the time series phase space as
  described in \cite{Sugihara 1994}. 
}
\usage{
SMap(dataFrame = NULL, columns, target, lib, pred, E = NULL, Tp = 1, knn = 0,
     tau = -1, theta = 0, exclusionRadius = 0, embedded = FALSE,
     validLib = logical(0), noTime = FALSE, ignoreNan = TRUE,
     backend = "RANN", pathIn = "./", dataFile = "", pathOut = "./",
     predictFile = "", parameterList = FALSE, includeState = FALSE,
     verbose = FALSE, showPlot = FALSE)
}
\arguments{
\item{dataFrame}{A data.frame of input data. The first column must be a time
index or time values unless \code{noTime = TRUE}. The columns must be named.}
\item{columns}{Column name(s) to build the embedding: character vector
or space-separated string, or vector
of column names used to create the library. If individual column names
contain whitespace place names in a vector, or, append ',' to the name.}
\item{target}{Target column name to predict.}
\item{lib}{string or vector with start and stop indices of input data
rows used to create the library from observations. Mulitple row index
pairs can be specified with each pair defining the first and last
rows of time series observation segments used to create the library.}
\item{pred}{string with start and stop indices of input data rows used for
predictions. A single contiguous range is supported.}
\item{E}{Embedding dimension. Required (no default); must be a positive integer unless \code{embedded = TRUE}, in which case it is inferred as the number of \code{columns}.}
\item{Tp}{Forecast interval (prediction horizon).}
\item{knn}{Number of nearest neighbours. \code{0} selects the method default.}
\item{tau}{Embedding delay (negative selects past lags).}
\item{theta}{S-map localisation parameter (\code{0} = global linear map).}
\item{exclusionRadius}{Temporal (Theiler) exclusion radius around each
prediction point. Excludes vectors from the search space of nearest 
neighbors if their relative time index is within exclusionRadius.}
\item{embedded}{If \code{TRUE}, \code{columns} already form the embedding.}
\item{validLib}{Logical vector marking admissible library rows (or length 0 for all).}
\item{noTime}{If \code{TRUE}, synthesise a 1..N time index instead of using column 1.}
\item{ignoreNan}{Remove rows with NaN in the embedding from the library and prediction sets.}
\item{pathIn}{File path for input \code{dataFile}.}
\item{dataFile}{Input dataFile, .csv format. The first column must be a time
index or time values unless noTime is TRUE. The first row must be column names.}
\item{pathOut}{Output file path for \code{predictFile}}
\item{predictFile}{Output file name, .csv format.}
\item{backend}{Nearest-neighbour backend: \code{"RANN"} (default) or \code{"brute"}.}
\item{parameterList}{Append named list of parameters/values to return.}
\item{includeState}{If \code{TRUE}, also return an \code{internal} list of engine state: \code{knn_neighbors}, \code{knn_distances}, \code{lib_i}, \code{pred_i}, \code{targetVec}, \code{embedding} (1-based row indices; 0 = neighbour sentinel).}
\item{verbose}{Emit diagnostic messages.}
\item{showPlot}{If \code{TRUE}, draw a base-graphics plot of the result.}
}
\value{
  A named list with three data.frames
  \code{[[predictions, coefficients, singularValues]]}.
  \code{predictions} has columns \code{Observations, Predictions}.
  The first column contains time or index values.

  \code{coefficients} data.frame has time or index values in the first column.
  Columns 2 through E+2 (E+1 columns) are the SMap coefficients.

  \code{singularValues} data.frame has time or index values in the first column.
  Columns 2 through E+2 (E+1 columns) are the SVD singularValues. The
  first value corresponds to the SVD bias (intercept) term.

  If \code{parameterList = TRUE} a named list "parameters" is added.

  If \code{includeState = TRUE} a named list "internal" is added.
}

\references{Sugihara G. 1994. Nonlinear forecasting for the classification of natural time series. Philosophical Transactions: Physical Sciences and Engineering, 348 (1688):477-495.}

\details{
  If \code{embedded} is \code{FALSE}, the data \code{column(s)} are embedded
  to dimension \code{E} with time lag \code{tau}. This embedding forms an
  n-columns * E-dimensional phase space for the \code{\link{SMap}} projection.
  If embedded is \code{TRUE}, the data are assumed to contain an
  E-dimensional embedding with E equal to the number of \code{columns}.
  See the Note below for proper use of multivariate data (number of
  \code{columns} > 1).

  If \code{ignoreNan} is \code{TRUE}, the library (\code{lib}) is
  internally redefined to exclude nan embedding vectors. If
  \code{ignoreNan} is \code{FALSE} no library adjustment is made. The
  (\code{lib}) can be explicitly specified to exclude nan library vectors.
  
  Predictions are made using leave-one-out cross-validation, i.e.
  observation rows are excluded from the prediction regression.

  In contrast to \code{\link{Simplex}}, \code{\link{SMap}} uses all
  available neighbors and weights them with an exponential decay
  in phase space distance with exponent \code{theta}. \code{theta}=0
  uses all neighbors corresponding to a global autoregressive model.
  As \code{theta} increases, neighbors closer in vicinity to the
  observation are considered. 
}

\note{
\code{\link{SMap}} should be called with columns explicitly corresponding to
dimensions E. In the univariate case (number of \code{columns} = 1) with
default \code{embedded = FALSE}, the time series will be time-delay
embedded to dimension E, SMap coefficients correspond to each dimension. 

If a multivariate data set is used (number of \code{columns} > 1) it
must use \code{embedded = TRUE} with E equal to the number of columns.
This prevents the function from internally time-delay embedding the
multiple columns to dimension E.  If the internal time-delay embedding
is performed, then state-space columns will not correspond to the
intended dimensions in the matrix inversion, coefficient assignment,
and prediction. In the multivariate case, the user should first prepare
the embedding (using \code{\link{Embed}} for time-delay embedding), then
pass this embedding to \code{\link{SMap}} with appropriately specified
\code{columns}, \code{E}, and \code{embedded = TRUE}.
}

\keyword{nonlinear}
