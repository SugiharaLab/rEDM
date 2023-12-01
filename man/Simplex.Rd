\name{Simplex}
\alias{Simplex}
\title{Simplex forecasting}
\usage{
Simplex(pathIn = "./", dataFile = "", dataFrame = NULL, pathOut = "./", 
  predictFile = "", lib = "", pred = "", E = 0, Tp = 1, knn = 0, tau = -1, 
  exclusionRadius = 0, columns = "", target = "", embedded = FALSE,
  verbose = FALSE, validLib = vector(), generateSteps = 0,
  parameterList = FALSE, showPlot = FALSE, noTime = FALSE)
}
\arguments{
\item{pathIn}{path to \code{dataFile}.}

\item{dataFile}{.csv format data file name. The first column must be a time
index or time values unless noTime is TRUE. The first row must be column names.}

\item{dataFrame}{input data.frame. The first column must be a time
index or time values unless noTime is TRUE. The columns must be named.}

\item{pathOut}{path for \code{predictFile} containing output predictions.}

\item{predictFile}{output file name.}

\item{lib}{string or vector with start and stop indices of input data
rows used to create the library from observations. Mulitple row index
pairs can be specified with each pair defining the first and last
rows of time series observation segments used to create the library.}

\item{pred}{string with start and stop indices of input data rows used for
predictions. A single contiguous range is supported.}

\item{E}{embedding dimension.}

\item{Tp}{prediction horizon (number of time column rows).}

\item{knn}{number of nearest neighbors. If knn=0, knn is set to E+1.}

\item{tau}{lag of time delay embedding specified as number of
time column rows.}

\item{exclusionRadius}{excludes vectors from the search space of nearest 
neighbors if their relative time index is within exclusionRadius.}

\item{columns}{string of whitespace separated column name(s), or vector
of column names used to create the library. If individual column names
contain whitespace place names in a vector, or, append ',' to the name.}

\item{target}{column name used for prediction.}

\item{embedded}{logical specifying if the input data are embedded.}

\item{verbose}{logical to produce additional console reporting.}

\item{validLib}{logical vector the same length as the number of data
 rows.  Any data row represented in this vector as FALSE, will not be
 included in the library.}

\item{generateSteps}{number of predictive feedback generative steps.}

\item{parameterList}{logical to add list of invoked parameters.}

\item{showPlot}{logical to plot results.}

\item{noTime}{logical to allow input data with no time column.}
}

\value{
A data.frame with columns \code{Observations, Predictions}.
The first column contains the time values.

If \code{parameterList = TRUE}, a named list with "predictions" holding the
data.frame, "parameters" with a named list of invoked parameters.
}

\references{Sugihara G. and May R. 1990. Nonlinear forecasting as a way
of distinguishing chaos from measurement error in time series.
Nature, 344:734-741.
}

\description{
\code{\link{Simplex}} performs time series forecasting based on
weighted nearest neighbors projection in the time series phase space as
described in \cite{Sugihara and May}.
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
\examples{
data( block_3sp )
smplx = Simplex( dataFrame = block_3sp, lib = "1 100", pred = "101 190",
E = 3, columns = "x_t", target = "x_t" )
ComputeError( smplx $ Predictions, smplx $ Observations )
}
