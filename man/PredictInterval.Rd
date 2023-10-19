\name{PredictInterval}
\alias{PredictInterval}
\title{Forecast interval accuracy}
\usage{
PredictInterval(pathIn = "./", dataFile = "", dataFrame = NULL, pathOut = "./", 
  predictFile = "", lib = "", pred = "", maxTp = 10, E = 1, tau = -1,
  exclusionRadius = 0, columns = "", target = "", embedded = FALSE,
  verbose = FALSE, validLib = vector(), numThreads = 4, showPlot = TRUE,
  noTime = FALSE)  
}
\arguments{
\item{pathIn}{path to \code{dataFile}.}

\item{dataFile}{.csv format data file name. The first column must be a time
index or time values. The first row must be column names.}

\item{dataFrame}{input data.frame. The first column must be a time
index or time values. The columns must be named.}

\item{pathOut}{path for \code{predictFile} containing output predictions.}

\item{predictFile}{output file name.}

\item{lib}{string with start and stop indices of input data rows used to
create the library of observations. A single contiguous range is supported.}

\item{pred}{string with start and stop indices of input data rows used for
predictions. A single contiguous range is supported.}

\item{maxTp}{maximum value of Tp to evalulate.}

\item{E}{embedding dimension.}

\item{tau}{lag of time delay embedding specified as number of
time column rows.}

\item{exclusionRadius}{excludes vectors from the search space of nearest 
neighbors if their relative time index is within exclusionRadius.}

\item{columns}{string of whitespace separated column name(s) in the
input data used to create the library.}

\item{target}{column name in the input data used for prediction.}

\item{embedded}{logical specifying if the input data are embedded.}

\item{verbose}{logical to produce additional console reporting.}

\item{validLib}{logical vector the same length as the number of data
 rows.  Any data row represented in this vector as FALSE, will not be
 included in the library.}

\item{numThreads}{number of parallel threads for computation.}

\item{showPlot}{logical to plot results.}

\item{noTime}{logical to allow input data with no time column.}
}

\value{
  A data.frame with columns \code{Tp, rho}.
}

\description{
\code{\link{PredictInterval}} uses \code{\link{Simplex}} to evaluate
prediction accuracy as a function of forecast interval Tp.
}

%\details{
%}

\examples{
data(TentMap)
Tp.rho = PredictInterval( dataFrame = TentMap, lib = "1 100",
pred = "201 500", E = 2, columns = "TentMap", target = "TentMap",
showPlot = FALSE )
}
