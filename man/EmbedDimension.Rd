\name{EmbedDimension}
\alias{EmbedDimension}
\title{Optimal embedding dimension}
\usage{
EmbedDimension(pathIn = "./", dataFile = "", dataFrame = NULL, pathOut = "", 
  predictFile = "", lib = "", pred = "", maxE = 10, Tp = 1, tau = -1,
  exclusionRadius = 0, columns = "", target = "", embedded = FALSE,
  verbose = FALSE, validLib = vector(), numThreads = 4, showPlot = TRUE,
  noTime = FALSE)  
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

\item{maxE}{maximum value of E to evalulate.}

\item{Tp}{prediction horizon (number of time column rows).}

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
  A data.frame with columns \code{E, rho}.
}

\description{
\code{\link{EmbedDimension}} uses \code{\link{Simplex}} to evaluate
prediction accuracy as a function of embedding dimension.
}

%\details{
%}

\examples{
data(TentMap)
E.rho = EmbedDimension( dataFrame = TentMap, lib = "1 100", pred = "201 500",
columns = "TentMap", target = "TentMap", showPlot = FALSE )
}
