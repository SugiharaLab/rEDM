\name{PredictNonlinear}
\alias{PredictNonlinear}
\title{Test for nonlinear dynamics}
\usage{
PredictNonlinear(pathIn = "./", dataFile = "", dataFrame = NULL,
  pathOut = "./",  predictFile = "", lib = "", pred = "", theta = "",
  E = 1, Tp = 1, tau = 1, columns = "", target = "", embedded = FALSE,
  verbose = FALSE, numThreads = 4, showPlot = TRUE)  
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

\item{theta}{A whitespace delimeted string with values of the S-map 
  localisation parameter. An empty string will use default values of
\code{[0.01 0.1 0.3 0.5 0.75 1 1.5 2 3 4 5 6 7 8 9]}.}

\item{E}{embedding dimension.}

\item{Tp}{prediction horizon (number of time column rows).}

\item{tau}{lag of time delay embedding specified as number of
time column rows.}

\item{columns}{string of whitespace separated column name(s) in the
input data used to create the library.}

\item{target}{column name in the input data used for prediction.}

\item{embedded}{logical specifying if the input data are embedded.}

\item{verbose}{logical to produce additional console reporting.}

\item{numThreads}{number of parallel threads for computation.}

\item{showPlot}{logical to plot results.}
}

\value{
  A data.frame with columns \code{Theta, rho}.
}

\description{
\code{\link{PredictNonlinear}} uses \code{\link{SMap}} to evaluate
prediction accuracy as a function of the localisation parameter
\code{theta}.
}

\details{The localisation parameter \code{theta} weights nearest
  neighbors according to exp( (-theta D / D_avg) ) where D is the
  distance between the observation vector and neighbor, D_avg the mean
  distance.  If theta = 0, weights are uniformally unity corresponding
  to a global autoregressive model.  As theta increases, neighbors in
  closer proximity to the observation are considered. 
}

\examples{
data(TentMapNoise)
theta.rho <- PredictNonlinear( dataFrame=TentMapNoise, E=2,lib="1 100",
pred="201 500", columns="TentMap", target="TentMap", showPlot = FALSE)
}
