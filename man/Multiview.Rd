\name{Multiview}
\alias{Multiview}
\title{Forecasting using multiview embedding}
\usage{
Multiview(pathIn = "./", dataFile = "", dataFrame = NULL, pathOut = "./", 
  predictFile = "", lib = "", pred = "", E = 0, Tp = 1, knn = 0, 
  tau = 1, columns = "", target = "", multiview = 0, exclusionRadius = 0,
  verbose = FALSE, numThreads = 4, showPlot = FALSE)
}
\arguments{
\item{pathIn}{path to \code{dataFile}.}

\item{dataFile}{.csv format data file name. The first column must be a time
index or time values. The first row must be column names.}

\item{dataFrame}{input data.frame. The first column must be a time
index or time values. The columns must be named.}

\item{pathOut}{path for \code{predictFile} containing output predictions.}

\item{predictFile}{prediction output file name.}

\item{lib}{string with start and stop indices of input data rows used to
create the library of observations. A single contiguous range is supported.}

\item{pred}{string with start and stop indices of input data rows used for
predictions. A single contiguous range is supported.}

\item{E}{multivariate dimension.}

\item{Tp}{prediction horizon (number of time column rows).}

\item{knn}{number of nearest neighbors. If knn=0, knn is set to E+1.}

\item{tau}{lag of time delay embedding specified as number of
time column rows.}

\item{columns}{string of whitespace separated column name(s) in the
input data used to create multivariable data sets.}

\item{target}{column name in the input data used for prediction.}

\item{multiview}{number of multiview ensembles to average for the final
prediction estimate.}

\item{exclusionRadius}{number of adjacent observation vector rows to exclude as
nearest neighbors in prediction.}

\item{verbose}{logical to produce additional console reporting.}

\item{numThreads}{number of CPU threads to use in multiview processing.}

\item{showPlot}{logical to plot results.}
}

\value{
Named list with data.frames \code{[[Combo_rho, Predictions]]}.

data.frame \code{Combo_rho} columns:
\tabular{ll}{
  Col_1 \tab column index\cr
  ...   \tab column index\cr
  Col_E \tab column index\cr
  rho   \tab Pearson correlation\cr
  MAE   \tab mean absolute error\cr
  RMSE  \tab root mean square error\cr
}
}

\references{Ye H., and G. Sugihara, 2016. Information leverage in 
interconnected ecosystems: Overcoming the curse of dimensionality.
Science 353:922-925.
}

\description{
\code{\link{Multiview}} applies the method of \cite{Ye & Sugihara}
to find optimal combinations of variables that best represent the
dynamics.
}

\details{Multiview embedding is a method to identify variables in a
  multivariate dynamical system that are most likely to contribute to
  the observed dynamics.  It is a multistep algorithm with these general
  steps:
  \enumerate{
    \item Evaluate E-dimensional variable combinations within the
    library itself (in-sample).
    \item Rank in-sample (library) forecasts.
    \item Compute Multiview averaged out-of-sample prediction
  }
}

\examples{
data(block_3sp)
L = Multiview( dataFrame = block_3sp, lib = "1 99", pred = "105 190",
E = 2, columns = "x_t y_t z_t", target = "x_t" )
}
