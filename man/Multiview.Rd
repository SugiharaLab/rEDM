\name{Multiview}
\alias{Multiview}
\title{Forecasting using multiview embedding}
\usage{
Multiview(pathIn = "./", dataFile = "", dataFrame = NULL,
  lib = "", pred = "", D = 0, E = 1, Tp = 1, knn = 0, 
  tau = -1, columns = "", target = "", multiview = 0, exclusionRadius = 0,
  trainLib = TRUE, excludeTarget = FALSE, parameterList = FALSE,
  verbose = FALSE, numThreads = 4, showPlot = FALSE)
}
\arguments{
\item{pathIn}{path to \code{dataFile}.}

\item{dataFile}{.csv format data file name. The first column must be a time
index or time values. The first row must be column names.}

\item{dataFrame}{input data.frame. The first column must be a time
index or time values. The columns must be named.}

\item{lib}{a 2-column matrix, data.frame, 2-element vector or string of 
  row indice pairs, where each pair specifies the first and last *rows* of
  the time series to create the library.}

\item{pred}{(same format as lib), but specifying the sections of the time 
series to forecast.}

\item{D}{multivariate dimension.}

\item{E}{embedding dimension.}

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

\item{trainLib}{logical to use in-sample (lib=pred) projections for the
  ranking of column combinations.}

\item{excludeTarget}{logical to exclude embedded target column from combinations.}

\item{parameterList}{logical to add list of invoked parameters.}

\item{verbose}{logical to produce additional console reporting.}

\item{numThreads}{number of CPU threads to use in multiview processing.}

\item{showPlot}{logical to plot results.}
}

\value{
Named list with data.frames \code{[[View, Predictions]]}.

data.frame \code{View} columns:
\tabular{ll}{
  Col_1 \tab column index\cr
  ...   \tab column index\cr
  Col_D \tab column index\cr
  rho   \tab Pearson correlation\cr
  MAE   \tab mean absolute error\cr
  RMSE  \tab root mean square error\cr
  name_1 \tab column name\cr
  ...    \tab column name\cr
  name_D \tab column name\cr
}

  If \code{parameterList = TRUE} a named list "parameters" is added.
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
    \item Compute D-dimensional variable combination forecasts.
    \item Rank forecasts.
    \item Compute predictions of top combinations.
    \item Compute multiview averaged prediction.
  }
  If \code{E>1}, all variables are embedded to dimension E. 
  If \code{trainLib} is \code{TRUE} initial forecasts and ranking are
  done in-sample (\code{lib=pred}) and predictions using the top ranked
  combinations use the specified \code{lib} and \code{pred}.
  If \code{trainLib} is \code{FALSE} initial forecasts and ranking use
  the specified \code{lib} and \code{pred}, the step of computing
  predictions of the top combinations is skipped. 
}

\examples{
data(block_3sp)
L = Multiview( dataFrame = block_3sp, lib = "1 100", pred = "101 190",
E = 2, columns = "x_t y_t z_t", target = "x_t" )
}
