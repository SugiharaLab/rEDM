\name{Multiview}
\alias{Multiview}
\title{Multiview embedding forecast}
\description{Multiview embedding forecast (Ye & Sugihara 2016): rank
D-column views and average the best to find optimal combinations of
variables that best represent the dynamics.}
\usage{
Multiview(dataFrame = NULL, columns, target, lib = numeric(0), pred = numeric(0),
          D = 0, E = 1, Tp = 1, knn = 0, tau = -1, multiview = 0,
          exclusionRadius = 0, trainLib = TRUE, excludeTarget = FALSE,
          ignoreNan = TRUE, numProcess = 4, backend = "RANN",
          pathIn = "./", dataFile = "", pathOut = "./", predictFile = "",
          parameterList = FALSE, showPlot = FALSE)
}
\arguments{
\item{dataFrame}{A data.frame of input data. The first column is time unless \code{noTime = TRUE}.}
\item{columns}{Column name(s) to build the embedding: character vector
or space-separated string. If individual column names
contain whitespace place names in a vector, or, append ',' to the name.}
\item{target}{Target column name to predict.}
\item{lib}{Library (training) index range as (start end) pairs.}
\item{pred}{Prediction index range as (start end) pairs.}
\item{D}{View dimension for Multiview (number of embedding columns combined).}
\item{E}{Embedding dimension.}
\item{Tp}{Forecast interval (prediction horizon).}
\item{knn}{Number of nearest neighbours. \code{0} selects the method default.}
\item{tau}{Embedding delay (negative selects past lags).}
\item{multiview}{Number of top-ranked views to average (default \code{floor(sqrt(nCombo))}).}
\item{exclusionRadius}{Temporal (Theiler) exclusion radius around each prediction point.}
\item{trainLib}{Rank views on the library (in-sample) rather than the prediction set.}
\item{excludeTarget}{Exclude the target variable from candidate views.}
\item{ignoreNan}{Remove rows with NaN in the embedding from the library and prediction sets.}
\item{numProcess}{Number of worker processes for the parameter sweep or task grid.}
\item{backend}{Nearest-neighbour backend: \code{"RANN"} (default) or \code{"brute"}.}
\item{pathIn}{File path for input \code{dataFile}.}
\item{dataFile}{Input dataFile, .csv format.}
\item{pathOut}{Output file path for \code{predictFile}}
\item{predictFile}{Output file name, .csv format.}
\item{parameterList}{Append named list of parameters/values to return.}
\item{showPlot}{If \code{TRUE}, draw a base-graphics plot of the result.}
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

\keyword{nonlinear,ts}
