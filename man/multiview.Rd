\name{multiview}
\alias{multiview}
\title{Perform forecasting using multiview embedding}
\usage{
multiview(block, lib = NULL, pred = NULL, norm = 2, E = 0, tau = -1, 
    tp = 1, max_lag = 3, num_neighbors = "e+1", k = "sqrt", na.rm = FALSE, 
    target_column = 1, stats_only = TRUE, save_lagged_block = FALSE, 
    first_column_time = FALSE, exclusion_radius = NULL, silent = FALSE)
}
\arguments{
\item{block}{either a vector to be used as the time series, or a 
data.frame or matrix where each column is a time series}

\item{lib}{a 2-column matrix, data.frame, 2-element vector or string of 
  row indice pairs, where each pair specifies the first and last *rows* of
  the time series to create the library. If not specified, all available
  rows are used}

\item{pred}{(same format as lib), but specifying the sections of the time 
series to forecast. If not specified, set equal to lib}

\item{norm}{the distance measure to use. see 'Details'}

\item{E}{the embedding dimensions to use for time delay embedding. The
  default value of 1 does not embed the data.}

\item{tau}{the time-delay offset to use for time delay embedding}

\item{tp}{the prediction horizon (how far ahead to forecast)}

\item{max_lag}{the maximum number of lags to use for variable combinations. 
If max_lag == 3, a variable X will be embedded with lags X[t], X[t + tau], 
X[t + 2*tau]}

\item{num_neighbors}{the number of nearest neighbors to use. Note that the 
default value will change depending on the method selected. (any of "e+1", 
"E+1", "e + 1", "E + 1" will set this parameter to E+1 for each run.)}

\item{k}{the number of embeddings to use for ensemble averaging.
  "sqrt" or 0 will use k = sqrt(m) where m is the number of
  multiview combinations of the set of input variables}

\item{na.rm}{logical. Should missing values (including `NaN`` be omitted 
from the calculations?)}

\item{target_column}{the name (index) of the column to forecast}

\item{stats_only}{specify whether to output just the forecast statistics or 
the raw predictions for each run}

\item{save_lagged_block}{specify whether to output the lagged block that 
is constructed as part of running \code{multiview}}

\item{first_column_time}{indicates whether the first column of the given 
block is a time column and excluded when building the library}

\item{exclusion_radius}{excludes vectors from the search space of nearest 
neighbors if their *time index* is within exclusion_radius (NULL turns 
this option off)}

\item{silent}{prevents warning messages from being printed to the R console}
}
\value{
A named list with items "View" and "Predictions".  View is a data.frame with components:
\tabular{ll}{
  col_i,... col_j \tab column indices of the embedding\cr
  name_i,... nam_j \tab column names of the embedding\cr
  rho \tab correlation of the projection\cr
  MAE \tab maximum absolute error of the projection\cr
  RMSE \tab root mean square error of the projection\cr
}
Predictions is a data.frame of the predictions from the best multivew ensemble.
}
\description{
\code{\link{multiview}} applies the method described in Ye & Sugihara (2016) for 
  forecasting, where multiple attractor reconstructions are tested, and a 
  single nearest neighbor is selected from each of the top \code{k} 
  reconstructions to produce final forecasts.
}
\details{
  \code{\link{multiview}} uses multiple time series given as input to
  generate an attractor 
  reconstruction, and then applies the simplex projection to make
  forecasts. This method generalizes the \code{\link{simplex}} routine,
  and allows for "mixed" embeddings, where multiple 
  time series can be used as different dimensions of an attractor 
  reconstruction.

The default parameters are set so that, given a matrix of time series, 
  forecasts will be produced for the first column. By default, all possible 
  combinations of the columns are used for the attractor construction, the 
  k = sqrt(m) heuristic will be used, forecasts will be one time step ahead. 
  If a time vector is not supplied, 1:NROW will be used.
  The default lib and pred are 
  to use the first half of the data for the "library" and to predict over the
  second half of the data. Unless otherwise set, the output will be just the 
  forecast statistics.

\code{norm = 2} (only option currently available) uses the "L2 norm",
Euclidean distance:
  \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
  }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
}
\examples{
block <- block_3sp[, c(2, 5, 8)]
multiview( block, k=10 )
}
