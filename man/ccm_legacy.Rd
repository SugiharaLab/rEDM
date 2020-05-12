\name{ccm}
\alias{ccm}
\title{Convergent cross mapping using simplex projection}
\usage{
ccm(block, lib = NULL, pred = NULL, norm = 2, E = 1, tau = -1, 
    tp = 0, num_neighbors = "e+1", lib_sizes = c(10, 75, 5), 
    random_libs = TRUE, num_samples = 100, replace = FALSE, lib_column = 1, 
    target_column = 2, first_column_time = FALSE, RNGseed = NULL, 
    exclusion_radius = NULL, epsilon = NULL, stats_only = TRUE, 
    silent = TRUE)
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

\item{E}{the embedding dimensions to use for time delay embedding}

\item{tau}{the time-delay offset to use for time delay embedding}

\item{tp}{the prediction horizon (how far ahead to forecast)}

\item{num_neighbors}{the number of nearest neighbors to use. Note that the 
default value will change depending on the method selected. (any of "e+1", 
"E+1", "e + 1", "E + 1" will set this parameter to E+1 for each run}

\item{lib_sizes}{three integers specifying the start, stop and increment
index of library sizes}

\item{random_libs}{indicates whether to use randomly sampled libs}

\item{num_samples}{is the number of random samples at each lib size (this 
parameter is ignored if random_libs is FALSE)}

\item{replace}{indicates whether to sample vectors with replacement}

\item{lib_column}{name (index) of the column to cross map from}

\item{target_column}{name (index) of the column to forecast}

\item{first_column_time}{indicates whether the first column of the given 
block is a time column}

\item{RNGseed}{will set a seed for the random number generator, enabling 
reproducible runs of ccm with randomly generated libraries}

\item{exclusion_radius}{excludes vectors from the search space of nearest 
neighbors if their *time index* is within exclusion_radius (NULL turns 
this option off)}

\item{epsilon}{not implemented}

\item{stats_only}{specify whether to output just the forecast statistics or 
the raw predictions for each run}

\item{silent}{prevents warning messages from being printed to the R console}
}
\value{
If stats_only = TRUE: a data.frame with forecast statistics for both the
forward and reverse mappings:
\tabular{ll}{
  LibSize \tab library length (number of vectors)\cr
  x:y \tab cross mapped correlation coefficient between observations x and
  predictions y\cr
  y:x \tab cross mapped correlation coefficient between observations y and
  predictions x\cr
  E \tab embedding dimension\cr
  tau \tab time delay offset\cr
  tp \tab forecast interval\cr
  nn \tab number nearest neighbors\cr
}
If stats_only = FALSE: a named list with the following items:
  settings:
\tabular{ll}{
  LibMeans \tab data.frame with the mean bidirectional forecast statistics\cr
  CCM1_PredictStat \tab data.frame with forward mapped prediction
  statistics for each prediction of the ensemble\cr
  CCM1_Predictions \tab list of prediction result data.frame each
  forward mapped prediction of the ensemble\cr
  CCM2_PredictStat \tab data.frame with reverse mapped prediction
  statistics for each prediction of the ensemble\cr
  CCM2_Predictions \tab  list of prediction result data.frame each
  reverse mapped prediction of the ensemble\cr
}
CCM1_PredictStat and CCM2_PredictStat data.frames have
columns:
\tabular{ll}{
  N \tab prediction number\cr
  E \tab embedding dimension\cr
  nn \tab number of nearest neighbors\cr
  tau \tab embedding time delay offset\cr
  LibSize \tab library size\cr
  rho \tab correlation coefficient\cr
  RMSE \tab root mean square error\cr
  MAE \tab maximum absolute error\cr
  lib \tab column name of the library vector\cr
  target \tab column name of the target vector\cr
}
}
\description{
\code{\link{ccm}} uses time delay embedding on one time series to generate an 
attractor reconstruction, and then applies the simplex projection algorithm 
to estimate concurrent values of another time series. This method is 
typically applied, varying the library sizes, to determine if one time series
contains the necessary dynamic information to recover the influence of 
another, causal variable.
}
\details{
\code{\link{ccm}} runs both forward and reverse cross maps in seperate
threads. Results are returned for both mappings. 
The default parameters are set so that passing a matrix as the only argument
will use E = 1 (embedding dimension), and leave-one-out cross-validation over
the whole time series to compute cross-mapping from the first column to the 
second column, letting the library size vary from 10 to 75 in increments of 
5.

\code{norm = 2} (only option currently available) uses the "L2 norm",
Euclidean distance:
  \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
  }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
}
\examples{
anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, 
  lib_column = "anchovy", target_column = "np_sst", 
  lib_sizes = c(10, 75, 5), num_samples = 100)
}
