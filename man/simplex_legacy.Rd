\name{simplex}
\alias{simplex}
\alias{s_map}
\title{Perform univariate forecasting}
\usage{
simplex(time_series, lib = NULL, pred = NULL, norm = 2, E = 1:10, 
    tau = -1, tp = 1, num_neighbors = "e+1", stats_only = TRUE, 
    exclusion_radius = NULL, epsilon = NULL, silent = TRUE)

s_map(time_series, lib = NULL, pred = NULL, norm = 2, E = 1, 
    tau = -1, tp = 1, num_neighbors = 0, theta = NULL, stats_only = TRUE, 
    exclusion_radius = NULL, epsilon = NULL, silent = TRUE,
    save_smap_coefficients = FALSE)
}
\arguments{
\item{time_series}{either a vector to be used as the time series, or a 
data.frame or matrix with at least 2 columns (in which case the first 
column will be used as the time index, and the second column as the time 
series)}

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
"E+1", "e + 1", "E + 1" will set this parameter to E+1.)}

\item{stats_only}{specify whether to output just the forecast statistics or 
the raw predictions for each run}

\item{exclusion_radius}{excludes vectors from the search space of nearest 
neighbors if their *time index* is within exclusion_radius (NULL turns 
this option off)}

\item{epsilon}{Deprecated.}

\item{silent}{prevents warning messages from being printed to the R console}

\item{theta}{the nonlinear tuning parameter (theta is only relevant if 
method == "s-map")}

\item{save_smap_coefficients}{specifies whether to include the s_map 
coefficients with the output}
}
\value{
For \code{\link{simplex}}, if stats_only = TRUE: a data.frame with components for the parameters and forecast statistics:
\tabular{ll}{
  E \tab embedding dimension\cr
  tau \tab embedding time offset\cr
  tp \tab prediction horizon\cr
  nn \tab number of neighbors\cr
  num_pred \tab number of predictions\cr
  rho \tab correlation coefficient between observations and 
    predictions\cr
  mae \tab mean absolute error\cr
  rmse \tab root mean square error\cr
  perc \tab percent correct sign\cr
  p_val \tab p-value that rho is significantly greater than 0 using 
    Fisher's z-transformation\cr
  const_pred_rho \tab same as \code{rho}, but for the constant predictor\cr
  const_pred_mae \tab same as \code{mae}, but for the constant predictor\cr
  const_pred_rmse \tab same as \code{rmse}, but for the constant predictor\cr
  const_pred_perc \tab same as \code{perc}, but for the constant predictor\cr
  const_p_val \tab same as \code{p_val}, but for the constant predictor\cr
}
For \code{\link{simplex}}, if stats_only = FALSE: a named list with
data.frame "stats" specified above, and named list "model_output":
\tabular{ll}{
  model_output \tab named list with data.frames for each model.
  Columns include the time index, observations, predictions, and
  estimated prediction variance\cr
}
For \code{\link{s_map}}, if stats_only = TRUE,
the same data.frame as for \code{\link{simplex}}, but with additional column:
\tabular{ll}{
  theta \tab the nonlinear tuning parameter\cr
}
For \code{\link{s_map}}, if save_smap_coefficients = TRUE, a named list with
data.frame "stats" specified above and the following list items:
\tabular{ll}{
  smap_coefficients \tab data.frame with columns for the s-map 
  coefficients\cr
  smap_coefficient_covariances \tab list of covariance matrices for 
  the s-map coefficients\cr
}
For \code{\link{s_map}}, if stats_only = FALSE, a named list with
data.frame "stats" specified above, and named list "model_output":
\tabular{ll}{
  model_output \tab named list with data.frames for each model.
  Columns include the time index, observations, predictions, and
  estimated prediction variance\cr
}
}
\description{
\code{\link{simplex}} uses time delay embedding on a single time 
  series to generate an attractor reconstruction, and then applies the 
  simplex projection algorithm to make forecasts.

\code{\link{s_map}} is similar to \code{\link{simplex}}, but uses the S-map 
  algorithm to make forecasts.
}
\details{
\code{\link{simplex}} is typically applied, and the embedding dimension 
varied, to find an optimal embedding dimension for the data. Thus, the 
default parameters are set so that passing a time series as the only 
argument will run over E = 1:10 (embedding dimension), using leave-one-out 
cross-validation over the whole time series, and returning just the forecast 
statistics.

\code{\link{s_map}} is typically applied, with fixed embedding
dimension, and theta 
varied, to test for nonlinear dynamics in the data. Thus, the default 
parameters are set so that passing a time series as the only argument will 
run over a default list of thetas (0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 
0.03, 0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6, and 8), using E = 1, 
leave-one-out cross-validation over the whole time series, and returning 
just the forecast statistics.

\code{norm = 2} (only option currently available) uses the "L2 norm",
Euclidean distance:
  \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
  }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}
}
\examples{
ts <- block_3sp$x_t
simplex(ts, lib = c(1, 100), pred = c(101, 190))

ts <- block_3sp$x_t
simplex(ts, stats_only = FALSE)
 
ts <- block_3sp$x_t
s_map(ts, E = 2)

ts <- block_3sp$x_t
s_map(ts, E = 2, theta = 1, save_smap_coefficients = TRUE)
}
