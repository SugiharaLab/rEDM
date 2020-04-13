\name{block_lnlp}
\alias{block_lnlp}
\title{Perform generalized forecasting using simplex projection or s-map}
\usage{
block_lnlp(block, lib = NULL, pred = NULL, norm = 2, method = c("simplex",
  "s-map"), tp = 1, num_neighbors = switch(match.arg(method), 
  simplex = "e+1", `s-map` = 0), columns = NULL, target_column = 1, 
  stats_only = TRUE, first_column_time = FALSE, exclusion_radius = NULL, 
  epsilon = NULL, theta = NULL, silent = TRUE, save_smap_coefficients = FALSE)
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

\item{method}{the prediction method to use. see 'Details'}

\item{tp}{the prediction horizon (how far ahead to forecast)}

\item{num_neighbors}{the number of nearest neighbors to use. Note that the 
default value will change depending on the method selected. (any of "e+1", 
"E+1", "e + 1", "E + 1" will set this parameter to E+1 for each run.)}

\item{columns}{either a vector with the columns to use (indices or names), 
or a list of such columns}

\item{target_column}{the index (or name) of the column to forecast}

\item{stats_only}{specify whether to output just the forecast statistics or 
to include the raw predictions for each run}

\item{first_column_time}{indicates whether the first column of the given 
block is a time column (and therefore excluded when building the library)}

\item{exclusion_radius}{excludes vectors from the search space of nearest 
neighbors if their *time index* is within exclusion_radius (NULL turns 
this option off)}

\item{epsilon}{Not implemented}

\item{theta}{the nonlinear tuning parameter (theta is only relevant if 
method == "s-map")}

\item{silent}{prevents warning messages from being printed to the R console}

\item{save_smap_coefficients}{specifies whether to include the s_map 
coefficients with the output}
}
\value{
A data.frame with components for the parameters and forecast 
  statistics:
\tabular{ll}{
  cols \tab embedding\cr
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
  model_output \tab data.frame with columns for the time index, 
    observations, predictions, and estimated prediction variance
    (if \code{stats_only == FALSE})\cr
}
If "s-map" is the method, then the same, but with additional columns:
\tabular{ll}{
  theta \tab the nonlinear tuning parameter\cr
  smap_coefficients \tab data.frame with columns for the s-map 
  coefficients (if \code{save_smap_coefficients == TRUE})\cr
  smap_coefficient_covariances \tab list of covariance matrices for 
  the s-map coefficients (if \code{save_smap_coefficients == TRUE})\cr
}
}
\description{
\code{\link{block_lnlp}} uses multiple time series given as input to generate 
  an attractor reconstruction, and then applies the simplex projection or 
  s-map algorithm to make forecasts. This method generalizes the 
  \code{\link{simplex}} and \code{\link{s_map}} routines, and allows for 
  "mixed" embeddings, where multiple time series can be used as different 
  dimensions of an attractor reconstruction.
}
\details{
The default parameters are set so that passing a vector as the only argument
  will use that vector to predict itself one time step ahead. If a matrix or 
  data.frame is given as the only argument, the first column will be 
  predicted (one time step ahead), using the remaining columns as the 
  embedding. If the first column is not a time vector, 1:NROW will be
  used as time values.

\code{norm = 2} (only option currently available) uses the "L2 norm",
Euclidean distance:
  \deqn{distance(a,b) := \sqrt{\sum_i{(a_i - b_i)^2}}
  }{distance(a, b) := \sqrt(\sum(a_i - b_i)^2)}

method "simplex" (default) uses the simplex projection forecasting algorithm

method "s-map" uses the s-map forecasting algorithm
}
\examples{
block <- block_3sp
block_lnlp(block[,2:4])

block <- block_3sp
block_lnlp(block[,1:4], first_column_time = TRUE)

block <- block_3sp
block_lnlp(block, target_column = "x_t", columns = c("y_t", "z_t"), first_column_time = TRUE)

block <- block_3sp
x_t_pred = block_lnlp(block, columns = c("x_t", "y_t"), first_column_time = TRUE,
stats_only = FALSE)

block <- block_3sp
x_t_pred = block_lnlp(block, method = "s-map", theta = 3, columns =
c("x_t", "y_t"), first_column_time = TRUE, stats_only = FALSE, save_smap_coefficients = TRUE)
}
