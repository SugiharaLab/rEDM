\name{compute_stats}
\alias{compute_stats}
\title{Compute performance metrics for predictions}
\arguments{
\item{observed}{a vector of the observed values}

\item{predicted}{a vector of corresponding predicted values}
}
\value{
A data.frame with components with various performance metrics:
\tabular{ll}{
  num_pred \tab number of predictions\cr
  rho \tab correlation coefficient between observations and predictions\cr
  mae \tab mean absolute error\cr
  rmse \tab root mean square error\cr
  perc \tab percent correct sign\cr
  p_val \tab p-value that rho is significantly greater than 0 using Fisher's 
}
}
\description{
Computes the rho, MAE, RMSE, perc, and p-val performance 
  metrics
}
\examples{
compute_stats(rnorm(100), rnorm(100))
}
