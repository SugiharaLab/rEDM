\name{ComputeError}
\alias{ComputeError}
\title{Prediction skill statistics}
\description{Pearson rho, MAE, CAE and RMSE between observations and predictions.}
\usage{
ComputeError(obs, pred, digits = 6)
}
\arguments{
\item{obs}{Numeric vector of observations.}
\item{pred}{Numeric vector of predictions aligned with \code{obs}.}
\item{digits}{Rounding applied to each statistic.}
}
\value{A named list with elements rho, MAE, CAE and RMSE.}
\examples{
ComputeError(1:10 + 0, c(1:9 + 0, 8))
}
\keyword{regression}
