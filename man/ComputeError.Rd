\name{ComputeError}
\alias{ComputeError}
\title{Compute error}
\usage{
ComputeError(obs, pred)  
}
\arguments{
\item{obs}{vector of observations.}

\item{pred}{vector of predictions.}
}

\value{
  A name list with components:
  \tabular{ll}{
    rho   \tab Pearson correlation\cr
    MAE   \tab mean absolute error\cr
    RMSE  \tab root mean square error\cr
  }
}

\description{
\code{\link{ComputeError}} evaluates the Pearson correlation
coefficient, mean absolute error and root mean square error between two
numeric vectors.
}

\examples{
data(block_3sp)
smplx <- Simplex( dataFrame=block_3sp, lib="1 99", pred="105 190", E=3,
columns="x_t", target="x_t")
err <- ComputeError( smplx$Observations, smplx$Predictions )
}
