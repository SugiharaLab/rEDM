\name{SurrogateData}
\alias{SurrogateData}
\title{Surrogate time series for significance testing}
\description{Generate surrogate series (random_shuffle, ebisuzaki phase randomisation, or seasonal).}
\usage{
SurrogateData(ts, method = c("random_shuffle", "ebisuzaki", "seasonal"),
              num_surr = 100, T_period = 1, alpha = 0)
}
\arguments{
\item{ts}{Numeric time-series vector to generate surrogates for.}
\item{method}{Surrogate method: \code{"random_shuffle"}, \code{"ebisuzaki"}, or \code{"seasonal"}.}
\item{num_surr}{Number of surrogate series to generate.}
\item{T_period}{Seasonal period (seasonal method).}
\item{alpha}{Gaussian-noise standard deviation (seasonal method).}
}
\value{A numeric matrix with \code{length(ts)} rows and \code{num_surr} columns.}
\examples{
SurrogateData(sin(1:100), method = "ebisuzaki", num_surr = 5)
}
\keyword{datagen}
