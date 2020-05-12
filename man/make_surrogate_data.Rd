\name{make_surrogate_data}
\alias{make_surrogate_data}
\title{Generate surrogate data for permutation/randomization tests}
\usage{
make_surrogate_data(ts, method = c("random_shuffle", "ebisuzaki",
  "seasonal"), num_surr = 100, T_period = 1, alpha = 0)
}
\arguments{
\item{ts}{the original time series}

\item{method}{which algorithm to use to generate surrogate data}

\item{num_surr}{the number of null surrogates to generate}

\item{T_period}{the period of seasonality for seasonal surrogates 
(ignored for other methods)}

\item{alpha}{standard deviation of seasonal cycle deviates.}
}
\value{
A matrix where each column is a separate surrogate with the same 
  length as `ts`.
}
\description{
This is a wrapper function for generating surrogate time series 
  using several different null models.
}
\examples{
data = make_surrogate_data(block_3sp$x_t)
}
