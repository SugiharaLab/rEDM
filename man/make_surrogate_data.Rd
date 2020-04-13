\name{make_surrogate_data}
\alias{make_surrogate_data}
\alias{make_surrogate_shuffle}
\alias{make_surrogate_ebisuzaki}
\alias{make_surrogate_seasonal}
\title{Generate surrogate data for permutation/randomization tests}
\usage{
make_surrogate_data(ts, method = c("random_shuffle", "ebisuzaki",
  "seasonal", "twin"), num_surr = 100, ...)

make_surrogate_shuffle(ts, num_surr = 100)

make_surrogate_ebisuzaki(ts, num_surr = 100)

make_surrogate_seasonal(ts, num_surr = 100, T_period = 12)
}
\arguments{
\item{ts}{the original time series}

\item{method}{which algorithm to use to generate surrogate data}

\item{num_surr}{the number of null surrogates to generate}

\item{T_period}{the period of seasonality for seasonal surrogates 
(ignored for other methods)}
}
\value{
A matrix where each column is a separate surrogate with the same 
  length as `ts`.
}
\description{
This is a wrapper function for generating surrogate time series 
  using several different null models.

\code{make_surrogate_shuffle()} creates surrogates by randomly 
  permuting the values of the original time series.

\code{make_surrogate_ebisuzaki()} creates surrogates by 
  randomizing the phases of a Fourier transform, preserving the power 
  spectra of the null surrogates

\code{make_surrogate_seasonal()} creates surrogates by 
  computing a mean seasonal trend of the specified period and shuffling the 
  residuals.
}
\examples{
data = make_surrogate_data(block_3sp$x_t)
}
