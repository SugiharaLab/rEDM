\name{SurrogateData}
\alias{SurrogateData}
\title{Generate surrogate data for permutation/randomization tests}
\usage{
SurrogateData( ts, method = c("random_shuffle", "ebisuzaki",
"seasonal"), num_surr = 100, T_period = 1, alpha = 0 )
}

\arguments{
\item{ts}{the original time series}

\item{method}{which algorithm to use to generate surrogate data}

\item{num_surr}{the number of null surrogates to generate}

\item{T_period}{the period of seasonality for seasonal surrogates
  (ignored for other methods)}

\item{alpha}{additive noise factor: N(0,alpha)}
}

\value{
  A matrix where each column is a separate surrogate with the same
  length as \code{ts}.
}

\description{
\code{SurrogateData} generates surrogate data under several different 
null models.
}

\details{
Method "random_shuffle" creates surrogates by randomly permuting the values 
of the original time series.

Method "Ebisuzaki" creates surrogates by randomizing the phases of a Fourier 
transform, preserving the power spectra of the null surrogates.

Method "seasonal" creates surrogates by computing a mean seasonal trend of 
the specified period and shuffling the residuals.  It is presumed that
the seasonal trend can be exracted with a smoothing spline.  Additive
Gaussian noise is included according to N(0,alpha). 
}

\examples{
data("block_3sp")
ts <- block_3sp$x_t
SurrogateData(ts, method = "ebisuzaki")
}
