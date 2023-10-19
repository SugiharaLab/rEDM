\docType{package}
\name{rEDM}
\alias{EDM}
\alias{rEDM-package}
\title{Empirical dynamic modeling}
\description{
  \pkg{rEDM} provides tools for data-driven time series analyses. It is
  based on reconstructing multivariate state space
  representations from uni or multivariate time series, then projecting
  state changes using various metrics applied to nearest neighbors.
  
  \pkg{rEDM} is a \pkg{Rcpp} interface to the
  \href{https://github.com/SugiharaLab/cppEDM}{cppEDM} library of
  Empirical Dynamic Modeling tools. Functionality includes:
  \itemize{
    \item Simplex projection (Sugihara and May 1990)
    \item Sequential Locally Weighted Global Linear Maps (S-map) (Sugihara 1994)
    \item Multivariate embeddings (Dixon et. al. 1999)
    \item Convergent cross mapping (Sugihara et. al. 2012)
    \item Multiview embedding (Ye and Sugihara 2016)
  }
}
\references{
  Sugihara G. and May R. 1990. Nonlinear forecasting as a way of
  distinguishing chaos from measurement error in time series.
  Nature, 344:734-741.

  Sugihara G. 1994. Nonlinear forecasting for the classification of
  natural time series. Philosophical Transactions: Physical Sciences
  and Engineering, 348 (1688) : 477-495.

  Dixon, P. A., M. Milicich, and G. Sugihara, 1999. Episodic
  fluctuations in larval supply. Science 283:1528-1530.

  Sugihara G., May R., Ye H., Hsieh C., Deyle E., Fogarty M.,
  Munch S., 2012. Detecting Causality in Complex Ecosystems.
  Science 338:496-500.

  Ye H., and G. Sugihara, 2016. Information leverage in
  interconnected ecosystems: Overcoming the curse of dimensionality.
  Science 353:922-925.
}
\details{

\strong{Main Functions}: 
  \itemize{
    \item \code{\link{Simplex}} - simplex projection
    \item \code{\link{SMap}} - S-map projection
    \item \code{\link{CCM}} - convergent cross mapping
    \item \code{\link{Multiview}} - multiview forecasting
  }
\strong{Helper Functions}: 
  \itemize{
    \item \code{\link{Embed}} - time delay embedding
    \item \code{\link{ComputeError}} - forecast skill metrics
    \item \code{\link{EmbedDimension}} - optimal embedding dimension
    \item \code{\link{PredictInterval}} - optimal prediction interval
    \item \code{\link{PredictNonlinear}} - evaluate nonlinearity
  }
}
\author{
\strong{Maintainer}: Joseph Park

\strong{Authors}: Joseph Park, Cameron Smith, Ethan Deyle, Erik
Saberski, George Sugihara

% \strong{Contributors}: 
}
\keyword{package}
