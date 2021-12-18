\name{CCM}
\alias{CCM}
\title{Convergent cross mapping using simplex projection}
\usage{
CCM(pathIn = "./", dataFile = "", dataFrame = NULL, pathOut = "./", 
  predictFile = "", E = 0, Tp = 0, knn = 0, tau = -1,
  exclusionRadius = 0, columns = "", target = "", 
  libSizes = "", sample = 0, random = TRUE, replacement = FALSE, seed = 0, 
  includeData = FALSE, parameterList = FALSE, verbose = FALSE,
  showPlot = FALSE)  
}
\arguments{
\item{pathIn}{path to \code{dataFile}.}

\item{dataFile}{.csv format data file name. The first column must be a time
index or time values. The first row must be column names.}

\item{dataFrame}{input data.frame. The first column must be a time
index or time values. The columns must be named.}

\item{pathOut}{path for \code{predictFile} containing output predictions.}

\item{predictFile}{output file name.}

\item{E}{embedding dimension.}

\item{Tp}{prediction horizon (number of time column rows).}

\item{knn}{number of nearest neighbors. If knn=0, knn is set to E+1.}

\item{tau}{lag of time delay embedding specified as number of
time column rows.}

\item{exclusionRadius}{excludes vectors from the search space of nearest 
neighbors if their relative time index is within exclusionRadius.}

\item{columns}{string of whitespace separated column name(s) in the
input data used to create the library.}

\item{target}{column name in the input data used for prediction.}

\item{libSizes}{string of 3 whitespace separated integer values
  specifying the intial library size, the final library size,
  and the library size increment.}

\item{sample}{integer specifying the number of random samples to draw at
each library size evaluation.}

\item{random}{logical to specify random (\code{TRUE}) or sequential
  library sampling.}

\item{replacement}{logical to specify sampling with replacement.}

\item{seed}{integer specifying the random sampler seed.  If
  \code{seed=0} then a random seed is generated.}

\item{includeData}{logical to include statistics and predictions for
  every prediction in the ensemble.}

\item{parameterList}{logical to add list of invoked parameters.}

\item{verbose}{logical to produce additional console reporting.}

\item{showPlot}{logical to plot results.}
}

\value{
  A data.frame with 3 columns. The first column is \code{LibSize}
  specifying the subsampled library size.  Columns 2 and 3 report
  Pearson correlation coefficients for the prediction of X from Y, and
  Y from X.

  if \code{includeData = TRUE} a named list with the following data.frames
  data.frame \code{Combo_rho} columns:
  \tabular{ll}{
      LibMeans\tab CCM mean correlations for each library size\cr
      CCM1_PredictStat\tab Forward cross map prediction statistics\cr
      CCM1_Predictions\tab Forward cross map prediction values\cr
      CCM2_PredictStat\tab Reverse cross map prediction statistics\cr
      CCM2_Predictions\tab cross map prediction values\cr
  }

  If \code{includeData = TRUE} and \code{parameterList = TRUE} a
  named list "parameters" is added.
}

\references{Sugihara G., May R., Ye H., Hsieh C., Deyle E., Fogarty M., Munch S., 2012. Detecting Causality in Complex Ecosystems. Science 338:496-500.
}

\description{
  The state-space of a multivariate dynamical system (not a purely
  stochastic one) encodes coherent phase-space variable trajectories. If
  enough information is available, one can infer the presence or absence
  of cross-variable interactions associated with causal links between
  variables. \code{\link{CCM}} measures the extent to which states of
  variable Y can reliably estimate states of variable X. This happens
  only if X is causally influencing Y. 

  If cross-variable state predictability converges as more state-space
  information is provided, this indicates a causal link. \code{\link{CCM}}
  performs this cross-variable mapping using Simplex, with convergence
  assessed across a range of observational library sizes as described in
  \cite{Sugihara et al. 2012}.
}

\details{
  \code{\link{CCM}} computes the X:Y and Y:X cross-mappings in parallel
  using threads. 
}

\examples{
data(sardine_anchovy_sst)
df <- CCM( dataFrame=sardine_anchovy_sst, E=3, Tp=0, columns="anchovy",
target="np_sst", libSizes="10 70 10", sample=100 )

}
