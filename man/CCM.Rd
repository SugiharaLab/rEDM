\name{CCM}
\alias{CCM}
\title{Convergent cross mapping}
\description{
  The state-space of a multivariate dynamical system (not a purely
  stochastic one) encodes coherent phase-space variable trajectories. If
  enough information is available, one can infer the presence or absence
  of cross-variable interactions associated with causal links between
  variables. \code{\link{CCM}} measures the extent to which states of
  variable Y can reliably estimate states of variable X. This can happen
  if X is causally influencing Y. 

  If cross-variable state predictability converges as more state-space
  information is provided, this indicates a causal link. \code{\link{CCM}}
  performs this cross-variable mapping using Simplex, with convergence
  assessed across a range of observational library sizes as described in
  \cite{Sugihara et al. 2012}.
}

\usage{
CCM(dataFrame = NULL, columns, target, E, Tp = 0, knn = 0, tau = -1,
    exclusionRadius = 0, libSizes, sample = 30, seed = NULL,
    embedded = FALSE, validLib = logical(0), includeData = FALSE,
    numProcess = 4, backend = "RANN", pathIn = "./", dataFile = "",
    pathOut = "./", predictFile = "", parameterList = FALSE, showPlot = FALSE)
}
\arguments{
\item{dataFrame}{A data.frame of input data. The first column must be a time
index or time values unless \code{noTime = TRUE}. The columns must be named.}
\item{columns}{Column name(s) to build the embedding: character vector
or space-separated string of columns to create the library. If
individual column names contain whitespace place names in a vector, or,
append ',' to the name.}
\item{target}{Target column name to predict.}
\item{E}{Embedding dimension.}
\item{Tp}{Forecast interval (prediction horizon).}
\item{knn}{Number of nearest neighbours. \code{0} selects the method default.}
\item{tau}{Embedding delay (negative selects past lags).}
\item{exclusionRadius}{Temporal (Theiler) exclusion radius around each prediction point.}
\item{libSizes}{Library sizes: string of 3 whitespace separated integer values
specifying the intial library size, the final library size,
and the library size increment. Can also be a list of strictly
 ncreasing library sizes.}
\item{sample}{Number of random subsamples per library size.}
\item{seed}{Random seed for reproducibility (not bit-comparable to pyEDM/NumPy).}
\item{embedded}{If \code{TRUE}, \code{columns} already form the embedding.}
\item{validLib}{Logical vector marking admissible library rows (or length 0 for all).}
\item{includeData}{Add per-library-size sample-variance columns to the result.}
\item{numProcess}{Number of worker processes for the parameter sweep or task grid.}
\item{backend}{Nearest-neighbour backend: \code{"RANN"} (default) or \code{"brute"}.}
\item{pathIn}{File path for input \code{dataFile}.}
\item{dataFile}{Input dataFile, .csv format. The first column must be a time
index or time values unless noTime is TRUE. The first row must be column names.}
\item{pathOut}{Output file path for \code{predictFile}}
\item{predictFile}{Output file name, .csv format.}
\item{parameterList}{Append named list of parameters/values to return.}
\item{showPlot}{If \code{TRUE}, draw a base-graphics plot of the result.}
}
\value{
  A data.frame with 3 columns. The first column is \code{LibSize}
  specifying the subsampled library size.  Columns 2 and 3 report
  Pearson correlation coefficients for the prediction of X from Y, and
  Y from X.

  If \code{includeData = TRUE} and \code{parameterList = TRUE} a
  named list "parameters" is added.
}

\references{Sugihara G., May R., Ye H., Hsieh C., Deyle E., Fogarty M., Munch S., 2012. Detecting Causality in Complex Ecosystems. Science 338:496-500.
}

\examples{
data(sardine_anchovy_sst)
df = CCM( dataFrame = sardine_anchovy_sst, E = 3, Tp = 0, columns = "anchovy",
target = "np_sst", libSizes = "10 70 10", sample = 100 )
}
\keyword{nonlinear,ts}
