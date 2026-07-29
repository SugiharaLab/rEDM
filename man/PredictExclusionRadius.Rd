\name{PredictExclusionRadius}
\alias{PredictExclusionRadius}
\title{Prediction skill versus exclusion radius}
\description{Simplex prediction skill (rho) swept over the temporal exclusion radius.}
\usage{
PredictExclusionRadius(dataFrame = NULL, columns, target, lib, pred,
                       exclusionRadius = NULL, E = 1, Tp = 1, tau = -1,
                       embedded = FALSE, validLib = logical(0),
                       noTime = FALSE, ignoreNan = TRUE, numProcess = 4,
                       backend = "RANN", pathIn = "./", dataFile = "",
                       pathOut = "./", predictFile = "",
                       parameterList = FALSE, showPlot = FALSE)
}
\arguments{
\item{dataFrame}{A data.frame of input data. The first column must be a time
index or time values unless \code{noTime = TRUE}.}
\item{columns}{Column name(s) to build the embedding: character vector or space-separated string.}
\item{target}{Target column name to predict.}
\item{lib}{Library (training) index range as (start end) pairs.}
\item{pred}{Prediction index range as (start end) pairs.}
\item{exclusionRadius}{Temporal (Theiler) exclusion radius around each prediction point.}
\item{E}{Embedding dimension.}
\item{Tp}{Forecast interval (prediction horizon).}
\item{tau}{Embedding delay (negative selects past lags).}
\item{embedded}{If \code{TRUE}, \code{columns} already form the embedding.}
\item{validLib}{Logical vector marking admissible library rows (or length 0 for all).}
\item{noTime}{If \code{TRUE}, synthesise a 1..N time index instead of using column 1.}
\item{ignoreNan}{Remove rows with NaN in the embedding from the library and prediction sets.}
\item{numProcess}{Number of worker processes for the parameter sweep or task grid.}
\item{backend}{Nearest-neighbour backend: \code{"RANN"} (default) or \code{"brute"}.}
\item{pathIn}{File path for input \code{dataFile}.}
\item{dataFile}{Input dataFile, .csv format.}
\item{pathOut}{Output file path for \code{predictFile}}
\item{predictFile}{Output file name, .csv format.}
\item{parameterList}{Append named list of parameters/values to return.}
\item{showPlot}{If \code{TRUE}, draw a base-graphics plot of the result.}
}
\value{A data.frame with columns ExclusionRadius and rho.}

\keyword{nonlinear}
