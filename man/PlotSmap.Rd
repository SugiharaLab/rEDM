\name{PlotSmap}
\alias{PlotSmap}
\title{Plot S-map predictions and coefficients}
\description{Base-graphics multi-panel plot of S-map predictions and coefficient series.}
\usage{
PlotSmap(SmapList, dataName = "", E = NULL, Tp = NULL)
}
\arguments{
\item{SmapList}{A list returned by \code{SMap} (predictions and coefficients).}
\item{dataName}{Optional label used in the plot title.}
\item{E}{Embedding dimension.}
\item{Tp}{Forecast interval (prediction horizon).}
}
\value{Invisibly \code{NULL}.}
\keyword{hplot}
