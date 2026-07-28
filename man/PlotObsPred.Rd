\name{PlotObsPred}
\alias{PlotObsPred}
\title{Plot observed versus predicted values}
\description{Base-graphics plot of Observations and Predictions against time.}
\usage{
PlotObsPred(df, dataName = "", E = NULL, Tp = NULL)
}
\arguments{
\item{df}{A data.frame with a time column plus Observations and Predictions.}
\item{dataName}{Optional label used in the plot title.}
\item{E}{Embedding dimension.}
\item{Tp}{Forecast interval (prediction horizon).}
}
\value{Invisibly \code{NULL}.}
\keyword{hplot}
