% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plotTimeSeries}
\alias{plotTimeSeries}
\title{Plotting Time Series Data}
\usage{
plotTimeSeries(
  ov_dat,
  colr = NULL,
  tp = NULL,
  x_label = "time",
  y_label = "contribution"
)
}
\arguments{
\item{ov_dat}{a numeric matrix consisting of all time points as columns and all barcode sequences as rows and the corresponding read counts as numerical values (see function \code{generateTimeSeriesData()}).}

\item{colr}{a vector of character strings identifying a certain color palette.}

\item{tp}{a numeric vector containing the time points of measurement (in case of unequally distributed time points).}

\item{x_label}{a character string providing the x-axis label.}

\item{y_label}{a character string providing the y-axis label.}
}
\value{
a ggplot2 object.
}
\description{
Uses the result of the generateTimeSeriesData function as inout and generates a visualisation of the clonal
contributions over a number of given time points (similar to a stacked barplot).
}
\examples{
ov_dat <- matrix(round(runif(1:100, min = 0, max = 1000)), ncol = 5)
rownames(ov_dat) <- paste("barcode", 1:20)
plotTimeSeries(ov_dat)
}
