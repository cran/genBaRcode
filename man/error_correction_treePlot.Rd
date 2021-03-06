% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{error_correction_treePlot}
\alias{error_correction_treePlot}
\title{Tree Plot}
\usage{
error_correction_treePlot(edges, vertices)
}
\arguments{
\item{edges}{a data frame containing edge definitions by two columns calles "from" and "to". Such data frame will be returned by the error_correction function with the EC_analysis parameter set to TRUE.}

\item{vertices}{a data frame with at least one column containing a list of nodes (also returned by the error_correction function with the EC_analysis parameter set to TRUE)}
}
\value{
a ggplot2 object.
}
\description{
creates a Tree Plot visualising of the barcode clustering as part of the error correction process.
}
