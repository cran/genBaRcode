% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{asBCdat}
\alias{asBCdat}
\title{Data Type Conversion}
\usage{
asBCdat(dat, label = "empty", BC_backbone = "none", resDir = getwd())
}
\arguments{
\item{dat}{a data.frame object with two columns containing read counts and barcode sequences.}

\item{label}{a optional character string used as label.}

\item{BC_backbone}{a optional character string, describing the barcode backbone structure.}

\item{resDir}{a optional character string, identifying the path to the results directory, default is current working directory.}
}
\value{
a BCdat object.
}
\description{
Converts a data.frame into a BCdat object.
}
