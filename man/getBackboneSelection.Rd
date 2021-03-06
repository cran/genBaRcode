% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{getBackboneSelection}
\alias{getBackboneSelection}
\title{Predefined Barcode Backbone Sequences}
\usage{
getBackboneSelection(id = NULL)
}
\arguments{
\item{id}{an integer or character value in order to choose a specific backbone.}
}
\value{
a character string.
}
\description{
allows the user to choose between predefined backbone sequences. Excecution of the function without any parameter
value will display all available backbone sequences. The id parameter will accept the name of the backbone or the
rownumber of the shown selection.
}
\examples{
getBackboneSelection()
getBackboneSelection(2)
getBackboneSelection("BC32-Venus")
}
