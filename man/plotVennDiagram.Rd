% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plotVennDiagram}
\alias{plotVennDiagram}
\title{Plotting a VennDiagram}
\usage{
plotVennDiagram(
  BC_dat,
  alpha_value = 0.4,
  colrs = NA,
  border_color = NA,
  plot_title = "",
  legend_sort = NULL,
  annotationSize = 5
)
}
\arguments{
\item{BC_dat}{a list of BCdat objects.}

\item{alpha_value}{color transparency value [0-1].}

\item{colrs}{a character vector containing the desired colors, if NA the colors will be chosen automatically.}

\item{border_color}{a character value specifying the desired border color, if NA no border will be drawn.}

\item{plot_title}{a character value.}

\item{legend_sort}{a character or factor vector in case the order of legend items needs to be changed.}

\item{annotationSize}{an integer value specifying the venn diagramm internal text size.}
}
\value{
ggplot2 object.
}
\description{
plotVennDiagramm will create a Venn Diagram ans is based on the VennDiagram package.
It accepts a list of BCdat objects and will return a ggplot2 output object.
}
