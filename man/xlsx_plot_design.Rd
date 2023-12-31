% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_gene.R
\name{xlsx_plot_design}
\alias{xlsx_plot_design}
\title{Designs the WB that will be displayed on a XLSX file afterwards}
\usage{
xlsx_plot_design(plots_to_xlsx)
}
\arguments{
\item{plots_to_xlsx}{Nested list of pileup plots. The main list is named after the mutated genes. Inside of them there is the gene plot and all mutation plots belonging to that gene.}
}
\value{
WB with all the arranged plots in it
}
\description{
Depending on the size of each plot, this function arranges the plot in the different cells of the WB
}
