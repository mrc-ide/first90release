% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inputs.R
\name{prepare_inputs}
\alias{prepare_inputs}
\title{Create model inputs from Spectrum PJNZ}
\usage{
prepare_inputs(pjnzlist)
}
\arguments{
\item{pjnzlist}{a vector of PJNZ file names to aggregate}
}
\description{
Create model inputs from Spectrum PJNZ
}
\details{
The aggregation makes a number of assumptions:
\itemize{
\item Progression parameters are the same in all files, and values frome the first file are used.
\item Special populations ART eligibility is the same in all files.
}
}
\examples{

pjnzlist <- list.files("~/Documents/Data/Spectrum files/2018 final/SSA/",
  "CotedIvoire.*PJNZ$", full.names=TRUE, ignore.case=TRUE)
pjnzlist <- "~/Documents/Data/Spectrum files/2018 final/SSA/Malawi_2018_version_8.PJNZ"

}
