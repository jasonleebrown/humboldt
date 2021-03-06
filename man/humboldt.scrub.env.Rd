% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humboldt66.R
\name{humboldt.scrub.env}
\alias{humboldt.scrub.env}
\title{Scrub environment data}
\usage{
humboldt.scrub.env(in_data)
}
\arguments{
\item{imported}{table for analysis}
}
\value{
This tool remove NAs and converts all input data to numeric values for analysis in Humboldt.
}
\description{
Scrub environment data
}
\examples{
library(humboldt)
##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
env1<-read.delim("env1.txt",h=T,sep="\t")

env1<-humboldt.scrub.env(env1)
}
