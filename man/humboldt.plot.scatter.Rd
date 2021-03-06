% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humboldt66.R
\name{humboldt.plot.scatter}
\alias{humboldt.plot.scatter}
\title{Density scatter plots with histograms on each axis}
\usage{
humboldt.plot.scatter(x = NA, color.ramp = 1, lhist = 20, ...)
}
\arguments{
\item{x}{A data frame with two columns of data for plotting}

\item{color.ramp}{An integer from 1-6 depict Humbodlt's six coloramps: 1= rgb 2= plasma, 3=viridis, 4=sunset, 5= rainbow, 6= greyscale. for visual depiction, see: https://github.com/jasonleebrown/humboldt/blob/master/HumboldtInputExp.pdf}

\item{lhist}{number of bins in histrograms}
}
\value{
A simple function to plot 2 dimensions of a data frame. Warmer colors represent higher densities. On each axis the scatter plot histograms are plotted.
}
\description{
Density scatter plots with histograms on each axis
}
\examples{

library(humboldt)
##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
env1<-read.delim("env1.txt",h=T,sep="\t")

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
humboldt.plot.scatter(env1[,3:4], xlab="Bio1", ylab="Bio2",main="environment")
}
