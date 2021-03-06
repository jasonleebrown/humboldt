% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{env2}
\alias{env2}
\title{Example Data: Bioclimatic data for Europe and surrounding areas}
\format{
An data frame with 23758 rows and 21 variables:
\describe{
  \item{x}{longitude}
  \item{y}{latittude}
  \item{var1}{any continious enviromental variable}
  \item{var2}{any continious enviromental variable}
}
}
\usage{
data(env2)
}
\description{
Example Data: Bioclimatic data for Europe and surrounding areas
}
\examples{

library(humboldt)

##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
data(env1)
## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
data(env2)

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

##load occurrence sites for the species at study area 1 (env1). Column names should be sp,x,y
data(sp1)
occ.sp1<-sp1

##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
data(sp2)
occ.sp2<-sp2
}
\keyword{datasets}
