% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humboldt66.R
\name{humboldt.top.env}
\alias{humboldt.top.env}
\title{Select top environmental variables for PCA}
\usage{
humboldt.top.env(
  env1,
  env2,
  sp1,
  sp2,
  rarefy.dist = 0,
  rarefy.units = "km",
  env.reso,
  learning.rt1 = 0.01,
  learning.rt2 = 0.01,
  e.var,
  pa.ratio = 4,
  steps1 = 50,
  steps2 = 50,
  method = "contrib",
  nvars.save = 5,
  contrib.greater = 5
)
}
\arguments{
\item{env1}{environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice}

\item{env2}{environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice}

\item{sp1}{occurrence sites for the species/population 1 at study area 1 (env1). Column names should be 'sp', 'x','y'}

\item{sp2}{occurrence sites for the species/population 2 at study area 2 (env2). Column names should be 'sp', 'x','y'}

\item{rarefy.dist}{removes occurrences within a minimum distance (specified here) to each other (this function uses the humboldt.occ.rarefy function). Values need to be in km[recommended] or decimal degrees. See associated parameter rarefy.units. Note: rarefy.dist=0 will remove no occurrences}

\item{rarefy.units}{the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees}

\item{env.reso}{the resolution of the input environmental data grid in decimal degrees}

\item{learning.rt1}{value from 0.01 to 0.001 for building SDM, start with 0.01 and if prompted, change to 0.001. the default value is 0.01}

\item{learning.rt2}{value from 0.01 to 0.001 for building SDM, start with 0.01 and if prompted, change to 0.001. The default value is 0.01}

\item{e.var}{Selection of variables to include in evaluation for each species}

\item{pa.ratio}{ratio of pseudoabsences to occurrence points, typically this is 4. The null value is 4}

\item{steps1}{numbers of trees to add at each cycle for modelling sp1. Start with 50 and if you run into problems gradually decrease, stopping at 1. The default value is 50}

\item{steps2}{numbers of trees to add at each cycle for modelling sp2. Start with 50 and if you run into problems gradually decrease, stopping at 1. The default value is 50}

\item{method}{this determines how important environmental variables are selected.  There are three options: "estimate", "contrib", "nvars". If method="estimate", the boosted regression tree algorithm will choose the number of variables to include by systematically removing variables until average change in the model exceeds the original standard error of deviance explained.  This is the most computationally intensive method. If method="contrib", variables above a relative influence value will be kept. See associated parameter 'contrib.greater'. If method="nvars", a fixed number of user specified variables will be kept. The kept variables are selected by their relative influence, selecting for the highest contributing variables. See associated parameter 'nvars.save'}

\item{nvars.save}{if method="nvars",this variable is required. It is the number of the top variables to save per species. The kept variables are selected by their relative influence in predicting the species distribution, selecting for the highest contributing variables. Often the total variables retained is lower due to identical variables select among both species. The default value is 5.  This value will be ignored if method="estimate" or "contrib"}

\item{contrib.greater}{if method="contrib", this variable is required. The kept variables are selected for their relative influence in predicting the species' distribution.  Here users select variables equal to or above an input model contribution value. The default value for this method is 5 (= variables with 5 percent or higher contribution to model of either species are kept). This value will be ignored if method="estimate" or "nvars"}
}
\value{
This function runs generalized boosted regression models (a machine learning SDM algorithm) to select top parameters for inclusion in PCA. This is important because you want the PC to reflect variables that are relevant to the species distribution. Alternatively you can run Maxent outside of R and manually curate the variables you include (also recommended).
}
\description{
Select top environmental variables for PCA
}
\examples{
library(humboldt)

##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
env1<-read.delim("env1.txt",h=T,sep="\t")

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
env2<-read.delim("env2.txt",h=T,sep="\t") 

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))

##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))

##perform modeling to determin imporant variables
reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=occ.sp1,sp2=occ.sp2,rarefy.dist=40, rarefy.units="km", env.reso=0.0833338,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)

##use new variables for env1 and evn2, use as you normally would do for env1/env2 (input above)

##for example, input into converted geographic space to espace
##zz<-humboldt.g2e(env1=reduc.vars$env1, env2=reduc.vars$env2....
}
