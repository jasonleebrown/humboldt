% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humboldt60.R
\name{humboldt.equivalency.test}
\alias{humboldt.equivalency.test}
\title{Niche equivalency test}
\usage{
humboldt.equivalency.test(z1, z2, rep = 100, correct.env = T,
  kern.smooth = 1, nae = "YES", thresh.espace.z = 0.001,
  run.silent.equ = F, ncores = 1)
}
\arguments{
\item{z1}{a grid of the density of a species occurrence in environmental space output from humboldt.grid.espace}

\item{z2}{a grid of the density of a second species occurrence in environmental space output from humboldt.grid.espace that you wish to compared to z1}

\item{rep}{is the number of iterations. Values higher than 100 is recommend for final analysis.}

\item{kern.smooth}{scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodal espace occupancy can be somewhat common when a species occupies an extreme aspect of habitat or when espace is not broadly accessible in both dimensions of espace (PCs 1 & 2)}

\item{nae}{do you include non-analogous environments in the niche similarity measurement? If nae="NO" (use capital letters), then non-analogous environments will be removed from both input environments during overlap measurement and only environments present in both will be used to measure overlap. If nae="YES" then no change will be made to input z1 and z2. Note: this is separate from trimming non-analogous environments from your input dataset (as done by humboldt.g2 specified by parameter non.analogous.environments). This parameter physically removes non-analogous environments from datasets pre-niche similarity measurement. Technically the removal of non-analogous environments via either way should result in similar overlap measurements (though they may not be identical). This because removing NAE from the dataset prior to gridding environments will resulting only non-analogous environments to be gridded (and typically finer grain applied to each grid cell). Whereas removing them only via this parameter (nae), which only removes non-analogous in the gridded environmental space for use in overlap measurements--- all the input environmental space is gridded (likely increasing the environmental space per gridded cell). A second cause of differences in values can result from rescaling of espace values during niche-overlap measurements so that the sum of the landscape equals one. If occupied non-analogous environmental are numerous in one of the datasets, this can theoretically cause overlap values to decrease in analogous environments (vs. nae) because differences in core niches are rescaled to 1 in both scenarios. The rescaling among fewer cells increases the values applied to highly suitable areas and, if not equivalently scaled in both datasets, differences among niches could increase, resulting a smaller overlap in non-analogous environments (again values should be similar). If you remove non-analogous environments in humboldt.g2e, I also suggest that you use this function (as it can remove any slight anomalies caused by gridding environments in humboldt.grid.clim due to the binning of values in the RxR grid).}

\item{thresh.espace.z}{this parameter is an experimental parameter and controls the level at which values below the kernel density z values are removed for creating areas of analogous environmental space. Higher values will increase value from which the low-density areas are removed from the environmental space of z1 and z2.  Basically values above this are retained and values below are removed. Default=0.001}

\item{ncores}{number of CPUs to use for tests. If unsure on the number of cores and want to use all but 1 CPU, input ncores="All"}
}
\value{
runs a modified niche equivalency test(see Warren et al 2008, Broenniman et al. 2012, but note these analyses are not identical) based on two species occurrence density grids and compares the observed niche similarity between z1 and z2 (created by humboldt.grid.espace) to overlap between random niches z1.sim and z2.sim. The z1.sim and z2.sim are built from random reallocations of occurrences of z1 and z2 (see 'sample' parameter). A significant value states that the two datasets are NOT statistically equivalent and rejects the NULL hypothesis that species niches are equivalent. \cr 
\cr
Output: $sim= simulation niche similarity values, $obs.D= Schoener's D value for two observed datasets, $obs.I= Hellinger's I value for observed datasets, $p.D=one-tail p-value of Schoener's D values (simulation vs. observed), $p.I=one-tail p-value of Hellinger's I value (simulation vs. observed)
}
\description{
Niche equivalency test
}
\examples{
library(humboldt)

##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
env1<-read.delim("env1.txt",h=T,sep="\\t")

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
env2<-read.delim("env2.txt",h=T,sep="\\t") 

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\\t"))

##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\\t"))

##convert geographic space to espace
zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim= T, e.var=c(3:21),  col.env = e.var, trim.buffer.sp1 = 200, trim.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)

##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
scores.env1<-zz$scores.env1[1:2]
scores.env2<-zz$scores.env2[1:2]
scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
scores.sp1<-zz$scores.sp1[1:2]
scores.sp2<-zz$scores.sp2[1:2]

## run create a grid of Environmental Space Function
z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)

## mesure niche equivalency
niche.equiv<- humboldt.equivalency.test(z1,z2,rep=10,kern.smooth=1, ncores=2)
}
\seealso{
\code{humboldt.niche.similarity, humboldt.doitall, humboldt.plot.density, humboldt.plot.histrogram} which use or depend on outputs of this function
}
