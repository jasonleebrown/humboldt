% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humboldt66.R
\name{humboldt.accessible.e.distance}
\alias{humboldt.accessible.e.distance}
\title{Measure the effects of accessible environment distance in divergence test}
\usage{
humboldt.accessible.e.distance(
  env1 = env1,
  env2 = env2,
  sp1 = occ.sp1,
  sp2 = occ.sp2,
  env.trim = T,
  env.trim.type = "RADIUS",
  e.var = c(3:21),
  rarefy.dist = 50,
  rarefy.units = "km",
  env.reso = 0.41666669,
  kern.smooth = 1,
  R = 100,
  rep = 100,
  min.d = 5,
  max.d = 1000,
  n.interval = 10,
  nae = "YES",
  thresh.espace.z = 0.001,
  correct.env = T,
  force.equal.sample = T,
  ncores = 1,
  run.silent = F
)
}
\arguments{
\item{env1}{environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice.}

\item{env2}{environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice.}

\item{sp1}{occurrence sites for the species/population 1 at study area 1 (env1). Column names should be 'sp', 'x','y'}

\item{sp2}{occurrence sites for the species/population 2 at study area 2 (env2). Column names should be 'sp', 'x', 'y'}

\item{env.trim}{Trim extent of environmental data in geographic space. Necessary for comparing if species have diffrent access to habitats.  If this env.trim=T, non-accesible environments will be removed.  See associated parameters "env.trim.type", "trim.buffer.sp1", "trim.buffer.sp2"}

\item{env.trim.type}{This parameter is only used if env.trim=TRUE.  There options exist, trim enviromental data by: a buffered minimum-convex-polygon (env.trim.type="MCP"), a buffer around all occurence localties (env.trim.type="RADIUS"), and an advanced input mask option (env.trim.type="MASK") that allows a user to input a shapefile. For env.trim.type="RADIUS" & "MCP", the parameters 'trim.buffer.sp1' and 'trim.buffer.sp2' specifices the buffer distance used to trim accessible enviroments in km.  Default= "MCP"}

\item{e.var}{selection of variables to include in all of the analyses of E-space. This is a separate parameter than col.env, but must contain all variables included in col.env.   Note that it can include more variables than those in col.env, as long as those in col.env are also included.}

\item{rarefy.dist}{remove occurrences closer than a minimum distance to each other (this function uses the humboldt.occ.rarefy function). Values need to be in km[recommended] or decimal degrees. See associated parameter rarefy.units. Note: rarefy.dist=0 will remove no occurrences}

\item{rarefy.units}{the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees}

\item{env.reso}{the resolution of the input environmental data grid in decimal degrees}

\item{kern.smooth}{scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodal espace occupancy can be somewhat common when a species occupies an extreme aspect of habitat or when espace is not broadly accessible in both dimensions of espace (PCs 1 & 2)}

\item{R}{resolution of grid in environmental space (RxR)}

\item{rep}{the number of iterations for the equivalence and background statistics. Values higher than 200 are recommend for final analysis}

\item{min.d}{the minimum distance to sample accessible environmental space}

\item{max.d}{the maximum distance to sample accessible environmental space}

\item{n.interval}{the number of distance intervals to measure acessible test accessible environmental space}

\item{nae}{do you include non-analogous environments in the niche similarity measurement? If nae="NO" (use captial letters), then non-analogous environments will be removed from both input environments during overlap measurement and only environments present in both datasets will be used to measure overlap. If nae="YES", then no change will be made to input z1 and z2. Note: this is separate from trimming non-analogous environments from your input dataset (as done by humboldt.g2 specified by parameter non.analogous.environments). This parameter physically removes non-analogous environments from datasets ONLY before the niche similarity measurement. Technically the removal of non-analogous environments via either way should result in similar overlap measurements (though they may not be identical). This because removing NAE from the dataset prior to gridding environments will resulting only non-analogous environments to be gridded (and typically finer grain applied to each grid cell). Whereas removing them only via this parameter (nae), which only removes non-analogous in the gridded environmental space for use in overlap measurements--- all the input environmental space is gridded (likely increasing the environmental space per gridded cell). A second cause of differences in values can result from rescaling of espace values during niche-overlap measurements so that the sum of the landscape equals one. If occupied non-analogous environmental are numerous in one of the datasets, this can theoretically cause overlap values to decrease in analogous environments (vs. nae) because differences in core niches are rescaled to 1 in both scenarios. The rescaling among fewer cells increases the values applied to highly suitable areas and, if not equivalently scaled in both datasets, differences among niches could increase, resulting a smaller overlap in non-analogous environments (again values should be similar). If you remove non-analogous environments in humboldt.g2e, I also suggest that you use this function (as it can remove any slight anomalies caused by gridding environments in humboldt.grid.clim due to the binning of values in the RxR grid).}

\item{thresh.espace.z}{this parameter is an experimental parameter and controls the level at which values below the kernel density z values are removed for creating areas of analogous environmental space. Higher values will increase value from which the low-density areas are removed from the environmental space of z1 and z2.  Basical values above this are retained and values below are removed. Default=0.001}

\item{correct.env}{if correct.env=T, the analysis corrects occurrence densities of each species by the prevalence of the environments in their range. If correct.env=F, the overlap measure does not correct occurrence densities of each species by the prevalence of the environments in their range.Default=T}

\item{ncores}{number of CPUs to use for tests. If unsure on the number of cores and want to use all but 1 CPU, input ncores="All"}

\item{run.silent}{if run.silent=T, texts boxes displaying 'sampling', 'rarefying', 'equivalence statistic', 'background statistic' progress will not be displayed}

\item{reduce.env}{the format to trim environmental space so that it is shared. If reduce.env=1, the second input environment (env2) will be trim the match the first input (env1). If reduce.env=2, both input environments trimmed so that extents of both are identical (the lower maximum value observed in env1 and env2 and the higher minimum value observed in env1 and env2 will be used to trim environmental space for each PC/environmental variable) If reduce.env=0, you will skip trimming environmental space}

\item{trim.buffer.sp1}{buffer distance (in km) for trimming available environmental space for sp1}

\item{trim.buffer.sp2}{buffer distance (in km) for trimming available environmental space for sp2}

\item{col.env}{if reductype="STANDARD", then parameter specifies the number of columns to trim environmental space on. This can be any number of columns. This can be a subset or all of the enviroment layers input.}
}
\value{
This tool runs the equivalence and background statistics at several intervals spanning a minimum and maximum distance.  This tool tests the role of the input accsisble distance values in your equivalence and backround statistics and helps test the robustness of your results and their sensitivty to this value.
}
\description{
Measure the effects of accessible environment distance in divergence test
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

##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y' 
occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))

##convert geographic space to espace
zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim= T, e.var=c(3:21),  col.env = e.var, trim.buffer.sp1 = 200, trim.buffer.sp2 = 200, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, R = 100, run.silent = F)

##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
env1<-read.delim("env1.txt",h=T,sep="\t")

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
env2<-read.delim("env2.txt",h=T,sep="\t") 

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))

##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))

##perform accessible distance test
humboldt.accessible.e.distance(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, env.trim= T, env.trim.type="RADIUS", e.var=c(3:21), rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, R = 100, rep=100, min.d=5, max.d=1000, n.interval=10, nae = "YES", thresh.espace.z = 0.001, correct.env = T, force.equal.sample=T, ncores=2, run.silent = F)
}
\seealso{
\code{humboldt.sample.spp,humboldt.g2e, humboldt.equivalence.stat, humboldt.background.stat, humboldt.niche.similarity, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function
}
