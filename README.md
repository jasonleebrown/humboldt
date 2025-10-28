![Alt text](https://raw.githubusercontent.com/jasonleebrown/humboldt/master/humboldt.jpg?raw=true "Title") 


## Humboldt 2.0 is finally here!
<p> September 9, 2025</p>
    I am pleased to announce that Humbodlt 2.0 is released!  
<p> Main features:</p>
<p> -Updated to work in the Terra infrastructure</p>
<p>     (rgeos, raster, dismo, rgdal... all have been retired)</p>
<p> -Predictive modeling is now done with MaxEnt in humboldt.top.env function (vs. GBM in v1.0)</p>

<p>Known bugs:</p>
<p>    -humboldt.background.stat: multicore support is not currently functional  </p>

## Installing ‘humboldt’ in R
<p>Skip first line if 'devtools' is already installed.</p>

Open R (or R studio) and type or paste:
```markdown
install.packages("devtools")
library(devtools)
install_github("jasonleebrown/humboldt")
```

## Supporting Visual Guides
Having troubles or confused? Check out our [visual guide to parameter input](https://github.com/jasonleebrown/humboldt/blob/master/HumboldtInputExp.pdf) and [visual guide to interpreting outputs](https://github.com/jasonleebrown/humboldt/blob/master/HumboldtFigsExp.pdf).

## Contact 
[We can help you sort out issues, maybe](https://www.jasonleebrown.org/get-in-touch)

## Software citation and associated manuscript:
**Brown, J.L.** & **Carnaval, A.C.** (2019) A tale of two niches: methods, concepts and evolution. _Frontiers of Biogeography_, 11, e44158. doi:10.21425/F5FBG44158

## Input Data formats
To explore data format for input data, see:
```markdown
library(humboldt)

##for environment data
data(env1)

##for species data
data(sp1)
```
If needed, see guide below to convert raster GIS data for use as environment data in 'humboldt'

## Analyses in 'humboldt'
The new methods introduced in 'humboldt' translate several important theoretical advances into tests of niche divergence that allow researchers to more accurately estimate whether species have actually evolved different niches, or if they occupy different environmental spaces as the result of differences in life history, their biological interactors, or in the variety and configuration of accessible environments.  The foundation of these improvements is based on the underlying assumption that most species contemporary distributions are in non-equilibrium states and, because of this, the geographic manifestations of their niches (occupied, potential, and available fundamental) are dynamic through time.  The new methods provide several quantitative advances that characterize the accessible climates in both species’ distributions and the corresponding relationship between non-analogous and analogous climates.  Overall, the new methods improve the accuracy of niche similarity quantifications and corresponding statistical tests, consistently outperforming similar tests in correctly quantifying niche equivalence and divergence in simulated data with known truths. 

## The r package 'humboldt' provides the following analyses: 
### Potential Niche Truncation Index
Inferring the fundamental niche from a species’ occupied niche remains a great challenge; most studies of niche divergence overlook how well (or how badly) the occupied niches characterized from contemporary distributions potentially characterize a species’ fundamental niche. To provide a first step towards understanding this relationship, ‘humboldt’ provides a way to quantify the potential for a species’ occupied E-space to be truncated by the available E-space in its environment (see Brown & Carnaval 2019, Fig. 5). The larger the proportion of the occupied niche that is truncated in E-space, the higher the risk that the occupied niche may poorly reflect the fundamental niche. Based on the relationship between the species’ E-space and that available in adjacent habitats, we can assess the risk that the observed E-space is truncated, and how likely we are underestimating the species’ fundamental niches (see Brown & Carnaval 2019, Fig. 5). Here we introduce a new quantitative method to measure this: the Potential Niche Truncation Index (PNTI). It describes the amount of observed E-space of the species that is truncated by the available E-space. Specifically, it is a measurement of the overlap between the 5% kernel density isopleth of the species’ E-space and the 10% kernel density isopleth of accessible environment E-space. The PNTI is the portion of the species’ isopleth that falls outside of the environmental isopleth. This value physically measures how much of the perimeter of the species' E-space abuts, overlaps, or is outside the margins of the environment’s E-space. If the value is large, there is moderate risk (PNTI= 0.15-0.3) or high risk (PNTI>0.3) that the measured occupied niche does not reflect the species' fundamental niche due to niche truncation driven by limited available E-space.

### Example 
```markdown
library(humboldt)

## load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
data(env1)

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
data(env2)

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

## load occurrence sites for the species at study area 1 (env1). Column names should be sp,x,y
data(sp1)

## load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
data(sp2)

## its highly recommended that you using the function "humboldt.top.env" to select only the important environmental variables in humboldt.doitall. This step can be skipped. If you downloaded tons of environmental data, you should use this step.  

reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=sp1,sp2=sp2,rarefy.dist=50, rarefy.units="km", env.reso=0.416669,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)

## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e<-ncol(reduc.vars$env1)

## convert geographic space to espace for measuring pnt.index
zz<-humboldt.g2e(env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, reduce.env = 0, reductype = "PCA", non.analogous.environments = "YES", env.trim= T, e.var=c(3:num.var.e),  col.env = e.var, trim.buffer.sp1 = 200, trim.buffer.sp2 = 200, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, R = 100, run.silent = F)

## store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
scores.env1<-zz$scores.env1[1:2]
scores.env2<-zz$scores.env2[1:2]
scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
scores.sp1<-zz$scores.sp1[1:2]
scores.sp2<-zz$scores.sp2[1:2]

## estimate the Potential Niche Truncation Index
pnt1<- humboldt.pnt.index(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
pnt2<- humboldt.pnt.index(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
```

## Niche Overlap and Niche Divergence Tests
In our recent paper (Brown & Carnaval 2009), we then proposed a Niche Divergence Test and a Niche Overlap Test, which allows assessment of whether differences between species emerge from true niche divergences. The new methods improve accuracy of niche similarity and associated tests – consistently outperforming other tests. Our methods characterize the relationships between non-analogous and analogous climates in the species’ distributions, something not available previously.  These improvements allow assessment of whether the different environmental spaces occupied by two taxa emerge from true niche evolution, as opposed to differences in life history and biological interactors, or differences in the variety and configuration of environments accessible to them. 
### Example 1 - using provided datasets
```markdown
library(humboldt)

## load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
data(env1)

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
data(env2)

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

## load occurrence sites for the species at study area 1 (env1). Column names should be sp,x,y
data(sp1)

## load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
data(sp2)

## its highly recommended that you using the function "humboldt.top.env" to select only the important environmental variables in humboldt.doitall. 
## This step can be skipped. If you downloaded tons of environmental data, you should use this step.  If you skip this step, input env1/env2 in place of reduc.vars$env1/reduc.vars$env2 
reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=sp1,sp2=sp2,rarefy.dist=50, rarefy.units="km", env.reso=0.416669,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)

## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e<-ncol(reduc.vars$env1)

## run it first with full environmental for background tests and equivalence statistic (total equivalence or divergence in current distributions)
full<-humboldt.doitall(inname="full_extent", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=0, reductype="PCA", non.analogous.environments="YES", correct.env=T, env.trim=T,  env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, pcx=1, pcy=2, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="YES",thresh.espace.z=0.0001, p.overlap=T, p.boxplot=F, p.scatter=F, run.silent=F, ncores=2)

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche divergence. For comparing results, change only the following model parameters: reduce.env, non.analogous.environmental, env.trim, nae
shared_ae<-humboldt.doitall(inname="shared_espace_ae", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=2, reductype="PCA", non.analogous.environments="NO", correct.env=T, env.trim=T, env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, pcx=1,pcy=2, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="NO",thresh.espace.z=0.0001, p.overlap=T, p.boxplot=F, p.scatter=T,run.silent=F, ncores=2)
```

### Example 2 - typical workflow
see below for help formating raster/environment data. 
```markdown
library(humboldt)
## load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
## in this example all input datasets are tab-delimited text files, if using '.csv' files change the parameters below for import steps from 'sep="\t"' to 'sep=","' 
env1<-read.delim("env1.txt",h=T,sep="\t")

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
env2<-read.delim("env2.txt",h=T,sep="\t") 

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

## load occurrence sites for the species at study area 1 (env1). Column names should be sp,x,y
sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))

## load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y 
sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))

## its highly recommended that you using the function "humboldt.top.env" to select only the important environmental variables in humboldt.doitall. 
## This step can be skipped. If you downloaded tons of environmental data, you should use this step.  If you skip this step, input env1/env2 in place of reduc.vars$env1/reduc.vars$env2 
reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=sp1,sp2=sp2,rarefy.dist=50, rarefy.units="km", env.reso=0.416669,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)

## Adjust the number of variables input for e.vars after reduction to only important variables
num.var.e<-ncol(reduc.vars$env1)

## run it first with full environmental for background tests and equivalence statistic (total equivalence or divergence in current distributions)
full<-humboldt.doitall(inname="full_extent", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=0, reductype="PCA", non.analogous.environments="YES", correct.env=T, env.trim=T,  env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, pcx=1, pcy=2, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="YES",thresh.espace.z=0.0001, p.overlap=T, p.boxplot=F, p.scatter=F, run.silent=F, ncores=2)

## run it a second time with a trimmed, shared-espace. Here the equivalence statistic tests for niche evolution or niche divergence. For comparing results, change only the following model parameters: reduce.env, non.analogous.environmental, env.trim, nae
shared_ae<-humboldt.doitall(inname="shared_espace_ae", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=2, reductype="PCA", non.analogous.environments="NO", correct.env=T, env.trim=T, env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, pcx=1,pcy=2, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="NO", thresh.espace.z=0.0001, p.overlap=T, p.boxplot=F, p.scatter=T,run.silent=F, ncores=2)
```
### Example 3 - extended workflow outside of the 'do it all' function
The doitall function is a large wrapper function for most the analyses in humbodlt.  As a coder, such functions are typically avoided, because they can be challenging to troubleshoot.  If you are having any issues with that function, consider doing your analyses outside of it, as shown below.  See below for help formating raster/environment data. 
```markdown
library(humboldt)

## load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
env1<-read.delim("env1.txt",h=T,sep="\t")

## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
env2<-read.delim("env2.txt",h=T,sep="\t") 

## remove NAs and make sure all variables are imported as numbers
env1<-humboldt.scrub.env(env1)
env2<-humboldt.scrub.env(env2)

## load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))

## load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))

## perform ecological niche modeling to determine important variables
reduc.vars<- humboldt.top.env.brt(env1=env1,env2=env2,sp1=occ.sp1,sp2=occ.sp2,rarefy.dist=40, rarefy.units="km", env.reso=0.0833338,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)

## adjust the number of variables input for e.vars after reduction to only important variables based on niche modeling
num.var.e<-ncol(reduc.vars$env1)

## convert geographic space to environmental space, this step determines the epsace and if you are running a niche similarity test (reduce.env=0, non.alalogous.environments = "YES") of niche divergence test (reduce.env=2, non.alalogous.environments = "NO")
zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim= T, e.var=num.var.e,  col.env = e.var, trim.buffer.sp1 = 200, trim.buffer.sp2 = 200, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, R = 100, run.silent = F)

## store e-space scores for species and environments environments combined output from humboldt.g2e
scores.env1<-zz$scores.env1[1:2]
scores.env2<-zz$scores.env2[1:2]
scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
scores.sp1<-zz$scores.sp1[1:2]
scores.sp2<-zz$scores.sp2[1:2]

##  estimate environmental space in 2D 
z.sp1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
z.sp2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
z.env1<- humboldt.grid.espace(scores.env12,scores.env1,scores.env1,kern.smooth=1,R=100)
z.env2<- humboldt.grid.espace(scores.env12,scores.env2,scores.env2,kern.smooth=1,R=100)

## plot species estimated fundamental niches and habitats in espace, increase 'kern.smooth' if niche is too patchy
humboldt.plot.niche(z.sp1,"Species 1","PC1","PC2")
humboldt.plot.niche(z.sp2,"Species 2","PC1","PC2")
humboldt.plot.niche(z.env1,"Env 1","PC1","PC2")
humboldt.plot.niche(z.env2,"Env 2","PC1","PC2")

## quantify niche similarity between species 1 and 2
niche.sim.uncorrect<-humboldt.niche.similarity(z.sp1,z.sp2,correct.env=F)
niche.sim.correct<-humboldt.niche.similarity(z.sp1,z.sp2,correct.env=T)

## view results, Sorensen's D
niche.sim.correct$D
niche.sim.uncorrect$D

## perform niche equivalence tests 
niche.equiv<- humboldt.equivalence.stat(z.sp1,z.sp2,rep=100,kern.smooth=1, ncores=2)

## plot results
humboldt.plot.histrogram(niche.equiv,"D","Equivalence") 

## measure e-space correction coeff
ee<- humboldt.espace.correction(Z.env1=z.env1,Z.env2=z.env2,Z.sp1=z.sp1,Z.sp2=z.sp2)

## plot differences between species' espaces
humboldt.plot.espace.diff(ee, correct.env=F, type="species")

## plot contour lines of environment 1 if env1 and env2 are not identical
if(ee$s.uncor.sum!=0){
contour(z.env1$x,(sort((z.env1$y))),z.env1$Z,add=T,levels=quantile(z.env1$Z[z.env1$Z>0],c(0.1,0.5,0.75)),drawlabels=F,lty=c(1,2,3), lwd=c(1,1,1), col="grey")}

## plot differences between environments' espaces
humboldt.plot.espace.diff(ee, correct.env=F, type="environments")

## plot contour lines of environmental 1 if env1 and env2 are not identical
if(ee$e.uncor.sum!=0){
contour(z.env1$x,(sort((z.env1$y))),z.env1$Z,add=T,levels=quantile(z.env1$Z[z.env1$Z>0],c(0.1,0.5,0.75)),drawlabels=F,lty=c(1,2,3), lwd=c(1,1,1), col="grey")}

## perform background tests
bg.sp1tosp2<-humboldt.background.stat(g2e=zz, rep = 100, sim.dir = 1, env.reso=0.41666669, kern.smooth = 1, correct.env = T, R = 100, run.silent.bak = F)
bg.sp2tosp1<-humboldt.background.stat(g2e=zz, rep = 100, sim.dir = 2, env.reso=0.41666669, kern.smooth = 1, correct.env = T, R = 100, run.silent.bak = F)

## plot background tests 
humboldt.plot.density(bg.sp1tosp2,"D","Background 1->2") 
humboldt.plot.density(bg.sp2tosp1,"D","Background 2->1") 

## estimate the Potential Niche Truncation Index
pnt1<- humboldt.pnt.index(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
pnt2<- humboldt.pnt.index(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)

## plot pca contributions
humboldt.plot.contrib(zz$pca.cal$co,zz$pca.cal$eig)

## additional plots
humboldt.plot.overlap(in.g2e=zz, pdfname="FullEspaceOverlapSp1and2.pdf") 
humboldt.plot.scatter(env1[,3:4], xlab="Bio1", ylab="Bio2",main="environment")
```

## Getting environmental data formatted for ‘humboldt’
The best practice is to start with climate data that is similar to or matches the desired analysis resolution (for most taxa we recommend 5 to 20 km pixel resolution because this usually represents the spatial grain of autocorrelation of environmental variables; if you use a smaller resolution, this likely will not impact results and only slow analyses down).

To format the climate/environment raster data so that it can be used in Humbolt, follow these steps:
1.	Importing raster data into R
2.	Clip extent down to area of analysis 
3.	Convert input rasters to 'humboldt' format

### Importing rasters into R
```markdown
library(terra)

##load rasters to sample, here Tiff rasters are in base working directory
BIO1 = rast("bio1.tif")
BIO2 = rast("bio2.tif")
BIO12 = rast("bio12.tif")
BIO15= rast("bio15.tif")
```
### Clipping the extent of rasters in R
```markdown
##create raster stack to speed process up 
raster_stack_full<-c(BIO1, BIO2, BIO12, BIO15)

##define extent for which you will clip data to.  
##I highly recommend determining these values from a GIS (ArcGIS or Google Earth).
##extent order input below is: xmin, xmax, ymin, yma
e.in <- terra::ext(-160, 10, 30, 60)

##perform crop function to clip to input extent
raster_stack <- terra::crop(raster_stack_full, e.in)
```

### Convert input rasters to 'humboldt' format for use as an environment file
```markdown
##convert one of the rasters to a point dataframe to sample.  Use any raster input.
env.points<-as.points(BIO1, values=TRUE, na.rm=TRUE, na.all=FALSE)

##rarefy points to the appropriate analysis resolution.  Note, it is best practice to start with climate data that has the desired resolution of the analysis.  This process can take a lot of time, particularly if starting with very high-resolution datasets.  If the climate layer used is exactly the resolution desired (we recommend 10-40km2 for most studies), this step can be skipped.  Other environmental layers to be sampled can remain high resolution.
env.sampling.res<-humboldt.occ.rarefy(env.points, colxy = 1:2, rarefy.dist = 40,  rarefy.units = "km", run.silent.rar = F)

##subset only the x and y data
env.sampling.res<- env.sampling.res[,1:2]

##Extract values to points from rasters
RAST_VAL<-data.frame(terra::extract(raster_stack, env.sampling.res))

##merge sampled data to input
Env1<-cbind(env.sampling.es,RAST_VAL)

##if necessary, repeat for  environment 2

```

[![Analytics](https://ga-beacon.appspot.com/UA-124588717-1/humboldt)](https://github.com/igrigorik/ga-beacon)
