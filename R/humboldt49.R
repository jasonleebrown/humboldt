##################################################################################################
##################################################################################################
#' Scrub environment data
#' @param imported table for analysis  
#' @return This tool remove NAs and converts all input data to numeric values for analysis in Humboldt. 
#' @export
#' @examples
#' library(humboldt)
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#' 
#' env1<-humboldt.scrub.env(env1)
humboldt.scrub.env <- function(in_data) {
    l <- list()
    r.in <- nrow(in_data)
    num_data <- data.frame(data.matrix(in_data))
    numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x))) < 0.7})
    ncc<-length(numeric_columns[numeric_columns==TRUE])
    nall<-ncol(in_data)
    if (ncc!=nall){
    final_data <- data.frame(num_data[, numeric_columns], char_data[, !numeric_columns])
    l <- na.exclude(final_data)
    }
    if (ncc==nall){
    final_data <- num_data
    l <- na.exclude(final_data)
    }
    r.out <- nrow(l)
    print(paste((r.in - r.out), " rows of data (of ", r.in, " input,", (round((r.in - r.out)/r.in, 
        4) * 100), "%) were removed due to the presence of missing data or text"))
    return(l)
}

##################################################################################################
##################################################################################################
#' Create a grid of environmental space
#' @param sp pca values in 2 dimensions for the occurrences of the species in the ordination. 
#' @param glob pca values in 2 dimensions for the whole study area
#' @param glob1 pca values in 2 dimensions for the range of species
#' @param R resolution of grid in environmental space (RxR)
#' @param kern.smooth scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodel espace occupancy can be somewhat common when a species occupies an extreem aspect of habitat or when espace is not broadly acessible in both dimensions of espace (PCs 1 & 2).
#' @return This tool uses the scores of principal components (or 2 environmental variables) and creates a grid of RxR pixels with occurrence densities for species input. \cr \cr 
#' Only two dimensions can be used. If performing a PCA, I strongly encourage you to run a species distribution model using lots of environmental data on the focal species and then load only the top contributing 
#' environmental variables to be included in the PCA. This way all included variables are known to be relevant in both species distributions. This can be done inside Humboldt (via humboldt.top.env) by importing many environmental variables and letting program select only those important. Alternatively, this can also be done using another method (MaxEnt) outside of R and then use only variables deemed important in the species' distributions. \cr \cr 
#' @seealso \code{humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall, humboldt.top.env} that use or depend on outputs of this function 
#' @import sp
#' @import raster
#' @importFrom adehabitatMA ascgen
#' @importFrom geosphere lengthLine
#' @importFrom adehabitatHR kernelUD
#' @import rgeos
#' @export
#' @examples
#' library(humboldt)
#' 
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#' 
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#' 
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run Create a grid of Environmental Space Function
#' z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)

humboldt.grid.espace <- function(glob, glob1, sp, R = 100, kern.smooth = 1) {
    if (kern.smooth == "auto") {
        Kern.SmZ <- "href"
        Kern.SmZZ <- 0.5
    }
    if (kern.smooth != "auto") {
        Kern.SmZ <- 0.03 * kern.smooth
        Kern.SmZZ <- 20 * Kern.SmZ
    }
    # glob: global background dataset for the whole study area, glob1: background for sp1 sp:
    # occurrence dataset R: resolution of the grid
    glob <- as.matrix(glob)
    glob1 <- as.matrix(glob1)
    sp <- as.matrix(sp)
    
    th.sp = 0
    th.env = 0
    l <- list()
    
    if (ncol(glob) > 2) 
        stop("cannot calculate overlap with more than two axes")
    if (ncol(glob) == 1) 
        stop("cannot calculate overlap on one axes")
    
    if (ncol(glob) == 2) {
        # if scores in two dimensions (e.g. PCA)
        xmin <- min(glob[, 1] - Kern.SmZZ)
        xmax <- max(glob[, 1] + Kern.SmZZ)
        ymin <- min(glob[, 2] - Kern.SmZZ)
        ymax <- max(glob[, 2] + Kern.SmZZ)  # data preparation\t
        glob1r <- data.frame(cbind((glob1[, 1] - xmin)/abs(xmax - xmin), (glob1[, 2] - ymin)/abs(ymax - 
            ymin)))
        spr <- data.frame(cbind((sp[, 1] - xmin)/abs(xmax - xmin), (sp[, 2] - ymin)/abs(ymax - 
            ymin)))
        mask <- ascgen(SpatialPoints(cbind((1:R)/R, (1:R)/R)), nrcol = R - 2, count = FALSE)
        sp.dens <- kernelUD(SpatialPoints(spr[, 1:2]), h = Kern.SmZ, grid = mask, kern = "bivnorm")
        sp.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, matrix(sp.dens$ud, 
            nrow = R))
        glob1.dens <- kernelUD(SpatialPoints(glob1r[, 1:2]), grid = mask, kern = "bivnorm")
        glob1.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, matrix(glob1.dens$ud, 
            nrow = R))
        x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), length.out = R)
        y <- seq(from = min(glob[, 2]), to = max(glob[, 2]), length.out = R)
        glob1r <- extract(glob1.dens, glob1)
        Z.th <- quantile(glob1r, th.env)
        glob1.dens[glob1.dens < Z.th] <- 0
        Z <- glob1.dens * nrow(glob1)/cellStats(glob1.dens, "sum")
        spr <- extract(sp.dens, sp)
        z.th <- quantile(spr, th.sp)
        sp.dens[Z == 0] <- 0
        sp.dens[sp.dens < z.th] <- 0
        z <- sp.dens * nrow(sp)/cellStats(sp.dens, "sum")
        z.uncor <- z/cellStats(z, "max")
        z.cor <- z/Z
        z.cor[is.na(z.cor)] <- 0
        z.cor <- z.cor/cellStats(z.cor, "max")
        Z2 <- matrix(Z, nrow = R, ncol = R, byrow = F)
        z.cor2 <- matrix(z.cor, nrow = R, ncol = R, byrow = F)
        z.uncor2 <- matrix(z.uncor, nrow = R, ncol = R, byrow = F)
        #store values
        l$x <- x
        l$y <- y
        l$z.uncor <- z.uncor2
        l$z.cor <- z.cor2
        l$Z <- Z2
        l$glob <- glob
        l$glob1 <- glob1
        l$sp <- sp
        l$z.uncor.raster <- z.uncor
        l$z.cor.raster <- z.cor
        l$Z.raster <- Z
        }
    return(l)  #as.matrix(Z@data@values,nrow=R,ncol=R,byrow=F)as.matrix(z.cor@data@values)as.matrix(z.uncor@data@values)
}


##################################################################################################
##################################################################################################
#' Measure potential niche truncation index 
#' @param sp pca values in 2 dimensions for the occurrences of the species in the ordination. 
#' @param glob pca values in 2 dimensions for the whole study area
#' @param glob1 pca values in 2 dimensions for the range of species
#' @param R resolution of grid in environmental space (RxR)
#' @param kern.smooth scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodel espace occupancy can be somewhat common when a species occupies an extreem aspect of habitat or when espace is not broadly acessible in both dimensions of espace (PCs 1 & 2).
#' @return This tool estimates the Potential Niche Truncation Index, stored as the "pnt.index".  The value describes the amount of observed E-space of the species that is truncated by available E-space.  It is a measurement of the overlap between the 5 percent Kernel Density isopleth of the species and the 10 percent Kernel Density isopleth of accessible environment E-space.   The value is the portion of the species’ isopleth that falls outside of the environmental isopleth.   This value physically measures how much of the perimeter of the species' E-space abuts or overlaps with the margins of the environment’s E-space.  If the PNT index is moderate or high, keep in mind that your realized niche likely does not adequately reflect the species' fundamental niche.
#' @seealso \code{humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall, humboldt.top.env} that use or depend on outputs of this function 
#' @import sp
#' @import raster
#' @importFrom adehabitatMA ascgen
#' @importFrom geosphere lengthLine
#' @importFrom adehabitatHR kernelUD
#' @import rgeos
#' @export
#' @examples
#' library(humboldt)
#' 
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#' 
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#' 
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##its highly recommened that you using the function "humboldt.top.env" to select only the important enviromnetal variables. This step can be skipped. If you downloaded tons of environmental data, you should use this step.  If you skip this step, input env1/env2 inplace of reduc.vars$env1/reduc.vars$env2 
#' reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=occ.sp1,sp2=occ.sp2,rarefy.dist=40, rarefy.units="km", env.reso=0.416669,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)
#' 
#' ##Adjust the number of variables input for e.vars after reduction to only important variables
#' num.var.e<-ncol(reduc.vars$env1)
#'
#' ##convert geographic space to espace for measuring pnt.index
#' zz<-humboldt.g2e(env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 0, reductype = "PCA", non.analogous.environments = "YES", env.trim.mcp = T, e.var=c(3:num.var.e),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## estimate the Potential Niche Truncation Index
#' pnt1<- humboldt.pnt.index(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' pnt2<- humboldt.pnt.index(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)

humboldt.pnt.index <- function(glob, glob1, sp, R = 100, kern.smooth = 1) {
    if (kern.smooth == "auto") {
        Kern.SmZ <- "href"
        Kern.SmZZ <- 0.5
    }
    if (kern.smooth != "auto") {
        Kern.SmZ <- 0.03 * kern.smooth
        Kern.SmZZ <- 20 * Kern.SmZ
    }
    # glob: global background dataset for the whole study area, glob1: background for sp1 sp:
    # occurrence dataset R: resolution of the grid
    glob <- as.matrix(glob)
    glob1 <- as.matrix(glob1)
    sp <- as.matrix(sp)
    
    th.sp = 0
    th.env = 0
    l <- list()
    
    if (ncol(glob) > 2) 
        stop("cannot calculate overlap with more than two axes")
    if (ncol(glob) == 1) 
        stop("cannot calculate overlap on one axes")
    
    if (ncol(glob) == 2) {
        # if scores in two dimensions (e.g. PCA)
        xmin <- min(glob[, 1] - Kern.SmZZ*10)
        xmax <- max(glob[, 1] + Kern.SmZZ*10)
        ymin <- min(glob[, 2] - Kern.SmZZ*10)
        ymax <- max(glob[, 2] + Kern.SmZZ*10)  # data preparation\t
        glob1r <- data.frame(cbind((glob1[, 1] - xmin)/abs(xmax - xmin), (glob1[, 2] - ymin)/abs(ymax - 
            ymin)))
        spr <- data.frame(cbind((sp[, 1] - xmin)/abs(xmax - xmin), (sp[, 2] - ymin)/abs(ymax - 
            ymin)))
        mask <- ascgen(SpatialPoints(cbind((1:R)/R, (1:R)/R)), nrcol = R - 2, count = FALSE)
        sp.dens <- kernelUD(SpatialPoints(spr[, 1:2]), h = Kern.SmZ, grid = mask, kern = "bivnorm")
        sp.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, matrix(sp.dens$ud, 
            nrow = R))
        glob1.dens <- kernelUD(SpatialPoints(glob1r[, 1:2]), grid = mask, kern = "bivnorm")
        glob1.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, matrix(glob1.dens$ud, 
            nrow = R))
        x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), length.out = R)
        y <- seq(from = min(glob[, 2]), to = max(glob[, 2]), length.out = R)
        glob1r <- extract(glob1.dens, glob1)
        Z.th <- quantile(glob1r, th.env)
        glob1.dens[glob1.dens < Z.th] <- 0
        Z <- glob1.dens * nrow(glob1)/cellStats(glob1.dens, "sum")
        spr <- extract(sp.dens, sp)
        z.th <- quantile(spr, th.sp)
        sp.dens[Z == 0] <- 0
        sp.dens[sp.dens < z.th] <- 0
        z <- sp.dens * nrow(sp)/cellStats(sp.dens, "sum")
        z.uncor <- z/cellStats(z, "max")
        z.cor <- z/Z
        z.cor[is.na(z.cor)] <- 0
        z.cor <- z.cor/cellStats(z.cor, "max")
        Z2 <- matrix(Z, nrow = R, ncol = R, byrow = F)
        z.cor2 <- matrix(z.cor, nrow = R, ncol = R, byrow = F)
        z.uncor2 <- matrix(z.uncor, nrow = R, ncol = R, byrow = F)
        ratio=0
        Benv=contourLines(x, (sort((y * -1))), Z2, nlevels=1, quantile(Z2[Z2 > 0], c(0.15)))
        for (i in 1:length(Benv)) {		
           Bx<-Benv[[i]]$x
           By<-Benv[[i]]$y
           xyB<-cbind(Bx,By)
           p = Polygon(xyB)
           ps1 = Polygons(list(p),1)
           if (i ==1) {sBenv = SpatialPolygons(list(ps1))}
           if (i !=1) {sBenv = SpatialPolygons(list(ps1))+sBenv}
           }
        Bsp=contourLines(x, (sort((y * -1))), z.uncor2, nlevels=1, levels = quantile(z.uncor2[z.uncor2 > 0], c(0.05)))
        for (i in 1:length(Bsp)) {
           Ax<-Bsp[[i]]$x
           Ay<-Bsp[[i]]$y
           xyA<-cbind(Ax,Ay)
           p = Polygon(xyA)
           ps1 = Polygons(list(p),1)
           if (i ==1) {sBsp = SpatialPolygons(list(ps1))}
           if (i !=1) {sBsp = SpatialPolygons(list(ps1))+sBsp}
           }
        sBsp=gUnaryUnion(sBsp)
        sBenv=gUnaryUnion(sBenv)
        sBf=sBsp+sBenv
        borders = gDifference(
           as(sBf,"SpatialLines"),
           as(gUnaryUnion(sBf),"SpatialLines"),
           byid=TRUE)
        bordersL<-lengthLine(borders)
        bordersF<-lengthLine(sBsp)
        nborders=length(borders)
        bordersD<-bordersL[1]-bordersL[2]
        if (nborders==3){
             ratio<-bordersL[1]/bordersF}
        if (nborders==2 & bordersD>=1){
             ratio<-1
             }
        if (nborders==2 & bordersD<1){
             ratio<-0
             }
        print("Potential Niche Truncation Index:")
        print(ratio)
        if (ratio>1){
             ratio<-1}
        if (ratio<0.05){
             print("Minimal potential niche truncation")}
        if (ratio>=0.05 & ratio<0.15){
             print("Some potential niche truncation")}
        if (ratio>=0.15 & ratio<0.3){
             print("Moderate potential niche truncation")}			 
        if (ratio>=0.3 ){
             print("High potential niche truncation")}
        print("****************")
        l$pnt.index <- ratio		
        }
    return(l)  #as.matrix(Z@data@values,nrow=R,ncol=R,byrow=F)as.matrix(z.cor@data@values)as.matrix(z.uncor@data@values)
}


##################################################################################################
##################################################################################################
##################################################################################################
#' Measure niche overlap
#' @param z1 a grid of the density of a species occurrence in environmental space output from humboldt.grid.espace
#' @param z2 a grid of the density of a second species occurrence in environmental space output from humboldt.grid.espace that you wish to compared to z1
#' @param correct.env if correct.env=T, the analysis corrects occurrence densities of each species by the prevalence of the environments in their range. If correct.env=F, the overlap measure does not correct occurrence densities of each species by the prevalence of the environments in their range
#' @param nae do you include non-analogous environments in the niche overlap measurement? If nae="NO" (use capital letters), then non-analogous environments will be removed from both input environments during overlap measurement and only environments present in both will be used to measure overlap. If nae="YES" then no change will be made to input z1 and z2. \cr \cr Note: this is separate from trimming non-analogous environments from your input dataset (as done by humbodlt.g2 specified by parameter non.analogous.environments). The latter parameter physically removes non-analogous environments from datasets pre-niche overlap measurement. Technically the removal of non-analogous environments via either way should result in similar overlap measurements (though they may not be identical). This because removing non-analogous environmental from the dataset prior to gridding environments will resulting only non-analogous environments to be gridded (and typically finer grain applied to each grid cell). Whereas removing them only via this parameter (nae), which only removes non-analogous in the gridded environmental space for use in overlap measurements,--- all the input environmental space is gridded (likely increasing the environmental space per gridded cell). \cr \cr A second cause of differences in values can result from rescaling of espace values during niche-overlap measurements so that the sum of the landscape equals one. If occupied, non-analogous environments are numerous in one of the datasets, this can theoretically cause overlap values to decrease in analogous environments (vs. non-analogous environments) because differences in gridded niches are rescaled to 1 in both scenarios. The rescaling among fewer cells increases the values applied to highly suitable areas and, if not equivalently scaled in both datasets, differences among niches could increase, resulting a smaller overlap in non-analogous environments (again values should be similar). If you remove non-analogous environments in humboldt.g2e, I also suggest that you use this function (as it can remove any slight anomalies caused by gridding environments in humbodlt.grid.clim due to the binning of values in the RxR grid).
#' @param thresh.espace.z this parameter is an experimental parameter and controls the level at which values below the kernel density z values are removed for creating areas of analogous environmental space. Higher values will increase value from which the low-density areas are removed from the environmental space of z1 and z2. Basically values above this are retained and values below are removed. Default=0.001
#' @return Compares two z1 and z2 created by humboldt.grid.espace by calculating the overlap metrics: Schoener's D (herein called 'D') and a measure derived from Hellinger called 'I'(herein called 'I'; see Warren et al. 2008 & Broennimann et al. 2012). Both measures are based on two species occurrence density grids z1 and z2 created by humboldt.grid.espace\cr
#'\cr Output are six items: $D= Schoener's D value, $I= Hellinger's I value, $remaining.espace.per= the remaining espace remaining after removal of non-analogous espace, $num.cells.analog.espace= the number of cells of analogous espace, $num.cells.total.espace= the total number of cells of espace present(both analogous and non-analogus espace)
#' @seealso \code{humboldt.sample.spp,humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.plot.niche, humboldt.doitall} which use or depend on outputs of this function 
#' @import sp
#' @import rgdal
#' @import raster
#' @export
#' @examples
#' #######################################################################################
#' ###################################    EXAMPLE 1    ###################################
#' #######################################################################################
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run create a grid of Environmental Space Function
#' z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
#' 
#' ## mesure niche overlap
#' humboldt.niche.overlap(z1,z2,correct.env=F)
#' 
#' #######################################################################################
#' ###################################    EXAMPLE 2    ###################################
#' #######################################################################################
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#' 
#' ##merge environment files
#' env12<-rbind(env1,env2)
#'
#' ##load occurrence sites for the species at study area 1. Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2. Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#'
#' ##remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
#' occ.sp1<-humboldt.occ.rarefy(in.pts=occ.sp1,colxy=2:3, rarefy.dist=40,rarefy.units="km")
#' occ.sp2<-humboldt.occ.rarefy(in.pts=occ.sp2,colxy=2:3, rarefy.dist=40,rarefy.units="km")
#'
#' ##sample environment using humboldt.sample.spp() function
#' ## env.reso should be the resolution of the environmental data grid
#' occ.sp1<-na.exclude(humboldt.sample.spp(dfsp=occ.sp1,colspxy=2:3,colspkept=NULL,dfvar=env1,colvarxy=1:2,colvar="all",resolution=0.16667))
#' occ.sp2<-na.exclude(humboldt.sample.spp(dfsp=occ.sp2,colspxy=2:3,colspkept=NULL,dfvar=env2,colvarxy=1:2,colvar="all",resolution=0.16667))
#'
#' ##row weighting of environment density for PCA
#' row.w.1.env<-1-(nrow(env1)/nrow(env12))  # prevalence of env1
#' row.w.2.env<-1-(nrow(env2)/nrow(env12))  # prevalence of env2
#' row.w.env <-c(rep(row.w.1.env, nrow(env1)),rep(row.w.2.env, nrow(env2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
#'
#' ##create a global dataset with all environments and species, include environmental variables e.var[x1-xN]
#' e.var<-c(3:21)
#' data.env.occ<-rbind(env1,env2,occ.sp1,occ.sp2)[e.var]
#'
#' ##perform PCA of environment and species 1 data which is weighted by density of environment types in both environment 1 and environment 2
#' pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
#'
#' ##store PCA scores for sp1 andenvironments 1,2 and both environments combined
#' ## specific exact locations of environment and sp. data in PCA for storage
#' row.env1<-1:nrow(env1)
#' row.env2<-(nrow(env1)+1):(nrow(env1)+nrow(env2))
#' row.env12<-1:(nrow(env1)+nrow(env2))
#' row.sp1<-(nrow(env1)+nrow(env2)+1):(nrow(env1)+nrow(env2)+nrow(occ.sp1))
#' row.sp2<-(nrow(env1)+nrow(env2)+nrow(occ.sp1)+1):(nrow(env1)+nrow(env2)+nrow(occ.sp1)+nrow(occ.sp2))
#' ##glob: global background dataset for the whole study area 
#' glob<- pca.cal$li[row.env12,]
#' ##glob1: background for sp1
#' glob1<- pca.cal$li[row.env1,]
#' ##glob2: background for sp2
#' glob2<- pca.cal$li[row.env2,]
#' ##sp1: occurrence dataset
#' sp1<- pca.cal$li[row.sp1,]
#' ##sp2: occurrence dataset
#' sp2<- pca.cal$li[row.sp2,]
#'
#' ## run Create a grid of Environmental Space Function
#' z1<- humboldt.grid.espace(glob,glob1,sp1,R=100,kern.smooth=1)
#' z2<- humboldt.grid.espace(glob,glob2,sp2,R=100,kern.smooth=1)
#' humboldt.niche.overlap(z1,z2,correct.env=F)

humboldt.niche.overlap <- function(z1, z2, correct.env = F, nae = "NO", thresh.espace.z = 0.001) {
    R <- length(z1$x)
    l <- list()
    
    z1a <- z1
    z2a <- z2
    z2a$Z[z2a$Z < thresh.espace.z] <- 0
    z1a$Z[z1a$Z < thresh.espace.z] <- 0
    z2a$Z[z2a$Z >= thresh.espace.z] <- 1
    z1a$Z[z1a$Z >= thresh.espace.z] <- 1
    za.sum <- z1a$Z + z2a$Z
    za.sum[za.sum > 0] <- 1
    za.sum.ci <- sum(za.sum)
    za.sum.mask <- z1a$Z + z2a$Z
    za.sum.mask[za.sum.mask <= 1] <- 0
    za.sum.mask[za.sum.mask > 1] <- 1
    za.sum.cf <- sum(za.sum.mask)
    
    if (nae == "NO" & correct.env == F) {
        z1$z.uncor[za.sum.mask == 0] <- 0
        z2$z.uncor[za.sum.mask == 0] <- 0
        Clc <- round(100 * (za.sum.cf/za.sum.ci), 2)
    }
    if (nae == "NO" & correct.env == T) {
        z1$z.cor[za.sum.mask == 0] <- 0
        z2$z.cor[za.sum.mask == 0] <- 0
        Clc <- round(100 * (za.sum.cf/za.sum.ci), 2)
    }
    if (nae == "YES") {
        Clc <- 100
    }
    
    if (correct.env == F) {
        z1$z.uncor[z1$z.uncor < thresh.espace.z] <- 0
        z2$z.uncor[z2$z.uncor < thresh.espace.z] <- 0
        p1 <- z1$z.uncor/sum(z1$z.uncor)  # rescale occurrence densities so that the sum of densities is the same for both species
        p2 <- z2$z.uncor/sum(z2$z.uncor)  # rescale occurrence densities so that the sum of densities is the same for both species
    }
    
    if (correct.env == T) {
        z1$z.cor[z1$z.cor < thresh.espace.z] <- 0
        z2$z.cor[z2$z.cor < thresh.espace.z] <- 0
        p1 <- z1$z.cor/sum(z1$z.cor)  # rescale occurrence densities so that the sum of densities is the same for both species
        p2 <- z2$z.cor/sum(z2$z.cor)  # rescale occurrence densities so that the sum of densities is the same for both species
    }
    
    D <- 1 - (0.5 * (sum(abs(p1 - p2))))  # overlap metric D
    Diff.I <- sqrt(sum((sqrt(p1) - sqrt(p2))^2))
    I <- 1 - (Diff.I^2)/2  # overlap metric I
    l$D <- D
    l$I <- I
    l$remaining.espace.per <- Clc
    l$num.cells.analog.espace <- za.sum.cf
    l$num.cells.total.espace <- za.sum.ci
    l$mask.analog.e <- za.sum.mask
    l$mask.total.e <- za.sum
    mask.analog.e.raster <- raster(matrix(za.sum.mask, nrow = R))
    crs(mask.analog.e.raster) <- crs(z1$Z.raster)
    extent(mask.analog.e.raster) <- extent(z1$Z.raster)
    l$mask.analog.e.raster <- mask.analog.e.raster
    mask.total.e.raster <- raster(matrix(za.sum, nrow = R))
    crs(mask.total.e.raster) <- crs(z1$Z.raster)
    extent(mask.total.e.raster) <- extent(z1$Z.raster)
    l$mask.total.e.raster <- mask.total.e.raster
    z1.raster <- raster(matrix(p1, nrow = R))
    crs(z1.raster) <- crs(z1$Z.raster)
    extent(z1.raster) <- extent(z1$Z.raster)
    l$z1.raster <- z1.raster
    z2.raster <- raster(matrix(p2, nrow = R))
    crs(z2.raster) <- crs(z1$Z.raster)
    extent(z2.raster) <- extent(z1$Z.raster)
    l$z2.raster <- z2.raster
    
    return(l)
}

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#' Niche equivalency test
#' @param z1 a grid of the density of a species occurrence in environmental space output from humboldt.grid.espace
#' @param z2 a grid of the density of a second species occurrence in environmental space output from humboldt.grid.espace that you wish to compared to z1
#' @param rep is the number of iterations. Values higher than 100 is recommend for final analysis.
#' @param sample parameter controls how espace is sampled.  There are two options: "points" or "density".  If sample="points", only input localities for both species are used for randomizations. If sample=density, the gridded density of species' espace is randomly resampled with higher sampling probabilities applied to areas of high density espace.  In both case the number of samples select for each species is equal to the observed datasets.  Again 'points' focuses only on input points, whereas 'density' samples the entire gridded density of each species. Generally I prefer using 'points' for datasets, because this maintains most of the nuances associated the original datasets (i.e. spatial autocorrelation).  Though sometimes I do use 'density' for datasets with limited localities (<20 points)
#' @param nae do you include non-analogous environments in the niche overlap measurement? If nae="NO" (use capital letters), then non-analogous environments will be removed from both input environments during overlap measurement and only environments present in both will be used to measure overlap. If nae="YES" then no change will be made to input z1 and z2. Note: this is separate from trimming non-analogous environments from your input dataset (as done by humbodlt.g2 specified by parameter non.analogous.environments). This parameter physically removes non-analogous environments from datasets pre-niche overlap measurement. Technically the removal of non-analogous environments via either way should result in similar overlap measurements (though they may not be identical). This because removing NAE from the dataset prior to gridding environments will resulting only non-analogous environments to be gridded (and typically finer grain applied to each grid cell). Whereas removing them only via this parameter (nae), which only removes non-analogous in the gridded environmental space for use in overlap measurements--- all the input environmental space is gridded (likely increasing the environmental space per gridded cell). A second cause of differences in values can result from rescaling of espace values during niche-overlap measurements so that the sum of the landscape equals one. If occupied non-analogous environmental are numerous in one of the datasets, this can theoretically cause overlap values to decrease in analogous environments (vs. nae) because differences in core niches are rescaled to 1 in both scenarios. The rescaling among fewer cells increases the values applied to highly suitable areas and, if not equivalently scaled in both datasets, differences among niches could increase, resulting a smaller overlap in non-analogous environments (again values should be similar). If you remove non-analogous environments in humboldt.g2e, I also suggest that you use this function (as it can remove any slight anomalies caused by gridding environments in humbodlt.grid.clim due to the binning of values in the RxR grid).
#' @param thresh.espace.z this parameter is an experimental parameter and controls the level at which values below the kernel density z values are removed for creating areas of analogous environmental space. Higher values will increase value from which the low-density areas are removed from the environmental space of z1 and z2.  Basically values above this are retained and values below are removed. Default=0.001
#' @param kern.smooth scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodal espace occupancy can be somewhat common when a species occupies an extreme aspect of habitat or when espace is not broadly accessible in both dimensions of espace (PCs 1 & 2)
#' @return runs a modified niche equivalency test(see Warren et al 2008, Broenniman et al. 2012, but note these analyses are not identical) based on two species occurrence density grids and compares the observed niche overlap between z1 and z2 (created by humboldt.grid.espace) to overlap between random niches z1.sim and z2.sim. The z1.sim and z2.sim are built from random reallocations of occurrences of z1 and z2 (see 'sample' parameter). A significant value states that the two datasets are NOT statistically equivalent and rejects the NULL hypothesis that species niches are equivalent. \cr 
#' \cr
#' Output: $sim= simulation niche overlap values, $obs.D= Schoener's D value for two observed datasets, $obs.I= Hellinger's I value for observed datasets, $p.D=one-tail p-value of Schoener's D values (simulation vs. observed), $p.I=one-tail p-value of Hellinger's I value (simulation vs. observed) 
#' @seealso \code{humboldt.niche.overlap, humboldt.doitall, humboldt.plot.density, humboldt.plot.histrogram} which use or depend on outputs of this function 
#' @import sp
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run create a grid of Environmental Space Function
#' z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
#' 
#' ## mesure niche equivalency
#' niche.equiv<- humboldt.equivalency.test(z1,z2,rep=10,kern.smooth=1)

humboldt.equivalency.test <- function(z1, z2, rep = 100, correct.env = T, kern.smooth = 1, nae = "YES", 
    thresh.espace.z = 0.001, sample = "points", run.silent.equ = F) {
    threshinZ <- thresh.espace.z
    kern.smoothinZ <- kern.smooth
    nacinZ <- nae
    Rin <- length(z1$x)
    l <- list()
    if (run.silent.equ == F) {
        x11(2, 2, pointsize = 12)
        par(mar = c(0, 0, 0, 0))
    }
    if (sample == "density") {
        if (correct.env == T) {
            obs.o <- humboldt.niche.overlap(z1, z2, correct.env = T, nae = nacinZ, thresh.espace.z = threshinZ)  #observed niche overlap
            sim.o <- data.frame(matrix(nrow = rep, ncol = 2))  #empty list of random niche overlap
            names(sim.o) <- c("D", "I")
            if (!is.null(z1$y)) {
                # overlap on two axes
                coordinates <- which(z1$z.cor > 0, arr.ind = T)  # array of cell coordinates ((1,1),(1,2)...)
                weight <- z1$z.cor[z1$z.cor > 0]  #densities in the same format as cells
                coordinates.sim1 <- coordinates[sample(1:nrow(coordinates), size = nrow(z1$sp), 
                  replace = F, prob = weight), ]  #random sampling of coordinates following z1$z.cor distribution
                occ1.sim <- cbind(z1$x[coordinates.sim1[, 1]], z1$y[coordinates.sim1[, 2]])  # random occurrences following the corrected densities distribution
                
                coordinates <- which(z2$z.cor > 0, arr.ind = T)  # array of cell coordinates ((1,1),(1,2)...)
                weight <- z2$z.cor[z2$z.cor > 0]  #densities in the same format as cells
                coordinates.sim2 <- coordinates[sample(1:nrow(coordinates), size = nrow(z2$sp), 
                  replace = F, prob = weight), ]  #random sampling of coordinates following z1$z.cor distribution
                occ2.sim <- cbind(z2$x[coordinates.sim2[, 1]], z2$y[coordinates.sim2[, 2]])  # random occurrences following the corrected densities distribution
                
                o.i <- humboldt.niche.overlap(z1.sim, z2.sim, correct.env = T, nae = nacinZ, 
                  thresh.espace.z = threshinZ)  # overlap between random and observed niches
                # print(o.i)
                sim.o$D[i] <- o.i$D  # storage of overlaps
                sim.o$I[i] <- o.i$I
            }
            
            if (run.silent.equ == F) {
                dev.off()
            }
        }
        if (correct.env == F) {
            obs.o <- humboldt.niche.overlap(z1, z2, correct.env = F, nae = nacinZ, thresh.espace.z = threshinZ)  #observed niche overlap
            sim.o <- data.frame(matrix(nrow = rep, ncol = 2))  #empty list of random niche overlap
            names(sim.o) <- c("D", "I")
            for (i in 1:rep) {
                if (run.silent.equ == F) {
                  plot.new()
                  text(0.5, 0.5, paste("equivalency tests:", "\n", "runs to go:", rep - i + 1))
                }
                
                if (is.null(z1$y)) {
                  # overlap on one axis
                  
                  occ1.sim <- sample(z1$x, size = nrow(z1$sp), replace = F, prob = z1$z.cor)  #random sampling of occurrences following the corrected densities distribution
                  occ2.sim <- sample(z2$x, size = nrow(z2$sp), replace = F, prob = z2$z.cor)
                  
                  occ.pool <- c(occ1.sim, occ2.sim)  # pool of random occurrences
                  rand.row <- sample(1:length(occ.pool), length(occ1.sim), replace = F)  # random reallocation of occurrences to datasets
                  sp1.sim <- occ.pool[rand.row]
                  sp2.sim <- occ.pool[-rand.row]
                  
                  z1.sim <- humboldt.grid.espace(z1$glob, z1$glob1, data.frame(sp1.sim), kern.smooth = kern.smoothinZ, 
                    R = Rin)  # gridding
                  z2.sim <- humboldt.grid.espace(z2$glob, z2$glob1, data.frame(sp2.sim), kern.smooth = kern.smoothinZ, 
                    R = Rin)
                }
                
                if (!is.null(z1$y)) {
                  # overlap on two axes
                  coordinates <- which(z1$z.cor > 0, arr.ind = T)  # array of cell coordinates ((1,1),(1,2)...)
                  weight <- z1$z.cor[z1$z.cor > 0]  #densities in the same format as cells
                  coordinates.sim1 <- coordinates[sample(1:nrow(coordinates), size = nrow(z1$sp), 
                    replace = F, prob = weight), ]  #random sampling of coordinates following z1$z.cor distribution
                  occ1.sim <- cbind(z1$x[coordinates.sim1[, 1]], z1$y[coordinates.sim1[, 2]])  # random occurrences following the corrected densities distribution
                  
                  coordinates <- which(z2$z.cor > 0, arr.ind = T)  # array of cell coordinates ((1,1),(1,2)...)
                  weight <- z2$z.cor[z2$z.cor > 0]  #densities in the same format as cells
                  coordinates.sim2 <- coordinates[sample(1:nrow(coordinates), size = nrow(z2$sp), 
                    replace = F, prob = weight), ]  #random sampling of coordinates following z1$z.cor distribution
                  occ2.sim <- cbind(z2$x[coordinates.sim2[, 1]], z2$y[coordinates.sim2[, 2]])  # random occurrences following the corrected densities distribution
                  occ.pool <- rbind(occ1.sim, occ2.sim)  # pool of random occurrences
                  rand.row <- sample(1:nrow(occ.pool), nrow(occ1.sim), replace = F)  # random reallocation of occurrences to datasets
                  sp1.sim <- occ.pool[rand.row, ]
                  sp2.sim <- occ.pool[-rand.row, ]
                  z1.sim <- humboldt.grid.espace(z1$glob, z1$glob1, data.frame(sp1.sim), kern.smooth = kern.smoothinZ, 
                    R = Rin)
                  z2.sim <- humboldt.grid.espace(z2$glob, z2$glob1, data.frame(sp2.sim), kern.smooth = kern.smoothinZ, 
                    R = Rin)
                }
                
                o.i <- humboldt.niche.overlap(z1.sim, z2.sim, correct.env = F, nae = nacinZ, 
                  thresh.espace.z = threshinZ)  #T\t\t\t\t\t\t# overlap between random and observed niches
                sim.o$D[i] <- o.i$D  # storage of overlaps
                sim.o$I[i] <- o.i$I
            }
            is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
            
            sim.o[is.nan(sim.o)] <- 0
            
            if (run.silent.equ == F) {
                dev.off()
            }
        }
    }
    
    if (sample == "points") {
        if (correct.env == T) {
            obs.o <- humboldt.niche.overlap(z1, z2, correct.env = T, nae = nacinZ, thresh.espace.z = threshinZ)  #observed niche overlap
            sim.o <- data.frame(matrix(nrow = rep, ncol = 2))  #empty list of random niche overlap
            names(sim.o) <- c("D", "I")
            for (i in 1:rep) {
                if (run.silent.equ == F) {
                  plot.new()
                  text(0.5, 0.5, paste("equivalency tests:", "\n", "runs to go:", rep - i + 1))
                }
                if (!is.null(z1$y)) {
                  # overlap on two axes
                  all.points <- rbind(z1$sp, z2$sp)
                  row.names(all.points) <- c()
                  coordinates.sim <- all.points[sample(1:nrow(all.points), size = nrow(all.points), 
                    replace = F), ]  #random sampling of coordinates following z1$z.cor distribut
                  coordinates.sim1 <- coordinates.sim[1:nrow(z1$sp), ]
                  coordinates.sim2 <- coordinates.sim[(nrow(z1$sp) + 1):(nrow(z1$sp) + nrow(z2$sp)), 
                    ]
                  z1.sim <- humboldt.grid.espace(z1$glob, z1$glob1, data.frame(coordinates.sim1), 
                    kern.smooth = kern.smoothinZ, R = Rin)
                  z2.sim <- humboldt.grid.espace(z2$glob, z2$glob1, data.frame(coordinates.sim2), 
                    kern.smooth = kern.smoothinZ, R = Rin)
                }
                o.i <- humboldt.niche.overlap(z1.sim, z2.sim, correct.env = T, nae = nacinZ, 
                  thresh.espace.z = threshinZ)  # overlap between random and observed niches
                # print(o.i)
                sim.o$D[i] <- o.i$D  # storage of overlaps
                sim.o$I[i] <- o.i$I
            }
            
            if (run.silent.equ == F) {
                dev.off()
            }
        }
        
        if (correct.env == F) {
            obs.o <- humboldt.niche.overlap(z1, z2, correct.env = F, nae = nacinZ, thresh.espace.z = threshinZ)  #observed niche overlap
            sim.o <- data.frame(matrix(nrow = rep, ncol = 2))  #empty list of random niche overlap
            names(sim.o) <- c("D", "I")
            for (i in 1:rep) {
                if (run.silent.equ == F) {
                  plot.new()
                  text(0.5, 0.5, paste("equivalency tests:", "\n", "runs to go:", rep - i + 1))
                }
                if (!is.null(z1$y)) {
                  # overlap on two axes
                  all.points <- rbind(z1$sp, z2$sp)
                  row.names(all.points) <- c()
                  coordinates.sim <- all.points[sample(1:nrow(all.points), size = nrow(all.points), 
                    replace = F), ]  #random sampling of coordinates following z1$z.cor distribution
                  coordinates.sim1 <- coordinates.sim[1:nrow(z1$sp), ]
                  coordinates.sim2 <- coordinates.sim[(nrow(z1$sp) + 1):(nrow(z1$sp) + nrow(z2$sp)), 
                    ]
                  z1.sim <- humboldt.grid.espace(z1$glob, z1$glob1, data.frame(coordinates.sim1), 
                    kern.smooth = kern.smoothinZ, R = Rin)
                  z2.sim <- humboldt.grid.espace(z2$glob, z2$glob1, data.frame(coordinates.sim2), 
                    kern.smooth = kern.smoothinZ, R = Rin)
                }
                
                o.i <- humboldt.niche.overlap(z1.sim, z2.sim, correct.env = F, nae = nacinZ, 
                  thresh.espace.z = threshinZ)  #T\t\t\t\t\t\t# overlap between random and observed niches
                sim.o$D[i] <- o.i$D  # storage of overlaps
                sim.o$I[i] <- o.i$I
            }
            # is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
            
            # sim.o[is.nan(sim.o)] <- 0
            
            if (run.silent.equ == F) {
                dev.off()
            }
        }
    }
    l$sim <- sim.o  # storage
    l$obs$D <- obs.o$D  # storage
    l$obs$I <- obs.o$I  # storage
    # l$p.D<- min((sum(sim.o$D <= obs.o$D ) + 1),(sum(sim.o$D >= obs.o$D ) +
    # 1))*1/(length(sim.o$D) + 1) # storage of p-values l$p.I<- min((sum(sim.o$I <= obs.o$I ) +
    # 1),(sum(sim.o$I >= obs.o$I ) + 1))*1/(length(sim.o$I) + 1) #
    
    p.D <- (sum(sim.o$D <= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values
    l$p.D <- p.D
    l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of 
    return(l)
}

##################################################################################################
#' Niche background test
#' @param rep is the number of iterations. Values higher than 200 are recommend for final analysis.
#' @param run.silent.bak if run.silent=T, texts boxes displaying progress will not be displayed
#' @param g2e an espace file output from humboldt.g2e

#' @param sim.dir if sim.dir=1, z2 will be randomly shifted the resulting niche will be compared to z1 (which is unchanged). If sim.dir=2, z1 will be randomly shifted the resulting niche will be compared to z2 (which is unchanged)
#' @param env.reso the resolution of the environmental data grid (typically in decimal degrees)
#' @param kern.smooth scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodal espace occupancy can be somewhat common when a species occupies an extreme aspect of habitat or when espace is not broadly accessible in both dimensions of espace (PCs 1 & 2)
#' @param correct.env if correct.env=T, the analysis corrects occurrence densities of each species by the prevalence of the environments in their range. If correct.env=F, the overlap measure does not correct occurrence densities of each species by the prevalence of the environments in their range.  Default value is FALSE
#' @param R resolution of grid in environmental space (RxR). This needs to match the value input into humboldt.g2e.  Default value is 100
#' @return Runs a modified niche background test(see Warren et al 2008) based on two species' occurrence density grids. The function compares the observed niche overlap between z1 and z2 (created by humboldt.grid.espace) to overlap between z1 and the random shifting of the spatial distribution of z2 in geographic space and then measuring how that shift in geography changes occupied environmental space (called z2.sim). This test maintains the spatial structure of all the localities and thereby retains all nuances associated with each datasets' spatial autocorrelation (vs. random sampling). This test asks if the two distributed populations/species are more different than would be expected given the underlying environmental differences between the regions in which they occur. \cr
#' \cr
#' If the observed values of the niche similarity measures obtained from the two original populations are significantly higher than expected from this null distribution, then the null hypothesis that similarity between species is solely due to differences in habitat availability is rejected. Or in other words, a significant value suggests that the two species are more divergent than expected solely on the habitat availability (pending a significant equivalency test). A non-significant test suggest that most occupied environments are similar among environments. \cr
#' \cr
#' IMPORTANT. This test measures the power of the equivalency test to detect differences based on the available e-space. If both the equivalency test and background test are non-significant, this means the two species occupied environmental spaces are not significantly different and resulting niche 'equivalence' is likely a result of the limited environmental space present in habitat(s).  Basically in these situations, there is limited power for the equivalency test to actually detect and significant differences among taxa, even if they existed. However, conversely it also doesn't provide any evidence that they in fact on not equivalent- simply there is little power to detect it in input environmental data. \cr 
#' \cr
#' If both background test are non-significant and equivalency test are non-significant, try to increase the spatial extent of input climate data.  If the equivalency test is significant and this is not, this means the two species occupied environmental space is significantly different despite the existence of largely similar habitats. This is strong evidence of niche divergence.\cr
#' \cr
#' Output: $sim= simulation values for D, I, nDg (number times D.sim is greater than D.obs), nDl (number times D.sim is less than D.obs), nIg (number times I.sim is greater than I.obs), nIl (number times I.sim is less than I.obs), n.pts=number of points in random sample (often when shifting the spatial distribution of z2 in geographic space, points are moved to areas with no environmental data and thereby are excluded for simulation niche calculation), x.shift= shift on longitude, y.shift= shift in latitude; $obs= D & I values in observed datasets; $p.D=one-tailed p-value of Schoener's D values (simulation vs. observed), $p.I=one-tailed p-value of Hellinger's I value (simulation vs. observed).\cr
#' @seealso \code{humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz=humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 0, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##perform background tests 
#' bg.sp1tosp2<-humboldt.background.test(g2e=zz, rep = 100, sim.dir = 1, env.reso=0.41666669, kern.smooth = 1, correct.env = F, R = 100, run.silent.bak = F)
#' bg.sp2tosp1<-humboldt.background.test(g2e=zz, rep = 100, sim.dir = 2, env.reso=0.41666669, kern.smooth = 1, correct.env = F, R = 100, run.silent.bak = F)
humboldt.background.test<- function(g2e, rep = 100, sim.dir = 1, env.reso, kern.smooth = 1, 
    correct.env = F, R = 100, force.equal.sample=F, run.silent.bak = F) {
    kern.smoothinZ <- kern.smooth
    Rin = R
    env.reso.sim <- env.reso * 2
    inE <- correct.env
    scores.env12a <- rbind(g2e$scores.env1, g2e$scores.env2)
    if (sim.dir == 1) {
        scores.env1a <- g2e$scores.env1
        scores.env2a <- g2e$scores.env2
        scores.sp1a <- g2e$scores.sp1
        scores.sp2a <- g2e$scores.sp2  #[1:2]
        scores.env2a <- g2e$scores.env2[1:2]
        scores.env2b <- g2e$scores.env2[3:4]
    }
    if (sim.dir == 2) {
        scores.env1a <- g2e$scores.env2
        scores.env2a <- g2e$scores.env1
        scores.sp1a <- g2e$scores.sp2
        scores.sp2a <- g2e$scores.sp1  #[1:2]
        scores.env2a <- g2e$scores.env1[1:2]
        scores.env2b <- g2e$scores.env1[3:4]
    }
    scores.env2c <- cbind(scores.env2b, scores.env2a)
    # calculation of occurrence density and test of niche equivalency and Background
    z1 <- humboldt.grid.espace(scores.env12a[1:2], scores.env1a[1:2], scores.sp1a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    z2 <- humboldt.grid.espace(scores.env12a[1:2], scores.env2a[1:2], scores.sp2a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    ####### R<-length(z1$x)
    if (run.silent.bak == F) 
        {
            x11(2, 2, pointsize = 12)
            par(mar = c(0, 0, 0, 0))
        }  # countdown window
    l <- list()
    sim.o <- data.frame(matrix(nrow = rep, ncol = 9))  #empty list of random niche overlap
    ## jlb added
    obs.o <- data.frame(matrix(nrow = rep, ncol = 2))  #empty list of 
    ## trim observed niche overlap to shared environmental space
    names(sim.o) <- c("D", "I", "nDg", "nDl", "nIg", "nIl", "n.pts", "x.shift", "y.shift")
    names(obs.o) <- c("D", "I")
    
    ## figure out sample distance
    AmaxFG <- max(scores.env2b[, 1])
    AminFG <- min(scores.env2b[, 1])
    BmaxFG <- max(scores.env2b[, 2])
    BminFG <- min(scores.env2b[, 2])
    # do math for background test sample distance
    Env2.d.x <- abs(AmaxFG - AminFG)/3
    Env2.d.y <- abs(BmaxFG - BminFG)/3
    Env2.d <- min(Env2.d.x, Env2.d.y)
    
    ### need to sample at clim Res Create sample windows
    leg.int <- as.integer(Env2.d/env.reso)
    Rana <- seq(-1 * Env2.d, Env2.d, length.out = (leg.int) + 1)
    o.i <- humboldt.niche.overlap(z1, z2, correct.env = inE, nae = "YES")
    if (force.equal.sample==T){
      for (k in 1:rep) {
        if (run.silent.bak == F) 
            {
                plot.new()
                text(0.5, 0.5, paste("Background tests:", "\n", "runs to go:", rep - k + 1))
            }  # countdown
        {
            Xran <- sample(Rana, 1)
            Yran <- sample(Rana, 1)
            sp2rX <- scores.sp2a[, 3] + Xran
            sp2rY <- scores.sp2a[, 4] + Yran
            sp2r <- cbind(sp2rX, sp2rY)
            sp2r <- as.data.frame.matrix(sp2r)
            nsp2t<-nrow(sp2r)
            sp2rPCAs1 <<- na.exclude(humboldt.sample.spp(dfsp = sp2r, colspxy = 1:2, colspkept = NULL, 
                dfvar = scores.env2c, colvarxy = 1:2, colvar = 3:4, resolution = env.reso.sim, 
                run.silent = T))
            nsp2r <- nrow(sp2rPCAs1)
            nsp2.n.miss<-nsp2t-nsp2r
            nsp2.mis<<-scores.env2c[sample(nrow(scores.env2c), nsp2.n.miss), ][,3:4]
            sp2rPCA<-rbind(sp2rPCAs1,nsp2.mis)
            ### weights store PCA scores for sp1 andenvironments 1,2 and both environments combined
            ### specific exact locations of environment and sp. data in PCA for storage run Create a grid
            ### of Environmental Space Function
            z2.ran <- humboldt.grid.espace(scores.env12a[1:2], scores.env2a[1:2], sp2rPCA, R = Rin)
            
        }
        # store simulation values of overlap
        o.ot <- humboldt.niche.overlap(z1, z2.ran, correct.env = inE, nae = "YES")
        # overlap between random and observed niches
        sim.o$D[k] <- o.ot$D  # storage of overlaps
        obs.o$D[k] <- o.i$D
        sim.o$I[k] <- o.ot$I
        obs.o$I[k] <- o.i$I
        
        sim.o$n.pts[k] <- nsp2r
        sim.o$x.shift[k] <- Xran
        sim.o$y.shift[k] <- Yran 
      }
    }
    
    if (force.equal.sample==F){
    for (k in 1:rep) {
        if (run.silent.bak == F) 
            {
                plot.new()
                text(0.5, 0.5, paste("Background tests:", "\n", "runs to go:", rep - k + 1))
            }  # countdown
        {
            Xran <- sample(Rana, 1)
            Yran <- sample(Rana, 1)
            sp2rX <- scores.sp2a[, 3] + Xran
            sp2rY <- scores.sp2a[, 4] + Yran
            sp2r <- cbind(sp2rX, sp2rY)
            sp2r <- as.data.frame.matrix(sp2r)
            sp2rPCA <- na.exclude(humboldt.sample.spp(dfsp = sp2r, colspxy = 1:2, colspkept = NULL, 
                dfvar = scores.env2c, colvarxy = 1:2, colvar = 3:4, resolution = env.reso.sim, 
                run.silent = T))
            nsp2r <- nrow(sp2rPCA)
            ### weights store PCA scores for sp1 andenvironments 1,2 and both environments combined
            ### specific exact locations of environment and sp. data in PCA for storage run Create a grid
            ### of Environmental Space Function
            if (nsp2r < 5){            
                Xran <- sample(Rana, 1)
                Yran <- sample(Rana, 1)
                sp2rX <- scores.sp2a[, 3] + Xran
                sp2rY <- scores.sp2a[, 4] + Yran
                sp2r <- cbind(sp2rX, sp2rY)
                sp2r <- as.data.frame.matrix(sp2r)
                sp2rPCA <- na.exclude(humboldt.sample.spp(dfsp = sp2r, colspxy = 1:2, colspkept = NULL, 
                dfvar = scores.env2c, colvarxy = 1:2, colvar = 3:4, resolution = env.reso.sim, 
                run.silent = T))
                nsp2r <- nrow(sp2rPCA)}
            if (nsp2r < 5){            
                Xran <- sample(Rana, 1)
                Yran <- sample(Rana, 1)
                sp2rX <- scores.sp2a[, 3] + Xran
                sp2rY <- scores.sp2a[, 4] + Yran
                sp2r <- cbind(sp2rX, sp2rY)
                sp2r <- as.data.frame.matrix(sp2r)
                sp2rPCA <- na.exclude(humboldt.sample.spp(dfsp = sp2r, colspxy = 1:2, colspkept = NULL, 
                dfvar = scores.env2c, colvarxy = 1:2, colvar = 3:4, resolution = env.reso.sim, 
                run.silent = T))}
                nsp2r <- nrow(sp2rPCA)
            if (nsp2r > 4) {
                z2.ran <- humboldt.grid.espace(scores.env12a[1:2], scores.env2a[1:2], sp2rPCA, 
                  R = Rin)
                }
        }
        # store simulation values of overlap
        if (nsp2r > 4) {
            o.ot <- humboldt.niche.overlap(z1, z2.ran, correct.env = inE, nae = "YES")
            # overlap between random and observed niches
            sim.o$D[k] <- o.ot$D  # storage of overlaps
            obs.o$D[k] <- o.i$D
            sim.o$I[k] <- o.ot$I
            obs.o$I[k] <- o.i$I
            sim.o$n.pts[k] <- nsp2r
            sim.o$x.shift[k] <- Xran
            sim.o$y.shift[k] <- Yran
            }
      }
    }
    if (run.silent.bak == F) {
        dev.off()
    }
    
    sim.o$nDg[obs.o$D <= sim.o$D] <- 0
    sim.o$nDg[obs.o$D > sim.o$D] <- 1
    sim.o$nIg[obs.o$I <= sim.o$I] <- 0
    sim.o$nIg[obs.o$I > sim.o$I] <- 1
    sim.o$nDl[obs.o$D > sim.o$D] <- 0
    sim.o$nDl[obs.o$D <= sim.o$D] <- 1
    sim.o$nIl[obs.o$I > sim.o$I] <- 0
    sim.o$nIl[obs.o$I <= sim.o$I] <- 1
    obs.o <- obs.o[complete.cases(obs.o), ]
    sim.o <- sim.o[complete.cases(sim.o), ]
    l$sim <- sim.o  # storage
    l$obs <- obs.o[1, ]  # storage
    l$p.D <- ((min((sum(sim.o$nDg)), (sum(sim.o$nDl))) * 2) + 1)/(length(sim.o$D) + 1)  # storage of p-values
    l$p.I <- ((min((sum(sim.o$nIg)), (sum(sim.o$nIl))) * 2) + 1)/(length(sim.o$I) + 1)  # storage of p-values
    #l$p.D<- (sum(sim.o$nDl) + 1)/(length(sim.o$D) + 1) # storage of p-values, on tail Obs> Sim
    #l$p.I<- (sum(sim.o$nIl) + 1)/(length(sim.o$I) + 1) # storage of p-values

    return(l)
}
##################################################################################################
#' Plot niche
#' @param z plots a grid of the density of a species occurrence in environmental space output from humboldt.grid.espace
#' @param title title of graph
#' @param name.axis1 name of axis 1
#' @param name.axis2 name of axis 2
#' @param correct.env If correct.env=T, the analysis corrects occurrence densities of each species by the prevalence of the environments in their range. If correct.env=T, the overlap measure does not correct occurrence densities of each species by the prevalence of the environments in their range.
#' @return Plots the density of a species occurrence in environmental space output from humboldt.grid.espace. Graph attributes include a string value for: title,name.axis1 and name.axis2 to be included in the plot.
#' @seealso \code{humboldt.sample.spp,humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @importFrom scales alpha
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run create a grid of Environmental Space Function
#' z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
#' #'
#' ## plot niche in espace
#' humboldt.plot.niche(z1,"Species 1","PC1","PC2")
#' humboldt.plot.niche(z2,"Species 2","PC1","PC2")

humboldt.plot.niche <- function(z, title = "", name.axis1 = "PC1", name.axis2 = "PC2", correct.env = F) {
    if (correct.env == F) 
        image(z$x, (sort((z$y * -1))), z$z.uncor, col = colorRampPalette(c("#FFFFFF", "#000099", 
            "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256), zlim = c(1e-06, max(z$z.uncor)), 
            xlab = name.axis1, ylab = name.axis2)
    if (correct.env == T) 
        image(z$x, (sort((z$y * -1))), z$z.cor, col = colorRampPalette(c("#FFFFFF", "#000099", 
            "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256), zlim = c(1e-06, max(z$z.cor)), 
            xlab = name.axis1, ylab = name.axis2)
    if (correct.env == F) 
        contour(z$x, (sort((z$y * -1))), z$z.uncor, add = T, levels = quantile(z$z.uncor[z$z.uncor > 
            0], c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)), 
            drawlabels = F, col = alpha("white", 0.75))
    if (correct.env == T) 
        contour(z$x, (sort((z$y * -1))), z$z.cor, add = T, levels = quantile(z$z.cor[z$z.cor > 
            0], c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)), 
            drawlabels = F, col = alpha("white", 0.75))
    contour(z$x, (sort((z$y * -1))), z$Z, add = T, levels = quantile(z$Z[z$Z > 0], c(0.1, 0.5, 
        0.75)), drawlabels = F, lty = c(1, 2, 3), lwd = c(1, 1, 1))
    title(title)
}

###################################################################################################
##################################################################################################
#' Plot PCA contribution to environmental space
#' @param contrib The contribution of the environmental variables included in PCA
#' @param eigen Eigen values output from ordinations
#' @return Plots the contribution of the environmental variables to the analysis. Typically these are eigen vectors and eigen values in ordinations.
#' @seealso \code{humboldt.sample.spp,humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @importFrom ade4 s.corcircle
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y' 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ## plot pca contributions
#' humboldt.plot.contrib(zz$pca.cal$co,zz$pca.cal$eig)
humboldt.plot.contrib <- function(contrib = pca.cal$co, eigen = pca.cal$eig) {
    
    if (ncol(contrib) == 1) {
        h <- c(unlist(contrib))
        n <- row.names(contrib)
        barplot(h, space = 0, names.arg = n)
        title(main = "variable contribution")
    }
    if (ncol(contrib) == 2) {
        s.corcircle(contrib[, 1:2]/max(abs(contrib[, 1:2])), grid = F)
        title(main = "Correlation Circle", sub = paste("PC1 = ", round(eigen[1]/sum(eigen) * 
            100, 2), "%", "PC2 = ", round(eigen[2]/sum(eigen) * 100, 2), "%"))
    }
}

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#' Plot equivalency or background tests: simulation histogram
#' @param x Output from niche equivalency test (humboldt.equivalency.test).
#' @param type Values "D" or "I", which correspond to niche overlap using Schoener's D and Hellinger's I measurements
#' @param title Title of plot
#' @return Plots the results output from the equivalency test (humboldt.equivalency.test).
#' @seealso \code{humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-na.exclude(read.delim("env1.txt",h=T,sep="\t"))
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-na.exclude(read.delim("env2.txt",h=T,sep="\t")) 
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace scores for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run create a grid of Environmental Space Function
#' z1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' z2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
#' 
#' ## run equivalency test
#' a<- humboldt.equivalency.test(z1,z2,rep=100,kern.smooth=1)
#' 
#' ##plot equivalency test
#' humboldt.plot.histrogram(a,"D","Equivalency") 

humboldt.plot.histrogram <- function(x, type = "D", title = "equivalency") {
    if (type == "D") {
        obs <- x$obs$D
        sim <- x$sim$D
        p <- x$p.D
    }
    if (type == "I") {
        obs <- x$obs$I
        sim <- x$sim$I
        p <- x$p.I
    }
    r0 <- c(sim, obs)
    l0 <- max(sim) - min(sim)
    w0 <- l0/(log(length(sim), base = 2) + 1)
    xlim0 <- range(r0) + c(-w0, w0)
    if (xlim0[1] < 0) {
        xlim0 <- c(0, xlim0[2])
    }
    if (xlim0[2] > 1) {
        xlim0 <- c(xlim0[1], 1)
    }
    h0 <- hist(sim, plot = FALSE, nclass = 10)
    y0 <- max(h0$counts)
    hist(sim, plot = TRUE, nclass = 4, xlim = xlim0, col = rgb(0, 0, 1, 1/4), main = title, xlab = type, 
        sub = paste("p.value = ", round(p, 5)))
    lines(c(obs, obs), c(y0/2, 0), col = "red")
    points(obs, y0/2, pch = 18, cex = 2, col = "red")
    invisible()
}

##################################################################################################
##################################################################################################
##################################################################################################
#' Plot equivalency or background test: simulation density curve
#' @param x Output from niche equivalency test (humboldt.equivalency.test).
#' @param type Values "D" or "I", which correspond to niche overlap using Schoener's D and Hellinger's I measurements
#' @param title Title of plot
#' @return Plots the results from the Background test (humboldt.background.test).
#' @seealso \code{humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##perform background tests 
#' bg.sp1tosp2<-humboldt.background.test(g2e=zz, rep = 100, sim.dir = 1, env.reso, kern.smooth = 1, correct.env = F, R = 100, run.silent.bak = F)
#' bg.sp2tosp1<-humboldt.background.test(g2e=zz, rep = 100, sim.dir = 2, env.reso, kern.smooth = 1, correct.env = F, R = 100, run.silent.bak = F)
#' 

#' ##plot background tests 
#' humboldt.plot.density(bg.sp1tosp2,"D","Background 1->2") 
#' humboldt.plot.density(bg.sp2tosp1,"D","Background 2->1") 

humboldt.plot.density <- function(x, type = "D", title = "Background") {
    if (type == "D") {
        obs <- x$obs$D
        sim <- x$sim$D
        p <- x$p.D
    }
    if (type == "I") {
        obs <- x$obs$I
        sim <- x$sim$I
        p <- x$p.I
    }
    # densObs <- density(obs)
    densSim <- density(sim, adjust = 0.5)
    # xlim <- c(0,1)
    r0 <- c(sim, obs)
    l0 <- max(sim) - min(sim)
    w0 <- l0/(log(length(sim), base = 2) + 1)
    xlim <- range(r0) + c(-w0, w0)
    if (xlim[1] < 0) {
        xlim <- c(0, xlim[2])
    }
    if (xlim[2] > 1) {
        xlim <- c(xlim[1], 1)
    }
    ylim <- range(0, densSim$y, obs)
    # pick the colours obsCol <- rgb(1,0,0,0.2)
    simCol <- rgb(0, 0, 1, 0.2)
    ## plot the carrots and set up most of the plot parameters
    plot(densSim, xlim = xlim, ylim = ylim, main = title, xlab = type, sub = paste("p.value = ", 
        round(p, 5)))
    # put our density plots in polygon(densObs, density = -1, col = obsCol)
    polygon(densSim, density = -1, col = simCol)
    y0 <- max(densSim$y, obs)
    lines(c(obs, obs), c(y0 * 0.7, 0), col = "red")
    points(obs, y0 * 0.7, pch = 18, cex = 2, col = "red")
    invisible()
}



##################################################################################################
##################################################################################################
##################################################################################################
#' Rarefy occurrence points
#' @param df input data frame
#' @param colxy columns corresponding to longitude and latitude
#' @param colvar columns of environmental data to include in output, if left blank (=NULL), all will be included. (default=NULL)
#' @param rarefy.dist distance to rarefy points (values need to be in km[recommended] or decimal degrees).  See associated parameter rarefy.units.
#' @param rarefy.units the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees
#' @return A script to systematically select localities within a specified area at specified spatial resolution.  The outcome is always the same and is not random.  This reduces sampling biases in downstream analyses- you should do it! Output is a reduced dataset with less spatial autocorrelation.
#' @importFrom spatstat nndist
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#'
#' ##remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
#' occ.sp1<-humboldt.occ.rarefy(in.pts=occ.sp1,colxy=2:3, rarefy.dist = 50, rarefy.units = "km")
#' occ.sp2<-humboldt.occ.rarefy(in.pts=occ.sp2,colxy=2:3, rarefy.dist = 50, rarefy.units = "km")

humboldt.occ.rarefy <- function(in.pts, colxy = 2:3, rarefy.dist = 0, rarefy.units = "km") {
    if (rarefy.units == "km") {
        min.dist <- rarefy.dist * 1000  #values in km
        sp2p1 <- SpatialPoints(in.pts[, colxy], CRS("+proj=longlat +datum=WGS84"))
        sp2p2 <- spTransform(sp2p1, CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
        xy <- data.frame(x = coordinates(sp2p2)[, 1], y = coordinates(sp2p2)[, 2])
    }
    
    if (rarefy.units == "dd") {
        yy <- colxy[2]
        maxLat <- max(in.pts[, yy])
        minLat <- min(in.pts[, yy])
        # estimage dd to km for study area
        rare.dd <- (mean(c(maxLat, minLat)))
        adjKm <- (-0.0139 * (rare.dd * rare.dd)) + (0.0898 * rare.dd) + 111.1
        min.dist <- adjKm * rarefy.dist * 1000  #values in km
        print(paste("Value used for rarefying:", round((min.dist/1000), 2), "km. Remember that the length of a decimal degrees changes latitudinally due to the convergence of the lines of longitude at the poles. The value used here is the average distance of decimal-degrees within your study area. Alternatively, simply input distance as km value and change rarefy.units='km'"))
        sp2p1 <- SpatialPoints(in.pts[, colxy], CRS("+proj=longlat +datum=WGS84"))
        sp2p2 <- spTransform(sp2p1, CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
        xy <- data.frame(x = coordinates(sp2p2)[, 1], y = coordinates(sp2p2)[, 2])
    }
    
    # setup env- sepaerate from distance
    spName <- in.pts[, 1][1]
    new.data <- NULL
    
    to <- 1
    del.min.dist <- xy
    
    repeat {
        nn1 <- nndist(del.min.dist[, "x"], del.min.dist[, "y"])  # calculate distance nearest neighbour
        if (sum(nn1 < min.dist) == 0) {
            break
        }
        
        # iteratively removing points starting with the one having the minimal distance to the
        # nearest neighbour
        nn2 <- nndist(del.min.dist[, "x"], del.min.dist[, "y"], k = 2)
        
        del1 <- nn1 == min(nn1)
        del2 <- nn2 == min(nn2[del1])
        delk <- del1 & del2
        if (sum(del2) > 1) {
            for (k in 3:6) {
                nn <- nndist(del.min.dist[, "x"], del.min.dist[, "y"], k = k)
                delk <- delk & nn == min(nn[delk])
                if (sum(nn[delk] == min(nn[delk])) > 1) {
                  break
                }
            }
        }
        # from the two points which are the nearest neighbours of the whole set, remove the one
        # closest to the second neighbour
        del.min.dist <- del.min.dist[-(which(delk)[1]), ]
    }
    
    new.data <- rbind(new.data, del.min.dist)
    # }
    finalp1 <- SpatialPoints(new.data[, 1:2], CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    finalp2 <- spTransform(finalp1, CRS("+proj=longlat +datum=WGS84"))
    finalp3 <- data.frame(x = coordinates(finalp2)[, 1], y = coordinates(finalp2)[, 2])
    spNameF <- rep(spName, nrow(finalp3))
    finalp4 <- cbind(spNameF, finalp3)
    print(paste("Starting points =", nrow(xy), ", Final rarefied points =", nrow(new.data)))
    return(finalp4)
}

##################################################################################################
##################################################################################################
##################################################################################################
#' Sample environment to occurrence points
#' @param dfsp Data frame of species with x, y and optional other variables 
#' @param colspxy The range of columns corresponding to longitude and latitude
#' @param colspkept Columns to include in output from species dataset (default="xy")
#' @param dfvar Environmental data frame with x, y and environmental variables
#' @param colvarxy The range of columns for x and y in dfvar
#' @param colvar The range of environmental variables columns in dfvar (excluding "xy", NULL=all)
#' @param resolution of environmental data (typically in decimal degrees)
#' @param run.silent.sam=F  If run.silent.sam=T, texts boxes displaying progress will not be displayed. Default value is FASLE.
#' @return A simple function to sample input environmental data corresponding sites of input occurrence points. Output is a table with environmental variables for each input locality (minus sites with no or missing environmental data).
#' @export
#' @examples
#' library(humboldt)

#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#' 
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y' 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#'
#' ##remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
#' occ.sp1<-humboldt.occ.rarefy(in.pts=occ.sp1,colxy=2:3, rarefy.dist=40, rarefy.units="km")
#' occ.sp2<-humboldt.occ.rarefy(in.pts=occ.sp2,colxy=2:3, rarefy.dist=40, rarefy.units="km")
#'
#' ##sample environment using humboldt.sample.sp() function
#' ## env.reso should be the resolution of the environmental data grid
#' env.reso=0.41666669 
#' occ.sp1<-na.exclude(humboldt.sample.spp(dfsp=occ.sp1,colspxy=2:3,colspkept=NULL,dfvar=env1,colvarxy=1:2,colvar="all",resolution=env.reso))
#' occ.sp2<-na.exclude(humboldt.sample.spp(dfsp=occ.sp2,colspxy=2:3,colspkept=NULL,dfvar=env2,colvarxy=1:2,colvar="all",resolution=env.reso))

humboldt.sample.spp <- function(dfsp, colspxy = 2:3, colspkept = "xy", dfvar, colvarxy = 1:2, 
    colvar = "all", resolution, run.silent.sam = F) {
    
    if (sum(colspkept == "xy") == 1) 
        colspkept <- colspxy
    if (sum(colvar == "all") == 1) {
        if (!is.null(colspkept)) 
            colvar <- (1:ncol(dfvar))[-colvarxy]
        if (is.null(colspkept)) 
            colvar <- (1:ncol(dfvar))
    }
    colspx <- colspxy[1]
    colspy <- colspxy[2]
    colvarx <- colvarxy[1]
    colvary <- colvarxy[2]
    
    x <- dfsp[, colspx]
    X <- dfvar[, colvarx]
    y <- dfsp[, colspy]
    Y <- dfvar[, colvary]
    
    train <- data.frame(matrix(nrow = nrow(dfsp), ncol = length(colvar)))
    names(train) <- names(dfvar)[colvar]
    
    if (run.silent.sam == F) {
        x11(2, 2, pointsize = 12)
        par(mar = c(0, 0, 0, 0))
    }
    for (i in 1:nrow(dfsp)) {
        dist <- sqrt((X - x[i])^2 + (Y - y[i])^2)
        min <- min(dist)
        if (min <= resolution) {
            if (length(colvar) > 1) 
                train[i, ] <- dfvar[dist == min, colvar][1, ]
            if (length(colvar) == 1) 
                train[i, ] <- dfvar[dist == min, colvar][1]
        }
        if (run.silent.sam == F) {
            plot.new()
            text(0.5, 0.5, paste(paste("sampling:", "\n", "runs to go: ", nrow(dfsp) - i)))
        }
    }
    if (run.silent.sam == F) {
        dev.off()
    }
    
    if (!is.null(colspkept)) 
        final <- cbind(dfsp[, colspkept], train)
    if (is.null(colspkept)) 
        final <- train
    return(final)
}


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#' Density scatter plots with histograms on each axis
#' @param x A data frame with two columns of data for plotting 
#' @return A simple function to plot 2 dimensions of a data frame. Warmer colors represent higher densities. On each axis the scatter plot histograms are plotted.
#' @export
#' @examples
#' 
#' library(humboldt)
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' humboldt.plot.scatter(env1[,3:4], xlab="Bio1", ylab="Bio2",main="environment")

humboldt.plot.scatter <- function(x = NA, dcol = "blue", lhist = 20, num.dnorm = 5 * lhist, ...) {
    ## check input
    stopifnot(ncol(x) == 2)
    ## set up layout and graphical parameters
    layMat <- matrix(c(2, 0, 1, 3), ncol = 2, byrow = TRUE)
    layout(layMat, widths = c(5/7, 2/7), heights = c(2/7, 5/7))
    ospc <- 0.5  # outer space
    pext <- 4  # par extension down and to the left
    bspc <- 1  # space between scatter plot and bar plots
    par. <- par(mar = c(pext, pext, bspc, bspc), oma = rep(ospc, 4))  # plot parameters
    # add COLORS TO SCATTER PLOT
    z1 <- x[, 1]
    z2 <- x[, 2]
    df <- data.frame(z1, z2)
    ## Use densCols() output to get density at each point
    z <- densCols(z1, z2, colramp = colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(z)[1, ] + 1L
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
    df$col <- cols[df$dens]
    ## Plot it, reordering rows so that densest points are plotted on top
    plot(z2 ~ z1, data = df[order(df$dens), ], pch = 20, col = col, cex = 2, ...)
    # plot(x, xlim=range(x[,1]), ylim=range(x[,2]),...)  3) determine barplot and height
    # parameter histogram (for barplot-ting the density)
    xhist <- hist(x[, 1], plot = FALSE, breaks = seq(from = min(x[, 1]), to = max(x[, 1]), length.out = lhist))
    yhist <- hist(x[, 2], plot = FALSE, breaks = seq(from = min(x[, 2]), to = max(x[, 2]), length.out = lhist))  # note: this uses probability=TRUE
    ## determine the plot range and all the things needed for the barplots and lines
    xx <- seq(min(x[, 1]), max(x[, 1]), length.out = num.dnorm)  # evaluation points for the overlaid density
    xy <- dnorm(xx, mean = mean(x[, 1]), sd = sd(x[, 1]))  # density points
    yx <- seq(min(x[, 2]), max(x[, 2]), length.out = num.dnorm)
    yy <- dnorm(yx, mean = mean(x[, 2]), sd = sd(x[, 2]))
    ## barplot and line for x (top)
    par(mar = c(0, pext, 0, 0))
    barplot(xhist$density, axes = FALSE, ylim = c(0, max(xhist$density, xy)), space = 0)  # barplot
    lines(seq(from = 0, to = lhist - 1, length.out = num.dnorm), xy, col = dcol)  # line
    ## barplot and line for y (right)
    par(mar = c(pext, 0, 0, 0))
    barplot(yhist$density, axes = FALSE, xlim = c(0, max(yhist$density, yy)), space = 0, horiz = TRUE)  # barplot
    lines(yy, seq(from = 0, to = lhist - 1, length.out = num.dnorm), col = dcol)  # line
    ## restore parameters
    par(par.)
}
##################################################################################################
##################################################################################################
##################################################################################################
#' Plot niche overlap of 2 species in environmental space
#' @param in.pca A single dataset with two PCs for both species (sp1 & sp2) and both environmental spaces (env1 & evn2). Note the input order in the matrix needs to be in this order: sp1, sp2, env1, env2 (seriously- do it!)
#' @param nsp1 the number of records in sp1 dataset
#' @param nsp2 the number of records in sp2 dataset
#' @param nenv1 the number of records in env1 dataset
#' @param nenv2 the number of records in env2 dataset
#' @param pdfname a string for the output pdf of the results
#' @param swap if swap=T, colors will be switched in plots and the order (from top to bottom) will be reversed. The default value is FALSE.
#' @return Plots the overlap of two species' environmental space based on PCs. This code is a derivative work based on a function 'NiceOverPlot' by  Javi Fernández-López and Irene Villa.
#' @return This tool uses the scores of principal components (or 2 environmental variables. Only two dimensions can be used. If performing a PCA, I strongly encourage you to run a species distribution model using lots of environmental data on the focal species and then load only the top contributing environmental variables to be included in the PCA. This way all included variables are known to be relevant in both species distributions. This can be done inside Humboldt (via humboldt.top.env) by importing many environmental variables and letting program select only those important. Alternatively, this can also be done using another method (MaxEnt) outside of R and then use only variables deemed important in the species' distributions 

#' @seealso \code{humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @importFrom ks kde
#' @import ggplot2
#' @importFrom gtable gtable
#' @import grid
#' @import gridExtra
#' @export
#' @examples
#' library(humboldt)
#' ##load environment variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##merge environment files
#' env12<-rbind(env1,env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp','x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp','x','y' 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#'
#' ##remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
#' occ.sp1<-humboldt.occ.rarefy(in.pts=occ.sp1,colxy=2:3, rarefy.dist=40, rarefy.units="km")
#' occ.sp2<-humboldt.occ.rarefy(in.pts=occ.sp2,colxy=2:3, rarefy.dist=40, rarefy.units="km")
#'
#' ## env.reso should be the resolution of the environmental data grid
#' ##sample environment using humboldt.sample.sp() function
#' env.reso=0.166667
#' occ.sp1<-na.exclude(humboldt.sample.spp(dfsp=occ.sp1,colspxy=2:3,colspkept=NULL,dfvar=env1,colvarxy=1:2,colvar="all",resolution=env.reso))
#' occ.sp2<-na.exclude(humboldt.sample.spp(dfsp=occ.sp2,colspxy=2:3,colspkept=NULL,dfvar=env2,colvarxy=1:2,colvar="all",resolution=env.reso))
#' 
#' row.sp1<-(nrow(env1)+nrow(env2)+1):(nrow(env1)+nrow(env2)+nrow(occ.sp1))
#' row.sp2<-(nrow(env1)+nrow(env2)+nrow(occ.sp1)+1):(nrow(env1)+nrow(env2)+nrow(occ.sp1)+nrow(occ.sp2))
#' 
#' ##row weighting of environment density for PCA
#' row.w.1.env<-1-(nrow(env1)/nrow(env12)) # prevalence of env1
#' row.w.2.env<-1-(nrow(env2)/nrow(env12)) # prevalence of env2
#' row.w.env <-c(rep(row.w.1.env, nrow(env1)),rep(row.w.2.env, nrow(env2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
#' 
#' ##create a global dataset with all environments and species, include environmental variables #' e.var[x1-xN]
#' e.var<-c(3:21)
#' data.env.occ2<-rbind(occ.sp1,occ.sp2,env1,env2)[e.var]
#' 
#' #data.env2.occ<-rbind(env1,env2,occ.sp1,occ.sp2)[e.var]
#' row.sp1<-1:nrow(occ.sp1)
#' row.sp2<-(nrow(occ.sp1)+1):(nrow(occ.sp1)+nrow(occ.sp2))
#' row.env1<-(nrow(occ.sp1)+nrow(occ.sp2)+1):(nrow(occ.sp1)+nrow(occ.sp2)+nrow(env1))
#' row.env2<-(nrow(occ.sp1)+nrow(occ.sp2)+nrow(env1)+1):(nrow(occ.sp1)+nrow(occ.sp2)+nrow(env1)+nrow(env2))
#' 
#' pca.cal <-dudi.pca(data.env.occ2, center = T, scale = T, scannf = F, nf = 2)
#' 
#' nSP1<-nrow(sp1<- pca.cal$li[row.sp1,])
#' nSP2<-nrow(sp2<- pca.cal$li[row.sp2,])
#' nCM1<-nrow(env1<- pca.cal$li[row.env1,])
#' nCM2<-nrow(env2<- pca.cal$li[row.env2,])
#' humboldt.plot.overlap(pca.cal, nsp1=nSP1, nsp2=nSP2, nenv1=nCM1, nenv2=nCM2,pdfname="OverlapSp1and2") 

humboldt.plot.overlap <- function(in.pca, nsp1, nsp2, nenv1, nenv2, pdfname = "OverLapPlot.pdf", 
    pdf.out = T, swap = F) {
    ## define inputs
    n1 <- nsp1
    n2 <- nsp2
    n3 <- nenv1
    n4 <- nenv2
    options(warn = -1)
    # prepare the data, depending of the type of input ('pca'/'dudi' object or raw scores)
    plot.axis = TRUE
    bw = NULL
    b = NULL
    a1cont = NULL
    a2cont = NULL
    # specify environmental data from PCA dev.off()
    if (swap == T) {
        sc_1 <- in.pca
        sc_2 <- in.pca
        cm_1 <- in.pca
        cm_2 <- in.pca
        sc2 <- sc_1$li[1:n1, ]
        sc1 <- sc_1$li[(n1 + 1):((n1 + 1) + n2), ]
        cm2 <- sc_1$li[(n2 + 1):((n2 + 1) + n3), ]
        cm1 <- sc_1$li[(n3 + 1):((n3 + 1) + n4), ]
    }
    
    if (swap == F) {
        sc_1 <- in.pca
        sc_2 <- in.pca
        cm_1 <- in.pca
        cm_2 <- in.pca
        sc1 <- sc_1$li[1:n1, ]
        sc2 <- sc_1$li[(n1 + 1):((n1 + 1) + n2), ]
        cm1 <- sc_1$li[(n2 + 1):((n2 + 1) + n3), ]
        cm2 <- sc_1$li[(n3 + 1):((n3 + 1) + n4), ]
    }
    # env1<- in.pca$li[row.env1,] env2<- in.pca$li[row.env2,] calculate kernel densities from
    # environmental data
    kd1 <- ks::kde(cm1, compute.cont = TRUE)
    kd2 <- ks::kde(cm2, compute.cont = TRUE)
    # calculated 99CI for each environmental space
    contour_99env1 <- with(kd1, contourLines(x = eval.points[[1]], y = eval.points[[2]], z = estimate, 
        levels = cont["1%"]))
    contour_99env2 <- with(kd2, contourLines(x = eval.points[[1]], y = eval.points[[2]], z = estimate, 
        levels = cont["1%"]))
    # prepare data for sorting - ONLY can plot a single contour line!
    contour_99env1x <- sapply(contour_99env1, "[", 2)
    contour_99env1y <- sapply(contour_99env1, "[", 3)
    contour_99env2x <- sapply(contour_99env2, "[", 2)
    contour_99env2y <- sapply(contour_99env2, "[", 3)
    
    # choose largest contour created
    contour_99env1xa <- contour_99env1x[order(sapply(contour_99env1x, length), decreasing = T)]
    contour_99env2xa <- contour_99env2x[order(sapply(contour_99env2x, length), decreasing = T)]
    contour_99env1ya <- contour_99env1y[order(sapply(contour_99env1y, length), decreasing = T)]
    contour_99env2ya <- contour_99env2y[order(sapply(contour_99env2y, length), decreasing = T)]
    
    ### merge and name sorted contours
    contour_99env1 <- c(contour_99env1xa[1], contour_99env1ya[1])
    contour_99env2 <- c(contour_99env2xa[1], contour_99env2ya[1])
    contour_99env2 <- data.frame(cbind(contour_99env2$x, contour_99env2$y))
    contour_99env1 <- data.frame(cbind(contour_99env1$x, contour_99env1$y))
    names(contour_99env1) <- c("x", "y")
    names(contour_99env2) <- c("x", "y")
    gc1 <- c(rep(0, nrow(contour_99env1)))
    gc2 <- c(rep(1, nrow(contour_99env2)))
    if (swap == T) {
        gc1 <- c(rep(1, nrow(contour_99env1)))
        gc2 <- c(rep(0, nrow(contour_99env2)))
    }
    # add species info so that plots match
    contour_99env2 <- data.frame(cbind(contour_99env2$x, contour_99env2$y, gc2))
    contour_99env1 <- data.frame(cbind(contour_99env1$x, contour_99env1$y, gc1))
    names(contour_99env1) <- c("x", "y", "g")
    names(contour_99env2) <- c("x", "y", "g")
    options(warn = -1)
    # prepare the data, depending of the type of input ('pca'/'dudi' object or raw scores)
    # recognize both species
    dd <- rbind(sc1, sc2, cm1, cm2)
    AmaxF <- max(dd[, 1])
    AminF <- min(dd[, 1])
    BmaxF <- max(dd[, 1])
    BminF <- min(dd[, 1])
    CmaxF <- max(dd[, 2])
    CminF <- min(dd[, 2])
    DmaxF <- max(dd[, 2])
    DminF <- min(dd[, 2])
    Xmin1 = pmin(AminF, BminF)
    Xmax1 = pmax(AmaxF, BmaxF)
    Ymin1 = pmin(CminF, DminF)
    Ymax1 = pmax(CmaxF, DmaxF)
    
    scores <- rbind(sc1, sc2)
    g <- c(rep(0, nrow(sc1)), rep(1, nrow(sc2)))
    df <- data.frame(cbind(scores$Axis1, scores$Axis2, g))
    names(df) <- c("x", "y", "g")
    df$g <- as.factor(df$g)
    
    # establish an empty plot to be placed at top-right corner (X)
    empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    # sp1
    options(warn = -1)
    p1 <- ggplot(data = df, aes(x, y, color = as.factor(g))) + stat_density2d(aes(fill = ..level..), 
        alpha = 0.15, bins = b, geom = "polygon", h = c(bw, bw)) + scale_fill_continuous(low = "#e94244", 
        high = "#d7191c", space = "Lab", name = "sp1") + scale_colour_discrete(guide = FALSE) + 
        scale_x_continuous(name = "PC1", limits = c(Xmin1, Xmax1)) + scale_y_continuous(name = "PC2", 
        limits = c(Ymin1, Ymax1)) + theme(legend.position = "none") + geom_path(aes(x, y), data = contour_99env2) + 
        geom_path(aes(x, y), data = contour_99env1)
    
    # sp2
    p2 <- ggplot(data = df, aes(x, y, color = as.factor(g))) + stat_density2d(aes(fill = ..level..), 
        alpha = 0.15, bins = b, geom = "polygon", h = c(bw, bw)) + scale_fill_continuous(low = "#6a78dd", 
        high = "#2b3cba", space = "Lab", name = "sp2") + scale_colour_discrete(guide = FALSE) + 
        scale_x_continuous(name = "PC1", limits = c(Xmin1, Xmax1)) + scale_y_continuous(name = "PC2", 
        limits = c(Ymin1, Ymax1)) + theme(legend.position = "none")
    pp1 <- ggplot_build(p1)
    ppp1 <- ggplot_build(p1 + aes(alpha = 0.15) + theme_classic() + theme(legend.position = "none") + 
        theme(text = element_text(size = 15)) + xlab("PC1") + ylab("PC2") + xlim(c(Xmin1 - 0.5, 
        Xmax1 + 0.5)) + ylim(c(Ymin1 - 0.5, Ymax1 + 0.5)))
    pp2 <- ggplot_build(p2 + aes(alpha = 0.15) + theme_classic() + theme(legend.position = "none") + 
        xlab("PC1") + ylab("PC2") + xlim(c(Xmin1 - 0.5, Xmax1 + 0.5)) + ylim(c(Ymin1 - 0.5, Ymax1 + 
        0.5)))$data[[1]]
    ppp1$data[[1]]$fill[grep(pattern = "^2", pp2$group)] <- pp2$fill[grep(pattern = "^2", pp2$group)]
    grob1 <- ggplot_gtable(ppp1)
    grob2 <- ggplotGrob(p2)
    grid.newpage()
    grid.draw(grob1)
    
    # marginal density of x - plot on top
    plot_top <- ggplot(df, aes(x, y = ..scaled.., fill = g)) + geom_density(position = "identity", 
        alpha = 0.5) + scale_x_continuous(name = "PC1", limits = c(Xmin1 - 0.5, Xmax1 + 0.5)) + 
        scale_fill_brewer(palette = "Set1") + theme_classic() + theme(legend.position = "none")
    
    
    # marginal density of y - plot on the right
    
    plot_right <- ggplot(df, aes(y, y = ..scaled.., fill = g)) + geom_density(position = "identity", 
        alpha = 0.5) + scale_x_continuous(name = "PC2", limits = c(Ymin1 - 0.5, Ymax1 + 0.5)) + 
        coord_flip() + scale_fill_brewer(palette = "Set1") + theme_classic() + theme(legend.position = "none")
    
    
    if (plot.axis == TRUE) 
        grid.arrange(plot_top, empty, grob1, plot_right, ncol = 2, nrow = 2, widths = c(4, 1), 
            heights = c(1, 4)) else grid.draw(grob1)
    if (pdf.out == TRUE) 
        ggsave(pdfname, grid.arrange(plot_top, empty, grob1, plot_right, ncol = 2, nrow = 2, 
            widths = c(4, 1), heights = c(1, 4)), device = "pdf", width = 12, height = 12, units = "in")
    # dev.off()
}


#################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#' Convert geographic space to shared environmental space
#' @param env1 environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice.
#' @param env2 environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice.
#' @param sp1 occurrence sites for the species/population 1 at study area 1 (env1). Column names should be 'sp','x','y'
#' @param sp2 occurrence sites for the species/population 2 at study area 2 (env2). Column names should be 'sp','x','y'
#' @param rarefy.dist remove occurrences closer than a minimum distance to each other (this function uses the humboldt.occ.rarefy function). Values need to be in km[recommended] or decimal degrees. See associated parameter rarefy.units. Note: rarefy.dist=0 will remove no occurrences 
#' @param rarefy.units the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees
#' @param env.reso the resolution of the input environmental data grid in decimal degrees
#' @param nae.window the spatial window from which non-analogous environments will be quantified.  The non-analogous environments are characterized by gridding the espace of env1 and env2 into a R x R grind (e.g. 100 x 100).  If nae.window=0, values absent from a cell in one environment will be removed from the other.   If nae.window>0, values absent from a window (or neighborhood of cells) in one environment will be trimmed from the other.  The nae.window value characterizes the number of cells to search from the focal cell of environmental space values in the other environment.  The larger the nae.window value, the fewer non-analogous environments removed. This parameter allows imperfect overlap of environments.   If areas of environmental space are a little patchy between environments—but generally present-- a larger nae.window value will retain more of the patch environments. The default value is a nae.window=5
#' @param reduce.env the format to trim environmental space so that it is shared. If reduce.env=1, the second input environment (env2) will be trim the match the first input (env1). If reduce.env=2, both input environments trimmed so that extents of both are identical (the lower maximum value observed in env1 and env2 and the higher minimum value observed in env1 and env2 will be used to trim environmental space for each PC/environmental variable) If reduce.env=0, you will skip trimming environmental space
#' @param reductype if reduce.env= 1 or 2, the 'reducetype' parameter specifies the format for how to reduce environmental space ("PCA" or "STANDARD"). If reductype="PCA", the environmental space will be trimmed based on two principal components. If reductype="STANDARD", the environmental space will be trimmed by each included variable specified in col.env. if reduce.env=0, do not include this parameter
#' @param non.analogous.environments allow non-analogous environment in environmental space? If non.analogous.environments="YES" non-analogous environments between env1 and env2 will retained. If non.analogous.environments="NO" non-analogous environments between env1 and env2 will be removed. This parameter is only usable under the combinations of reductype="PCA" & reduce.env=1 or reductype="PCA" & reduce.env=2
#' @param env.trim.mcp trim extent of environmental data by buffered MCP in geographic space. Necessary for comparing different species or different extents, and for removing non-accesible environments
#' @param mcp.buffer.sp1 buffer distance (in km) for trimming available environmental space for sp1
#' @param mcp.buffer.sp2 buffer distance (in km) for trimming available environmental space for sp2
#' @param col.env if reductype="STANDARD", then parameter specifies the number of columns to trim environmental space on. This can be any number of columns
#' @param PROJ If PROJ =F, the models are calibrated on both ranges. If PROJ =T, the models are calibrated on species 1 range only and projected to range 2
#' @param e.var selection of variables to include in all of the analyses- this is a separate parameter than col.env, but must contain all variables included in col.env 
#' @param R Resolution of grid in environmental space (RxR)
#' @param thresh.espace.z this parameter is an experimental parameter and controls the level at which values below the kernel density z values are removed for creating areas of analogous environmental space. Higher values will increase value from which the low-density areas are removed from the environmental space of z1 and z2.  Basically values above this are retained and values below are removed. Default=0.001
#' @param kern.smooth scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making a espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Note this value is only used in this script as a secondary test of espace trimming efficiency. Default value is 1
#' @param env.pca calibrate PCA on only occurrence data (env.pca=F) or environment data (env.pca=T), recommended is env.pca=T
#' @param run.silent If run.silent=T, texts boxes displaying progress will not be displayed
#' @return This function converts geographic distributions in environmental space (e-space). Several factors need to be considered when doing this. Do you include or exclude non-analogous environments? Do you trim environment space so that only shared environmental space is being analyzed? These factors directly effects the questions asked in humboldt.equivalency.test & humboldt.background.test.  For example, if run on full environment input for background tests and equivalency test, you are testing the total equivalence between species (or divergence) in current distributions. However if ran with trimmed, shared espace. The equivalency tests for niche evolution or niche divergence (as these measures should only be evaluated in shared analogous espace)

#' @seealso \code{humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @importFrom adehabitatHR mcp
#' @importFrom raster buffer
#' @import sp
#' @importFrom ade4 dudi.pca
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be species,x,y)
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be species,x,y). 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' zz=humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)

humboldt.g2e <- function(env1, env2, sp1, sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", nae.window=5,
    env.trim.mcp = T, e.var, col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, 
     rarefy.dist  = 0, rarefy.units = "km" , env.reso, kern.smooth = 1, PROJ = F, R = 100, run.silent = F) {
    if (reduce.env == 0) {
        REDUC = 0
    }
    if (reduce.env == 1 & reductype == "STANDARD") {
        REDUC = 1
    }
    if (reduce.env == 2 & reductype == "STANDARD") {
        REDUC = 2
    }
    if (reduce.env == 1 & reductype == "PCA" & non.analogous.environments == "YES") {
        REDUC = 3
    }
    if (reduce.env == 2 & reductype == "PCA" & non.analogous.environments == "YES") {
        REDUC = 4
    }
    if (reduce.env == 2 & reductype == "PCA" & non.analogous.environments == "NO") {
        REDUC = 5
    }
    if (reduce.env == 1 & reductype == "PCA" & non.analogous.environments == "NO") {
        REDUC = 6
    }
    
    l <- list()
    ########################################################################### 
    kern.smoothinZ <- kern.smooth
    env1FULL <- env1
    env2FULL <- env2
    
    ###########################################################################
    ## consider if data are rarefied
    if (rarefy.units=="dd"){rde <- rarefy.dist/env.reso}
    if (rarefy.units=="km"){env.reso1<-rarefy.dist/111;rde <- env.reso1/env.reso}
    if (rde <= 0.97 & rarefy.units=="dd" & env.reso >= env.reso) {
        print(paste("Warning!! Please consider increasing the spatial rarefying distance. Your input enviromenment resolution is ", 
            env.reso, "decimal degrees, and a rarefy distance=", rarefy.dist, " decimal degrees."))}
    if (rde <= 0.97 & rarefy.units=="dd" & env.reso < env.reso) {
        print(paste("Warning!! Please consider spatially rarefying your input occurrence localities. Equivalency and background tests can result in type I errors if spatial data have not been properly rarefied. In most situations, data should be spatially rarefied at a distance of 10-80km (~0.0836-0.6667 decimal degrees).  If you have done this outside of Humboldt- ignore this warning. Input: environment resolution=", 
            env.reso, "decimal degrees, and a rarefy distance=", rarefy.dist, " decimal degrees."))}
	if (rde <= 0.97 & rarefy.units=="km" & env.reso >= env.reso) {
        print(paste("Warning!! Please consider increasing the spatial rarefying distance. Your input enviromenment resolution is ", 
            env.reso, "decimal degrees, and a rarefy distance=", rarefy.dist, " km (~",env.reso1," decimal degrees)"))}
	if (rde <= 0.97 & rarefy.units=="km" & env.reso < env.reso) {
        print(paste("Warning!! Please consider spatially rarefying your input occurrence localities. Equivalency and background tests can result in type I errors if spatial data have not been properly rarefied. In most situations, data should be spatially rarefied at a distance of 10-80km (~0.0836-0.6667 decimal degrees).  If you have done this outside of Humboldt- ignore this warning. Input: environment resolution=", 
            env.reso, "decimal degrees, and a rarefy distance=", rarefy.dist, " km (~",env.reso1," decimal degrees)"))}

    ###########################################################################
    ##################Trim the occurrence locs by buffered MCP#################
    ###########################################################################

    if (env.trim.mcp == T) {
        options(warn = -1)
        
        CmaxF <- max(env1[, 2])
        CminF <- min(env1[, 2])
        DmaxF <- max(env2[, 2])
        DminF <- min(env2[, 2])
        
        # estimage dd to km for study area
        Env1.dd.y <- mean(c(CmaxF,CminF))
        Env2.dd.y <- mean(c(DmaxF,DminF))
        
        ### Jason's dirty dd to km for lat
        AvgKm1 <- abs(((-0.0139 * (Env1.dd.y * Env1.dd.y)) + (0.0898 * Env1.dd.y) + 111.1))
        AvgKm2 <- abs(((-0.0139 * (Env2.dd.y * Env2.dd.y)) + (0.0898 * Env2.dd.y) + 111.1))
        
        # Sp1
        env1in <- env1
        mcp.buffer.sp1v <- mcp.buffer.sp1/AvgKm1  #values in km
        mcp.buffer.sp2v <- mcp.buffer.sp2/AvgKm2  #values in km
        OcSp1T <- SpatialPoints(sp1[, 2:3], CRS("+proj=longlat +datum=WGS84"))
        OcSp1Tmcp <- mcp(OcSp1T, percent = 100)
        OcSp1Tmcp <- spTransform(OcSp1Tmcp, CRS("+proj=longlat +datum=WGS84"))
        OcSp1TmcpB <- raster::buffer(OcSp1Tmcp, width = mcp.buffer.sp1v, dissolve = T)
        env1Tmcp <- SpatialPoints(env1in, CRS("+proj=longlat +datum=WGS84"))
        env1sdf <- env1Tmcp[OcSp1TmcpB, ]
        env1 <- as.data.frame(env1sdf)
        # Sp2
        env2in <- env2
        OcSp2T <- SpatialPoints(sp2[, 2:3], CRS("+proj=longlat +datum=WGS84"))
        OcSp2Tmcp <- mcp(OcSp2T, percent = 100)
        OcSp2Tmcp <- spTransform(OcSp2Tmcp, CRS("+proj=longlat +datum=WGS84"))
        OcSp2TmcpB <- raster::buffer(OcSp2Tmcp, width = mcp.buffer.sp2v, dissolve = T)
        env2Tmcp <- SpatialPoints(env2in, CRS("+proj=longlat +datum=WGS84"))
        env2sdf <- env2Tmcp[OcSp2TmcpB, ]
        env2 <- as.data.frame(env2sdf)
        options(warn = 0)
    }
    # env2
    
    ########################################################################### 
    if (REDUC == 1) {
        for (i in col.env) {
            y <- i
            Amax <- max(env1[, y])
            Bmin <- min(env1[, y])
            Cmax <- max(env2[, y])
            Dmin <- min(env2[, y])
            env2 <- env2[which(env2[, y] <= Amax & env2[, y] >= Bmin), ]
            print(paste("***Environmental Variable:", y))
            print("BEFORE REDUCTION")
            print(paste("Maximum values: Environment 1", Amax, "and Environment 2", Cmax))
            print(paste("Minimum values: Environment 1", Bmin, "and Environment 2", Dmin))
            Amax2 <- max(env1[, y])
            Bmin2 <- min(env1[, y])
            Cmax2 <- max(env2[, y])
            Dmin2 <- min(env2[, y])
            print(paste("AFTER REDUCTION"))
            print(paste("Maximum values for Environment 2: ", Cmax2))
            print(paste("Minimum values for Environment 2: ", Dmin2))
            print("*******************")
        }
        print(paste(nrow(env2FULL) - nrow(env2), " sites removed from environment 2"))
    }
    if (REDUC == 2) {
        for (i in col.env) {
            y <- i
            Amax <- max(env1[, y])
            Bmin <- min(env1[, y])
            env2 <- env2[which(env2[, y] <= Amax & env2[, y] >= Bmin), ]
            env1 <- env1[which(env1[, y] <= Cmax & env1[, y] >= Dmin), ]
            print(paste("***Environmental Variable:", y))
            print("BEFORE REDUCTION")
            print(paste("Maximum values: Environment 1", Amax, "and Environment 2", Cmax))
            print(paste("Minimum values: Environment 1", Bmin, "and Environment 2", Dmin))
            Amax2 <- max(env1[, y])
            Bmin2 <- min(env1[, y])
            Cmax2 <- max(env2[, y])
            Dmin2 <- min(env2[, y])
            print(paste("AFTER REDUCTION"))
            print(paste("Maximum values for Environment 1", Amax2, "and Environment 2: ", Cmax2))
            print(paste("Minimum values for Environment 1", Bmin2, "and Environment 2: ", Dmin2))
            print("*******************")
        }
        print(paste(nrow(env1FULL) - nrow(env1), " sites removed from environment 1"))
        print(paste(nrow(env2FULL) - nrow(env2), " sites removed from environment 2"))
    }
    
    
    # Environment for both ranges
    env12 <- rbind(env1, env2)
    names(env12)
    ########################################################################### 
    if (rarefy.dist != 0) {
        print("Rarefying occurence data")
        print("Sp1:")
        occ.sp1 <- humboldt.occ.rarefy(in.pts = sp1, colxy = 2:3, rarefy.dist = rarefy.dist, 
            rarefy.units=rarefy.units)
        print("Sp2:")
        occ.sp2 <- humboldt.occ.rarefy(in.pts = sp2, colxy = 2:3, rarefy.dist = rarefy.dist, 
            rarefy.units=rarefy.units)
    }
    if (rarefy.dist == 0) {
        occ.sp1 <- sp1[1:3]
        occ.sp2 <- sp2[1:3]
    }
    # create sp occurrence dataset by adding environmental variables from the global
    # environmental datasets resolution should be the resolution of the environmental data grids
    print("Sampling environmental data to occurence data")
    occ.sp1FULL <- na.exclude(humboldt.sample.spp(dfsp = occ.sp1, colspxy = 2:3, colspkept = NULL, 
        dfvar = env1FULL, colvarxy = 1:2, colvar = "all", resolution = env.reso, run.silent.sam = run.silent))
    occ.sp2FULL <- na.exclude(humboldt.sample.spp(dfsp = occ.sp2, colspxy = 2:3, colspkept = NULL, 
        dfvar = env2FULL, colvarxy = 1:2, colvar = "all", resolution = env.reso, run.silent.sam = run.silent))
    
    occ.sp1 <- na.exclude(humboldt.sample.spp(dfsp = occ.sp1, colspxy = 2:3, colspkept = NULL, 
        dfvar = env1, colvarxy = 1:2, colvar = "all", resolution = env.reso, run.silent.sam = run.silent))
    occ.sp2 <- na.exclude(humboldt.sample.spp(dfsp = occ.sp2, colspxy = 2:3, colspkept = NULL, 
        dfvar = env2, colvarxy = 1:2, colvar = "all", resolution = env.reso, run.silent.sam = run.silent))
    
    #################################################################################################
    ################### row weighting and grouping factors for ade4 functions ######################
    #################################################################################################
    print("Converting G-space to E-space")
    
    # if PROJ = F
    if (PROJ == F) {
        row.w.1.occ <- 1 - (nrow(occ.sp1)/nrow(rbind(occ.sp1, occ.sp2)))  # prevalence of occ1
        row.w.2.occ <- 1 - (nrow(occ.sp2)/nrow(rbind(occ.sp1, occ.sp2)))  # prevalence of occ2
        row.w.occ <- c(rep(0, nrow(env1)), rep(0, nrow(env2)), rep(row.w.1.occ, nrow(occ.sp1)), 
            rep(row.w.2.occ, nrow(occ.sp2)), rep(0, nrow(env1FULL)), rep(0, nrow(env2FULL)), 
            rep(0, nrow(occ.sp1FULL)), rep(0, nrow(occ.sp2FULL)))
        row.w.1.env <- 1 - (nrow(env1)/nrow(env12))  # prevalence of env1
        row.w.2.env <- 1 - (nrow(env2)/nrow(env12))  # prevalence of env2
        row.w.env <- c(rep(row.w.1.env, nrow(env1)), rep(row.w.2.env, nrow(env2)), rep(0, nrow(occ.sp1)), 
            rep(0, nrow(occ.sp2)), rep(0, nrow(env1FULL)), rep(0, nrow(env2FULL)), rep(0, nrow(occ.sp1FULL)), 
            rep(0, nrow(occ.sp2FULL)))
        fac <- as.factor(c(rep(1, nrow(env1)), rep(2, nrow(env2)), rep(1, nrow(occ.sp1)), rep(2, 
            nrow(occ.sp2))))
    }
    # if PROJ = T
    if (PROJ == T) {
        row.w.occ.PROJT <- c(rep(0, nrow(env1)), rep(0, nrow(env2)), rep(1, nrow(occ.sp1)), rep(0, 
            nrow(occ.sp2)))
        row.w.env.PROJT <- c(rep(1, nrow(env1)), rep(0, nrow(env2)), rep(0, nrow(occ.sp1)), rep(0, 
            nrow(occ.sp2)))
    }
    # global dataset for the analysis and rows for each sub dataset
    data.env.occ <- rbind(env1, env2, occ.sp1, occ.sp2, env1FULL, env2FULL, occ.sp1FULL, occ.sp2FULL)[e.var]
    data.xy.occ <- rbind(env1, env2, occ.sp1, occ.sp2, env1FULL, env2FULL, occ.sp1FULL, occ.sp2FULL)[1:2]
    # data.env2.occ<-rbind(env1,env2,occ.sp1,occ.sp2)[e.var]
    row.env1 <- 1:nrow(env1)
    row.env2 <- (nrow(env1) + 1):(nrow(env1) + nrow(env2))
    row.env12 <- 1:(nrow(env1) + nrow(env2))
    row.sp1 <- (nrow(env1) + nrow(env2) + 1):(nrow(env1) + nrow(env2) + nrow(occ.sp1))
    row.sp2 <- (nrow(env1) + nrow(env2) + nrow(occ.sp1) + 1):(nrow(env1) + nrow(env2) + nrow(occ.sp1) + 
        nrow(occ.sp2))
    row.env1FULL <- (nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + 1):(nrow(env1) + 
        nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + nrow(env1FULL))
    row.env2FULL <- (nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + nrow(env1FULL) + 
        1):(nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + nrow(env1FULL) + nrow(env2FULL))
    row.sp1FULL <- (nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + nrow(env1FULL) + 
        nrow(env2FULL) + 1):(nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + nrow(env1FULL) + 
        nrow(env2FULL) + nrow(occ.sp1FULL))
    row.sp2FULL <- (nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + nrow(env1FULL) + 
        nrow(env2FULL) + nrow(occ.sp1FULL) + 1):(nrow(env1) + nrow(env2) + nrow(occ.sp1) + nrow(occ.sp2) + 
        nrow(env1FULL) + nrow(env2FULL) + nrow(occ.sp1FULL) + nrow(occ.sp2FULL))
    #################################################################################################
    #################################### PCA     ####################################################
    #################################################################################################
    # fit of the analyse using occurrences from both ranges
    if (env.pca == T & PROJ == F) {
        row.w.IN <- row.w.env
    }
    if (env.pca == T & PROJ == T) {
        # fit of the analyse using occurrences from range 1
        row.w.IN <- row.w.env.PROJT
    }
    if (env.pca == F & PROJ == F) {
        # fit of the analyse using occurrences from both ranges
        row.w.IN <- row.w.occ
    }
    if (env.pca == F & PROJ == T) {
        # fit of the analyse using occurrences from range
        row.w.IN <- row.w.occ.PROJT
    }
    
    ## approrpiate PCA
    pca.cal <- dudi.pca(data.env.occ, row.w = row.w.IN, center = T, scale = T, scannf = F, nf = 2)
    # predict the scores on the axes
    
    scores.env12 <- pca.cal$li[row.env12, ]
    scores.env1 <- pca.cal$li[row.env1, ]
    scores.env2 <- pca.cal$li[row.env2, ]
    scores.sp1 <- pca.cal$li[row.sp1, ]
    scores.sp2 <- pca.cal$li[row.sp2, ]
    scores.env1FULL <- pca.cal$li[row.env1FULL, ]
    scores.env2FULL <- pca.cal$li[row.env2FULL, ]
    scores.sp1FULL <- pca.cal$li[row.sp1FULL, ]
    scores.sp2FULL <- pca.cal$li[row.sp2FULL, ]
    ### CURRENT
    scores.sp1 <- cbind(pca.cal$li[row.sp1, ], data.xy.occ[row.sp1, ])
    scores.sp2 <- cbind(pca.cal$li[row.sp2, ], data.xy.occ[row.sp2, ])
    scores.env1 <- cbind(pca.cal$li[row.env1, ], data.xy.occ[row.env1, ])
    scores.env2 <- cbind(pca.cal$li[row.env2, ], data.xy.occ[row.env2, ])
    scores.sp1FULL <- cbind(pca.cal$li[row.sp1FULL, ], data.xy.occ[row.sp1FULL, ])
    scores.sp2FULL <- cbind(pca.cal$li[row.sp2FULL, ], data.xy.occ[row.sp2FULL, ])
    scores.env1FULL <- cbind(pca.cal$li[row.env1FULL, ], data.xy.occ[row.env1FULL, ])
    scores.env2FULL <- cbind(pca.cal$li[row.env2FULL, ], data.xy.occ[row.env2FULL, ])
    
    # flag same environments
    z3 <- humboldt.grid.espace(glob = scores.env1[1:2], glob1 = scores.env1[1:2], sp = scores.env1[1:2], 
        kern.smooth = kern.smoothinZ, R = R)
    z4 <- humboldt.grid.espace(glob = scores.env2[1:2], glob1 = scores.env2[1:2], sp = scores.env2[1:2], 
        kern.smooth = kern.smoothinZ, R = R)
    DvalClimP <- round(as.numeric(humboldt.niche.overlap(z3, z4, correct.env = F)[1]), 3)
    if (DvalClimP > 0.97 & reduce.env == 1 & reductype == "PCA") {
        REDUC = 3
    }
    if (DvalClimP > 0.97 & reduce.env == 2 & reductype == "PCA") {
        REDUC = 4
    }
    AmaxPCA <- max(scores.env1[, 1])
    BminPCA <- min(scores.env1[, 1])
    PCAvL <- AmaxPCA - BminPCA
    PCAvLi <- PCAvL/(R * 1)
    
    print("Processing of E-space")
    col.envPCA = c(1:2)
    if (REDUC == 3) {
        for (i in col.envPCA) {
            y <- i
            AmaxPCA <- max(scores.env1[, y])
            BminPCA <- min(scores.env1[, y])
            CmaxPCA <- max(scores.env2[, y])
            DminPCA <- min(scores.env2[, y])
            XminPCA = pmin(DminPCA, BminPCA)
            XmaxPCA = pmax(AmaxPCA, CmaxPCA)
            scores.env2 <- scores.env2[which(scores.env2[, y] <= AmaxPCA & scores.env2[, y] >= 
                BminPCA), ]
            scores.sp2 <- scores.sp2[which(scores.sp2[, y] <= AmaxPCA & scores.sp2[, y] >= BminPCA), 
                ]
            print(paste("***PC Environment Variable:", y))
            print("BEFORE REDUCTION")
            print(paste("Maximum values: Environment 1", AmaxPCA, "and Environment 2", CmaxPCA))
            print(paste("Minimum values: Environment 1", BminPCA, "and Environment 2", DminPCA))
            Amax2PCA <- max(scores.env1[, y])
            Bmin2PCA <- min(scores.env1[, y])
            Cmax2PCA <- max(scores.env2[, y])
            Dmin2PCA <- min(scores.env2[, y])
            print(paste("AFTER REDUCTION"))
            print(paste("Maximum values for Environment 2: ", Cmax2PCA))
            print(paste("Minimum values for Environment 2: ", Dmin2PCA))
            print("*******************")
        }
        print(paste(nrow(scores.env2FULL) - nrow(scores.env2), " sites removed from Environment 2"))
        print(paste(nrow(scores.sp2FULL) - nrow(scores.sp2), " localties removed from Sp 2 dataset"))
        scores.env12 <- rbind(scores.env1, scores.env2)
        AmaxPCA <- max(scores.env1[, 1])
        BminPCA <- min(scores.env1[, 1])
        PCAvL <- AmaxPCA - BminPCA
        PCAvLi <- PCAvL/(R * 1)
    }
    if (REDUC == 4) {
        for (i in col.envPCA) {
            y <- i
            AmaxPCA <- max(scores.env1[, y])
            BminPCA <- min(scores.env1[, y])
            CmaxPCA <- max(scores.env2[, y])
            DminPCA <- min(scores.env2[, y])
            XminPCA = pmin(DminPCA, BminPCA)
            XmaxPCA = pmax(AmaxPCA, CmaxPCA)
            scores.env2 <- scores.env2[which(scores.env2[, y] <= AmaxPCA & scores.env2[, y] >= 
                BminPCA), ]
            scores.sp2 <- scores.sp2[which(scores.sp2[, y] <= AmaxPCA & scores.sp2[, y] >= BminPCA), 
                ]
            scores.env1 <- scores.env1[which(scores.env1[, y] <= CmaxPCA & scores.env1[, y] >= 
                DminPCA), ]
            scores.sp1 <- scores.sp1[which(scores.sp1[, y] <= CmaxPCA & scores.sp1[, y] >= DminPCA), 
                ]
            print(paste("***PC Environment Variable:", y))
            print("BEFORE REDUCTION")
            print(paste("Maximum values: Environment 1", AmaxPCA, "and Environment 2", CmaxPCA))
            print(paste("Minimum values: Environment 1", BminPCA, "and Environment 2", DminPCA))
            Amax2PCA <- max(scores.env1[, y])
            Bmin2PCA <- min(scores.env1[, y])
            Cmax2PCA <- max(scores.env2[, y])
            Dmin2PCA <- min(scores.env2[, y])
            print(paste("AFTER REDUCTION"))
            print(paste("Maximum values for environment 1: ", Amax2PCA))
            print(paste("Minimum values for environment 1: ", Bmin2PCA))
            print(paste("Maximum values for environment 2: ", Cmax2PCA))
            print(paste("Minimum values for environment 2: ", Dmin2PCA))
            print("*******************")
        }
        print(paste("Step1. Reduce environmental space to max/mins of Environment 1:", nrow(scores.env2FULL) - 
            nrow(scores.env2), " sites (from a total of", nrow(scores.env2FULL), ") removed from Environment 2"))
        RedP1C <- nrow(scores.env2FULL) - nrow(scores.env2)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 1:", nrow(scores.sp2FULL) - 
            nrow(scores.sp2), " localities (from a total of", nrow(scores.sp2FULL), ") removed from sp 2 dataset"))
        RedP1sp <- nrow(scores.sp2FULL) - nrow(scores.sp2)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 2:", nrow(scores.env1FULL) - 
            nrow(scores.env1), " sites (from a total of", nrow(scores.env1FULL), ") removed from Environment 1"))
        RedP2C <- nrow(scores.env1FULL) - nrow(scores.env1)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 2:", nrow(scores.sp1FULL) - 
            nrow(scores.sp1), " localities (from a total of", nrow(scores.sp1FULL), ") removed from sp 1 dataset"))
        scores.env12 <- rbind(scores.env1, scores.env2)
        AmaxPCA <- max(scores.env1[, 1])
        BminPCA <- min(scores.env1[, 1])
        PCAvL <- AmaxPCA - BminPCA
        PCAvLi <- PCAvL/(R * 1)
    }
    options(warn = -1)
    ########################################################################################
    if (REDUC == 5) {
        for (i in col.envPCA) {
            y <- i
            AmaxPCA <- max(scores.env1[, y])
            BminPCA <- min(scores.env1[, y])
            CmaxPCA <- max(scores.env2[, y])
            DminPCA <- min(scores.env2[, y])
            XminPCA <- pmin(DminPCA, BminPCA)
            XmaxPCA <- pmax(AmaxPCA, CmaxPCA)
            scores.env2 <- scores.env2[which(scores.env2[, y] <= AmaxPCA & scores.env2[, y] >= 
                BminPCA), ]
            scores.sp2 <- scores.sp2[which(scores.sp2[, y] <= AmaxPCA & scores.sp2[, y] >= BminPCA), 
                ]
            scores.env1 <- scores.env1[which(scores.env1[, y] <= CmaxPCA & scores.env1[, y] >= 
                DminPCA), ]
            scores.sp1 <- scores.sp1[which(scores.sp1[, y] <= CmaxPCA & scores.sp1[, y] >= DminPCA), 
                ]
            print(paste("***PC Environment Variable:", y))
            print("BEFORE REDUCTION")
            print(paste("Maximum values: Environment 1", AmaxPCA, "and Environment 2", CmaxPCA))
            print(paste("Minimum values: Environment 1", BminPCA, "and Environment 2", DminPCA))
            Amax2PCA <- max(scores.env1[, y])
            Bmin2PCA <- min(scores.env1[, y])
            Cmax2PCA <- max(scores.env2[, y])
            Dmin2PCA <- min(scores.env2[, y])
            print(paste("AFTER REDUCTION"))
            print(paste("Maximum values for environment 1: ", Amax2PCA))
            print(paste("Minimum values for environment 1: ", Bmin2PCA))
            print(paste("Maximum values for environment 2: ", Cmax2PCA))
            print(paste("Minimum values for environment 2: ", Dmin2PCA))
            print("*******************")
        }
        print(paste("Step1. Reduce environmental space to max/mins of Environment 1:", nrow(scores.env2FULL) - 
            nrow(scores.env2), " sites (from a total of", nrow(scores.env2FULL), ") removed from Environment 2"))
        RedP1C <- nrow(scores.env2FULL) - nrow(scores.env2)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 1:", nrow(scores.sp2FULL) - 
            nrow(scores.sp2), " localities (from a total of", nrow(scores.sp2FULL), ") removed from sp 2 dataset"))
        RedP1sp <- nrow(scores.sp2FULL) - nrow(scores.sp2)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 2:", nrow(scores.env1FULL) - 
            nrow(scores.env1), " sites (from a total of", nrow(scores.env1FULL), ") removed from Environment 1"))
        RedP2C <- nrow(scores.env1FULL) - nrow(scores.env1)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 2:", nrow(scores.sp1FULL) - 
            nrow(scores.sp1), " localities (from a total of", nrow(scores.sp1FULL), ") removed from sp 1 dataset"))
        RedP2sp <- nrow(scores.sp1FULL) - nrow(scores.sp1)
        AmaxPCA <- max(scores.env1[, 1])
        BminPCA <- min(scores.env1[, 1])
        PCAvL <- AmaxPCA - BminPCA
        PCAvLi <- PCAvL/(R * 1)
        Rana <- seq(BminPCA, AmaxPCA, length.out = (R * 1))
        
        if(nae.window==0){
        for (i in Rana) {
            y <- i
            # scores.env2b <- scores.env2c
            scores.env1a <- scores.env1[which(scores.env1[, 1] <= y + PCAvLi & scores.env1[, 
                1] >= y), ]
            Maxi <- max(scores.env1a[, 2])
            Mini <- min(scores.env1a[, 2])
            scores.env2a <- scores.env2[which(scores.env2[, 1] <= y + PCAvLi & scores.env2[, 
                1] >= y), ]
            scores.env2a <- scores.env2a[which(scores.env2a[, 2] <= Maxi & scores.env2a[, 2] >= 
                Mini), ]
            if (!exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- scores.env2a
            }
            if (exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- rbind(scores.env2a, scores.env2b)
            }
            scores.sp2a <- scores.sp2[which(scores.sp2[, 1] <= y + PCAvLi & scores.sp2[, 1] >= 
                y), ]
            scores.sp2a <- scores.sp2a[which(scores.sp2a[, 2] <= Maxi & scores.sp2a[, 2] >= Mini), 
                ]
            if (!exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- scores.sp2a
            }
            if (exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- rbind(scores.sp2a, scores.sp2b)
            }
        }
        scores.sp2 <- scores.sp2b
        rm(scores.sp2b)
        scores.env2 <- scores.env2b
        rm(scores.env2b)
        }
        if(nae.window!=0){
        for (i in Rana) {
            y <- i
            # scores.env2b <- scores.env2c
            scores.env1a <- scores.env1[which(scores.env1[, 1] <= y + ((nae.window+1)*PCAvLi) & scores.env1[, 
                1] >= y - ((nae.window)*PCAvLi)), ]
            Maxi <- max(scores.env1a[, 2])
            Mini <- min(scores.env1a[, 2])
            scores.env2a <- scores.env2[which(scores.env2[, 1] <= y + PCAvLi & scores.env2[, 
                1] >= y), ]
            scores.env2a <- scores.env2a[which(scores.env2a[, 2] <= Maxi & scores.env2a[, 2] >= 
                Mini), ]
            if (!exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- scores.env2a
            }
            if (exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- rbind(scores.env2a, scores.env2b)
            }
            scores.sp2a <- scores.sp2[which(scores.sp2[, 1] <= y + PCAvLi & scores.sp2[, 1] >= 
                y), ]
            scores.sp2a <- scores.sp2a[which(scores.sp2a[, 2] <= Maxi & scores.sp2a[, 2] >= Mini), 
                ]
            if (!exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- scores.sp2a
            }
            if (exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- rbind(scores.sp2a, scores.sp2b)
            }
        }
        scores.sp2 <- scores.sp2b
        rm(scores.sp2b)
        scores.env2 <- scores.env2b
        rm(scores.env2b)
        }
        print(paste("Step2. Remove non-analogous environments:", nrow(scores.env2FULL) - nrow(scores.env2) - 
            RedP1C, " sites removed from environment 2 (from a total of ", nrow(scores.env2FULL), 
            ")"))
        print(paste("Step2. Remove non-analogous environments:", nrow(scores.sp2FULL) - nrow(scores.sp2) - 
            RedP1sp, " localities removed from sp 2 dataset (from a total of ", nrow(scores.sp2FULL), 
            ")"))
        ######################################################################################################################### 
        AmaxPCA <- max(scores.env2[, 1])
        BminPCA <- min(scores.env2[, 1])
        PCAvL <- AmaxPCA - BminPCA
        PCAvLi <- PCAvL/(R * 1)
        Rana <- seq(BminPCA, AmaxPCA, length.out = (R * 1))
        
        if(nae.window==0){
        for (i in Rana) {
            y <- i
            # scores.env2b <- scores.env2c
            scores.env2a <- scores.env2[which(scores.env2[, 1] <= y + PCAvLi & scores.env2[, 
                1] >= y), ]
            Maxi <- max(scores.env2a[, 2])
            Mini <- min(scores.env2a[, 2])
            scores.env1a <- scores.env1[which(scores.env1[, 1] <= y + PCAvLi & scores.env1[, 
                1] >= y), ]
            scores.env1a <- scores.env1a[which(scores.env1a[, 2] <= Maxi & scores.env1a[, 2] >= 
                Mini), ]
            if (!exists("scores.env1b") & exists("scores.env1a")) {
                scores.env1b <- scores.env1a
            }
            if (exists("scores.env1b") & exists("scores.env1a")) {
                scores.env1b <- rbind(scores.env1a, scores.env1b)
            }
            scores.sp1a <- scores.sp1[which(scores.sp1[, 1] <= y + PCAvLi & scores.sp1[, 1] >= 
                y), ]
            scores.sp1a <- scores.sp1a[which(scores.sp1a[, 2] <= Maxi & scores.sp1a[, 2] >= Mini), 
                ]
            if (!exists("scores.sp1b") & exists("scores.sp1a")) {
                scores.sp1b <- scores.sp1a
            }
            if (exists("scores.sp1b") & exists("scores.sp1a")) {
                scores.sp1b <- rbind(scores.sp1a, scores.sp1b)
            }
        }
        scores.sp1 <- scores.sp1b
        rm(scores.sp1b)
        scores.env1 <- scores.env1b
        rm(scores.env1b)
        }
        if(nae.window!=0){
        for (i in Rana) {
            y <- i
            # scores.env2b <- scores.env2c
            scores.env2a <- scores.env2[which(scores.env2[, 1] <= y + ((nae.window+1)*PCAvLi) & scores.env2[, 
                1] >= y - ((nae.window)*PCAvLi)), ]
            Maxi <- max(scores.env2a[, 2])
            Mini <- min(scores.env2a[, 2])
            scores.env1a <- scores.env1[which(scores.env1[, 1] <= y + PCAvLi & scores.env1[, 
                1] >= y), ]
            scores.env1a <- scores.env1a[which(scores.env1a[, 2] <= Maxi & scores.env1a[, 2] >= 
                Mini), ]
            if (!exists("scores.env1b") & exists("scores.env1a")) {
                scores.env1b <- scores.env1a
            }
            if (exists("scores.env1b") & exists("scores.env1a")) {
                scores.env1b <- rbind(scores.env1a, scores.env1b)
            }
            scores.sp1a <- scores.sp1[which(scores.sp1[, 1] <= y + PCAvLi & scores.sp1[, 1] >= 
                y), ]
            scores.sp1a <- scores.sp1a[which(scores.sp1a[, 2] <= Maxi & scores.sp1a[, 2] >= Mini), 
                ]
            if (!exists("scores.sp1b") & exists("scores.sp1a")) {
                scores.sp1b <- scores.sp1a
            }
            if (exists("scores.sp1b") & exists("scores.sp1a")) {
                scores.sp1b <- rbind(scores.sp1a, scores.sp1b)
            }
        }
        scores.sp1 <- scores.sp1b
        rm(scores.sp1b)
        scores.env1 <- scores.env1b
        rm(scores.env1b)
        }
        
        print(paste("Step2. Remove non-analogous environments:", nrow(scores.env1FULL) - nrow(scores.env1) - 
            RedP2C, " sites removed from environment 1 (from a total of ", nrow(scores.env1FULL), 
            ")"))
        print(paste("Step2. Remove non-analogous environments:", nrow(scores.sp1FULL) - nrow(scores.sp1) - 
            RedP2sp, " localities removed from sp 1 dataset (from a total of ", nrow(scores.sp1FULL), 
            ")"))
        scores.env12 <- rbind(scores.env1, scores.env2)
        # print(nrow(scores.env12))
    }
    if (REDUC == 6) {
        for (i in col.envPCA) {
            
            y <- i
            AmaxPCA <- max(scores.env1[, y])
            BminPCA <- min(scores.env1[, y])
            CmaxPCA <- max(scores.env2[, y])
            DminPCA <- min(scores.env2[, y])
            XminPCA <- pmin(DminPCA, BminPCA)
            XmaxPCA <- pmax(AmaxPCA, CmaxPCA)
            scores.env2 <- scores.env2[which(scores.env2[, y] <= AmaxPCA & scores.env2[, y] >= 
                BminPCA), ]
            scores.sp2 <- scores.sp2[which(scores.sp2[, y] <= AmaxPCA & scores.sp2[, y] >= BminPCA), 
                ]
            print(paste("***PC Environment Variable:", y))
            print("BEFORE REDUCTION")
            print(paste("Maximum values: Environment 1", AmaxPCA, "and Environment 2", CmaxPCA))
            print(paste("Minimum values: Environment 1", BminPCA, "and Environment 2", DminPCA))
            Amax2PCA <- max(scores.env1[, y])
            Bmin2PCA <- min(scores.env1[, y])
            Cmax2PCA <- max(scores.env2[, y])
            Dmin2PCA <- min(scores.env2[, y])
            print(paste("AFTER REDUCTION"))
            print(paste("Maximum values for environment 2: ", Cmax2PCA))
            print(paste("Minimum values for environment 2: ", Dmin2PCA))
            print("*******************")
        }
        print(paste("Step1. Reduce environmental space to max/mins of Environment 1:", nrow(scores.env2FULL) - 
            nrow(scores.env2), " sites (from a total of", nrow(scores.env2FULL), ") removed from Environment 2"))
        RedP1C <- nrow(scores.env2FULL) - nrow(scores.env2)
        print(paste("Step1. Reduce environmental space to max/mins of Environment 1:", nrow(scores.sp2FULL) - 
            nrow(scores.sp2), " localities (from a total of", nrow(scores.sp2FULL), ") removed from sp 2 dataset"))
        RedP1sp <- nrow(scores.sp2FULL) - nrow(scores.sp2)
        AmaxPCA <- max(scores.env1[, 1])
        BminPCA <- min(scores.env1[, 1])
        PCAvL <- AmaxPCA - BminPCA
        PCAvLi <- PCAvL/(R * 1)
        Rana <- seq(BminPCA, AmaxPCA, length.out = (R * 1))
        
        if(nae.window==0){
        for (i in Rana) {
            y <- i
            # scores.env2b <- scores.env2c
            scores.env1a <- scores.env1[which(scores.env1[, 1] <= y + PCAvLi & scores.env1[, 
                1] >= y), ]
            Maxi <- max(scores.env1a[, 2])
            Mini <- min(scores.env1a[, 2])
            scores.env2a <- scores.env2[which(scores.env2[, 1] <= y + PCAvLi & scores.env2[, 
                1] >= y), ]
            scores.env2a <- scores.env2a[which(scores.env2a[, 2] <= Maxi & scores.env2a[, 2] >= 
                Mini), ]
            if (!exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- scores.env2a
            }
            if (exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- rbind(scores.env2a, scores.env2b)
            }
            scores.sp2a <- scores.sp2[which(scores.sp2[, 1] <= y + PCAvLi & scores.sp2[, 1] >= 
                y), ]
            scores.sp2a <- scores.sp2a[which(scores.sp2a[, 2] <= Maxi & scores.sp2a[, 2] >= Mini), 
                ]
            if (!exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- scores.sp2a
            }
            if (exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- rbind(scores.sp2a, scores.sp2b)
            }
        }
        scores.sp2 <- scores.sp2b
        rm(scores.sp2b)
        scores.env2 <- scores.env2b
        rm(scores.env2b)
        }
        if(nae.window!=0){
        for (i in Rana) {
            y <- i
            # scores.env2b <- scores.env2c
            scores.env1a <- scores.env1[which(scores.env1[, 1] <= y + ((nae.window+1)*PCAvLi) & scores.env1[, 
                1] >= y - ((nae.window)*PCAvLi)), ]
            Maxi <- max(scores.env1a[, 2])
            Mini <- min(scores.env1a[, 2])
            scores.env2a <- scores.env2[which(scores.env2[, 1] <= y + PCAvLi & scores.env2[, 
                1] >= y), ]
            scores.env2a <- scores.env2a[which(scores.env2a[, 2] <= Maxi & scores.env2a[, 2] >= 
                Mini), ]
            if (!exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- scores.env2a
            }
            if (exists("scores.env2b") & exists("scores.env2a")) {
                scores.env2b <- rbind(scores.env2a, scores.env2b)
            }
            scores.sp2a <- scores.sp2[which(scores.sp2[, 1] <= y + PCAvLi & scores.sp2[, 1] >= 
                y), ]
            scores.sp2a <- scores.sp2a[which(scores.sp2a[, 2] <= Maxi & scores.sp2a[, 2] >= Mini), 
                ]
            if (!exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- scores.sp2a
            }
            if (exists("scores.sp2b") & exists("scores.sp2a")) {
                scores.sp2b <- rbind(scores.sp2a, scores.sp2b)
            }
        }
        scores.sp2 <- scores.sp2b
        rm(scores.sp2b)
        scores.env2 <- scores.env2b
        rm(scores.env2b)
        }
        print(paste("Step2. Remove non-analogous environments:", nrow(scores.env2FULL) - nrow(scores.env2) - 
            RedP1C, " sites removed from environment 2 (from a total of ", nrow(scores.env2FULL), 
            ")"))
        print(paste("Step2. Remove non-analogous environments:", nrow(scores.sp2FULL) - nrow(scores.sp2) - 
            RedP1sp, " localities removed from sp 2 dataset (from a total of ", nrow(scores.sp2FULL), 
            ")"))
        scores.env12 <- rbind(scores.env1, scores.env2)
        ## recreate background after adjustments
    }
    # print(REDUC)
    RedP1C <- nrow(scores.env2FULL) - nrow(scores.env2)
    RedP1sp <- nrow(scores.sp2FULL) - nrow(scores.sp2)
    RedP2C <- nrow(scores.env1FULL) - nrow(scores.env1)
    RedP2sp <- nrow(scores.sp1FULL) - nrow(scores.sp1)
    
    # scores.env12<- rbind(scores.env1,scores.env2)
    options(warn = 0)
    l$env1 <- env1
    l$env2 <- env2
    l$occ.sp1 <- occ.sp1
    l$occ.sp2 <- occ.sp2
    l$row.env1 <- row.env1
    l$row.env2 <- row.env2
    l$row.env12 <- row.env12
    l$row.sp1 <- row.sp1
    l$row.sp2 <- row.sp2
    l$row.env1FULL <- row.env1FULL
    l$row.env2FULL <- row.env2FULL
    l$row.sp1FULL <- row.sp1FULL
    l$row.sp2FULL <- row.sp2FULL
    l$pca.cal <- pca.cal
    l$scores.env12 <- scores.env12
    l$scores.env1 <- scores.env1
    l$scores.env2 <- scores.env2
    l$scores.sp1 <- scores.sp1
    l$scores.sp2 <- scores.sp2
    l$scores.env1FULL <- scores.env1FULL
    l$scores.env2FULL <- scores.env2FULL
    l$scores.sp1FULL <- scores.sp1FULL
    l$scores.sp2FULL <- scores.sp2FULL
    l$DvalClimP <- DvalClimP
    l$RedP1C <- RedP2C
    l$RedP1sp <- RedP2sp
    l$RedP2C <- RedP1C
    l$RedP2sp <- RedP1sp
    l$AmaxPCA <- AmaxPCA
    l$BminPCA <- BminPCA
    l$PCAvL <- PCAvL
    l$PCAvLi <- PCAvLi
    
    invisible(l)
}



##################################################################################################
##################################################################################################
############################## do it all #################################
#################################################################################################
##################################################################################################
#' 'Do it all'. Performs most analyses in Humboldt and outputs the results as Jason likes!
#' @param inname string for labelling results output from function
#' @param env1 environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice.
#' @param env2 environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice.
#' @param sp1 occurrence sites for the species/population 1 at study area 1 (env1). Column names should be 'sp', 'x','y'
#' @param sp2 occurrence sites for the species/population 2 at study area 2 (env2). Column names should be 'sp', 'x', 'y'
#' @param rarefy.dist remove occurrences closer than a minimum distance to each other (this function uses the humboldt.occ.rarefy function). Values need to be in km[recommended] or decimal degrees. See associated parameter rarefy.units. Note: rarefy.dist=0 will remove no occurrences 
#' @param rarefy.units the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees
#' @param env.reso the resolution of the input environmental data grid in decimal degrees
#' @param reduce.env the format to trim environmental space so that it is shared. If reduce.env=1, the second input environment (env2) will be trim the match the first input (env1). If reduce.env=2, both input environments trimmed so that extents of both are identical (the lower maximum value observed in env1 and env2 and the higher minimum value observed in env1 and env2 will be used to trim environmental space for each PC/environmental variable) If reduce.env=0, you will skip trimming environmental space
#' @param reductype if reduce.env= 1 or 2, the 'reducetype' parameter specifies the format for how to reduce environmental space ("PCA" or "STANDARD"). If reductype="PCA", the environmental space will be trimmed based on two principal components. If reductype="STANDARD", the environmental space will be trimmed by each included variable specified in col.env. If reduce.env=0, do not include this parameter
#' @param non.analogous.environments allow non-analogous environments in environmental space? If non.analogous.environments="YES", non-analogous environments between env1 and env2 will retained. If non.analogous.environments="NO", non-analogous environments between env1 and env2 will be removed. This parameter is only usable under the combinations of reductype="PCA" & reduce.env=1 or reductype="PCA" & reduce.env=2. If something else, this parameter will not be used.
#' @param correct.env if correct.env=T, the analysis corrects occurrence densities of each species by the prevalence of the environments in their range. If correct.env=F, the overlap measure does not correct occurrence densities of each species by the prevalence of the environments in their range.Default=T
#' @param env.trim.mcp trim extent of environmental data by buffered MCP in geographic space. Necessary for comparing different species or different extents and for removing non-accesible environments
#' @param mcp.buffer.sp1 buffer distance (in km) for trimming available environmental space for sp1
#' @param mcp.buffer.sp2 buffer distance (in km) for trimming available environmental space for sp2
#' @param col.env if reductype="STANDARD", then parameter specifies the number of columns to trim environmental space on. This can be any number of columns
#' @param PROJ If PROJ =F, the models are calibrated on both ranges. If PROJ =T, the models are calibrated on species 1 range only and projected to range 2
#' @param e.var selection of variables to include in all of the analyses- this is a separate parameter than col.env but must contain all variables included in col.env 
#' @param R resolution of grid in environmental space (RxR)
#' @param thresh.espace.z this parameter is an experimental parameter and controls the level at which values below the kernel density z values are removed for creating areas of analogous environmental space. Higher values will increase value from which the low-density areas are removed from the environmental space of z1 and z2.  Basical values above this are retained and values below are removed. Default=0.001
#' @param kern.smooth scale at which kernel smoothing occurs on environmental data, larger values (i.e. 2) increase scale (making espace transitions smoother and typically larger) and smaller values (i.e. 0.5) decrease scale (making occupied espace clusters more dense and irregular). Default value is 1.  You can also input: "auto", which estimates the kernel parameter by calculating the standard deviation of rescaled PC1 and PC2 coordinates divided by the sixth root of the number of locations. This method can be unreliable when used on multimodal espace distributions as it results in over-smoothing of home ranges.  Multimodal espace occupancy can be somewhat common when a species occupies an extreme aspect of habitat or when espace is not broadly accessible in both dimensions of espace (PCs 1 & 2)
#' @param env.pca calibrate PCA on only occurrence data (env.pca=F) or environment data (env.pca=T), recommended is env.pca=T
#' @param e.reps the number of iterations for the equivalency test (humboldt.equivalency.test). Values higher than 200 are recommend for final analysis
#' @param b.reps the number of iterations for the Background tests (humboldt.background.test). Values higher than 200 are recommend for final analyses
#' @param p.overlap if p.overlap=T, the niche overlap plot will be output (humboldt.plot.overlap). Turn this to 'F' if you want to speed of analyses. This plot takes a lot of CPU time because it calculates kernel densities for all environment and species datasets input
#' @param p.boxplot if p.boxplot=T, a boxplot niche overlap plot will be created. The whisker of boxplot depict the environmental space of environment. Whereas the bars depict environmental space of each species within that environment. Dots on bars depict density of species localities in environmental space. This is a quick plot and doesn't requires lots of CPU time
#' @param p.overlap if p.scatter=T, this will plot scatter plots with histrograms of your environment and species datasets using the humboldt.plot.scatter function. This is a quick plot and doesn't requires lots of CPU time
#' @param run.silent if run.silent=T, texts boxes displaying 'sampling', 'rarefying', 'equivalency test', 'Background test' progress will not be displayed
#' @param nae.window the spatial window from which non-analogous environments will be quantified.  The non-analogous environments are characterized by gridding the espace of env1 and env2 into a R x R grind (e.g. 100 x 100).  If nae.window=0, values absent from a cell in one environment will be removed from the other.   If nae.window>0, values absent from a window (or neighborhood of cells) in one environment will be trimmed from the other.  The nae.window value characterizes the number of cells to search from the focal cell of environmental space values in the other environment.  The larger the nae.window value, the fewer non-analogous environments removed. This parameter allows imperfect overlap of environments.   If areas of environmental space are a little patchy between environments—but generally present-- a larger nae.window value will retain more of the patch environments. The default value is a nae.window=5
#' @param e.sample.method this parameter controls how espace is sampled.  There are two options: "points" or "density".  If e.sample.method="points", only input localities for both species are used for randomizations. If e.sample.method="density", the gridded density of species' espace is randomly resampled with higher sampling probabilities applied to areas of high density espace.  In both cases the number of samples selected for each species is equal to the observed datasets.  Again 'points' focuses only on input points, whereas 'density' samples the entire gridded density of each species' espace. Generally I prefer using 'points' for most datasets, because it maintains most of the nuances associated the original datasets (i.e. spatial autocorrelation).  Though sometimes I do use 'density' for datasets with limited localities (<20 points)
#' @param nae do you include non-analogous environments in the niche overlap measurement? If nae="NO" (use captial letters), then non-analogous environments will be removed from both input environments during overlap measurement and only environments present in both datasets will be used to measure overlap. If nae="YES", then no change will be made to input z1 and z2. Note: this is separate from trimming non-analogous environments from your input dataset (as done by humbodlt.g2 specified by parameter non.analogous.environments). This parameter physically removes non-analogous environments from datasets ONLY before the niche overlap measurement. Technically the removal of non-analogous environments via either way should result in similar overlap measurements (though they may not be identical). This because removing NAE from the dataset prior to gridding environments will resulting only non-analogous environments to be gridded (and typically finer grain applied to each grid cell). Whereas removing them only via this parameter (nae), which only removes non-analogous in the gridded environmental space for use in overlap measurements--- all the input environmental space is gridded (likely increasing the environmental space per gridded cell). A second cause of differences in values can result from rescaling of espace values during niche-overlap measurements so that the sum of the landscape equals one. If occupied non-analogous environmental are numerous in one of the datasets, this can theoretically cause overlap values to decrease in analogous environments (vs. nae) because differences in core niches are rescaled to 1 in both scenarios. The rescaling among fewer cells increases the values applied to highly suitable areas and, if not equivalently scaled in both datasets, differences among niches could increase, resulting a smaller overlap in non-analogous environments (again values should be similar). If you remove non-analogous environments in humboldt.g2e, I also suggest that you use this function (as it can remove any slight anomalies caused by gridding environments in humbodlt.grid.clim due to the binning of values in the RxR grid).
#' @return Performs almost all the analyses in Humboldt and outputs the results as Jason likes! See my example for my recommended running. For each comparison I recommend running the 'doitall' function twice. In the first run, I highly recommend running the analysis on the full environment input for both background and equivalency tests. In this case you are testing the total equivalence between species (or divergence) in their current realized distributions. It asks the question how equivalent are the two species realized niches?  This test is also the best gauge for background test and the power to discriminate differences.  I only use background values from this.  However, I do re-run the analyses with trimmed, shared espace. This is extremely important if you are interested in quantifying niche evolution/divergence among two populations/species . This test how two species diverge in shared analogous e-space and is the only actual test of niche divergence/evolution. Note that non-significant background test are not uncommon in shared-analogous espace and should not be given priority overlap total espace measurements from run 1. Under the situation where there is very little or no analogous e-space ---you may not be able conclude much from this situation regarding niche divergence. In this situation, if the species' 'niches' are not truncated by espace boundaries, then this suggests that the g2e conversion is a decent approximation of the species' fundamental niches. And because of this, the results should be accurate and species divergence can also be concluded. However if the core densities on one or both of species are on margins or near margins of espace, then this suggests that the fundamental niche may be much larger than available climate space provides. In this latter case, you cannot conclude that the species have diverged at all. In the former case, where there is no truncation of espace of each species in their habitats, but little overlap between, this is evidence the species have diverged.  In this case, simply report the niche overlap values in the total environmental space and the amount analogous espace.
#' @seealso \code{humboldt.g2e, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche, humboldt.doitall, humboldt.top.env} which use or depend on outputs of this function 
#' @import raster
#' @export
#' @examples
#' #######################################################################################
#' ###################################    EXAMPLE 1    ###################################
#' #######################################################################################
#' library(humboldt)
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)

#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be sp,x,y
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#'
#' ##its highly recommened that you using the function "humboldt.top.env" to select only the important enviromnetal variables in humbodlt.doitall. This step can be skipped. If you downloaded tons of environmental data, you should use this step.  If you skip this step, input env1/env2 inplace of reduc.vars$env1/reduc.vars$env2 
#' reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=occ.sp1,sp2=occ.sp2,rarefy.dist=40, rarefy.units="km", env.reso=0.416669,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)
#' 
#' ##Adjust the number of variables input for e.vars after reduction to only important variables
#' num.var.e<-ncol(reduc.vars$env1)

#' ##run it first with full environmental for backgroud tests and equivalency test (total equivalence or divergence in current distributions)
#' full<-humboldt.doitall(inname="full_extent", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=occ.sp1, sp2=occ.sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=0, reductype="PCA", non.analogous.environments="YES", correct.env=T, env.trim.mcp=F, mcp.buffer.sp1=200, mcp.buffer.sp2=200, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, PROJ=F, e.reps=100, e.sample.method="points", b.reps=100, env.pca=T, nae="YES",thresh.espace.z=0.001, p.overlap=T, p.boxplot=F, p.scatter=F, run.silent=F)
#'
#' ##run it a second time with a trimmed, shared-espace. Here the equivalency test tests for niche evolution or niche divergence. For comparing results, change only the following model parameters: reduce.env, non.analogous.environmental,env.trim.mcp
#' shared_ae<-humboldt.doitall(inname="shared_espace_ae", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=occ.sp1, sp2=occ.sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=2, reductype="PCA", non.analogous.environments="NO", correct.env=T, env.trim.mcp=T, mcp.buffer.sp1=200, mcp.buffer.sp2=200, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, PROJ=F, e.reps=100, e.sample.method="points", b.reps=100, env.pca=T, nae="YES",thresh.espace.z=0.001, p.overlap=T, p.boxplot=F, p.scatter=T,run.silent=F)

#' #######################################################################################
#' ###################################    EXAMPLE 2    ###################################
#' #######################################################################################
#' ############################  Using Provided Example Data   ###########################
#' #######################################################################################
#' library(humboldt)
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' data(env1)
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' data(env2)
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be sp,x,y
#' data(sp1)
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
#' data(sp2)
#'
#' ##its highly recommened that you using the function "humboldt.top.env" to select only the important enviromnetal variables in humbodlt.doitall. This step can be skipped. If you downloaded tons of environmental data, you should use this step.  If you skip this step, input env1/env2 inplace of reduc.vars$env1/reduc.vars$env2 
#' reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=sp1,sp2=sp2,rarefy.dist=50, rarefy.units="km", env.reso=0.416669,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)
#' 
#' ##Adjust the number of variables input for e.vars after reduction to only important variables
#' num.var.e<-ncol(reduc.vars$env1)

#' ##run it first with full environmental for backgroud tests and equivalency test (total equivalence or divergence in current distributions)
#' full<-humboldt.doitall(inname="full_extent", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=0, reductype="PCA", non.analogous.environments="YES", correct.env=T, env.trim.mcp=F, mcp.buffer.sp1=200, mcp.buffer.sp2=200, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, PROJ=F, e.reps=100, e.sample.method="points", b.reps=100, env.pca=T, nae="YES",thresh.espace.z=0.001, p.overlap=T, p.boxplot=F, p.scatter=F, run.silent=F)
#'
#' ##run it a second time with a trimmed, shared-espace. Here the equivalency test tests for niche evolution or niche divergence. For comparing results, change only the following model parameters: reduce.env, non.analogous.environmental,env.trim.mcp
#' shared_ae<-humboldt.doitall(inname="shared_espace_ae", env1=reduc.vars$env1, env2=reduc.vars$env2, sp1=sp1, sp2=sp2, rarefy.dist=50, rarefy.units="km", env.reso=0.416669, reduce.env=2, reductype="PCA", non.analogous.environments="NO", correct.env=T, env.trim.mcp=T, mcp.buffer.sp1=200, mcp.buffer.sp2=200, col.env=e.var, e.var=c(3:num.var.e), R=100, kern.smooth=1, PROJ=F, e.reps=100, e.sample.method="points", b.reps=100, env.pca=T, nae="YES",thresh.espace.z=0.001, p.overlap=T, p.boxplot=F, p.scatter=T,run.silent=F)


humboldt.doitall <- function(inname = "DoItAll", env1, env2, sp1, sp2, rarefy.dist = 0, rarefy.units = "dd", env.reso, 
    reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", nae.window=5, env.trim.mcp = T, mcp.buffer.sp1 = 200, 
    mcp.buffer.sp2 = 200, correct.env = T, col.env = e.var, e.var, R = 100, kern.smooth = 1, 
    PROJ = F, e.reps = 100, e.sample.method = "points", b.reps = 100, b.force.equal.sample=F, env.pca = T, nae = "YES", 
    thresh.espace.z = 0.001, p.overlap = T, p.boxplot = F, p.scatter = F, run.silent = F) {
    
    l <- list()
    
    #################################################################################################
    #################################### Gspace to Espace ###########################################
    #################################################################################################
    zz <<- humboldt.g2e(env1, env2, sp1, sp2, reduce.env, reductype, non.analogous.environments, nae.window,
        env.trim.mcp, e.var, col.env = e.var, mcp.buffer.sp1, mcp.buffer.sp2, env.pca, rarefy.dist,
        rarefy.units, env.reso, kern.smooth, PROJ, R, run.silent = run.silent)
    ##################################################################################################
    ############################# plots of PCA and Espace ############################################
    ############################# ####################################################################
    ptlss1 <- nrow(zz$scores.sp1)/nrow(zz$scores.sp1FULL)
    ptlss2 <- nrow(zz$scores.sp2)/nrow(zz$scores.sp2FULL)
    if (reduce.env != 0 & ptlss1 > 0.5 & nrow(zz$scores.sp2) < 10) {
        print(paste("Only ", nrow(zz$scores.sp2), " points remain in species dataset 2. In this analysis, it is recommended to not trim shared e-space or remove non-analogous climates because input habitats are too divergent.  Basically input environments are non-equivalent and have little e-space overlap relevant to species 2. Because e-space cannot be reduced to shared/non-analogous environments: report shared analogos space percentage, niche overlap values and equivalency test significance for full climate space (see below for correct.env parameter recommendations)"))
    }
    if (reduce.env != 0 & ptlss2 > 0.5 & nrow(zz$scores.sp1) < 10) {
        print(paste("Only ", nrow(zfz$scores.sp1), " points remain in species dataset 1. In this analysis, it is recommended to not trim shared e-space or remove non-analogous climates because input habitats are too divergent.  Basically input environments are non-equivalent and have little e-space overlap relevant to species 1. Because e-space cannot be reduced to shared/non-analogous environments: report shared analogos space percentage, niche overlap values and equivalency test significance for full climate space (see below for correct.env parameter recommendations)"))
    }
    
    threshinZ <- thresh.espace.z
    nacinZ <- nae
    kern.smoothinZ <- kern.smooth
    Rin <- R
    scores.env12 <- rbind(zz$scores.env1, zz$scores.env2)
    scores.env1a <- zz$scores.env1
    scores.env2a <- zz$scores.env2
    scores.env12a <- scores.env12
    scores.sp1a <- zz$scores.sp1
    scores.sp2a <- zz$scores.sp2
    scores.env1a$Axis2 <- (scores.env1a$Axis2 * -1)
    scores.env2a$Axis2 <- (scores.env2a$Axis2 * -1)
    scores.env12a$Axis2 <- (scores.env12a$Axis2 * -1)
    scores.sp1a$Axis2 <- (scores.sp1a$Axis2 * -1)
    scores.sp2a$Axis2 <- (scores.sp2a$Axis2 * -1)
    # calculation of occurrence density and test of niche equivalency and Background
    z1 <- humboldt.grid.espace(scores.env12a[1:2], scores.env1a[1:2], scores.env1a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    z2 <- humboldt.grid.espace(scores.env12a[1:2], scores.env2a[1:2], scores.env2a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    z3 <- humboldt.grid.espace(scores.env1a[1:2], scores.env1a[1:2], scores.env1a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    z4 <- humboldt.grid.espace(scores.env2a[1:2], scores.env2a[1:2], scores.env2a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    z5 <- humboldt.grid.espace(scores.env12a[1:2], scores.env1a[1:2], scores.sp1a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)
    z6 <- humboldt.grid.espace(scores.env12a[1:2], scores.env2a[1:2], scores.sp2a[1:2], kern.smooth = kern.smoothinZ, 
        R = Rin)   
    ## cacluate environment overlaps for adjusted Dvalues
    DvalClim <- round(as.numeric(humboldt.niche.overlap(z3, z4, correct.env = F, nae = "YES", 
        thresh.espace.z = threshinZ)[1]), 3)
    pdf(file = paste(inname, "_ENV.pdf", sep = ""))  # create a pdf file named from the names of the 2 species
    layout(matrix(c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4), 4, 4, byrow = TRUE))
    humboldt.plot.niche(z1, correct.env = F, title = "E-space Environment 1", name.axis1 = "PC1", 
        name.axis2 = "PC2")
    humboldt.plot.niche(z2, correct.env = F, title = "E-space Environment 2", name.axis1 = "PC1", 
        name.axis2 = "PC2")
    ## measure differences
    ee <- humboldt.espace.correction(Z.env1 = z1, Z.env2 = z2, Z.sp1 = z5, Z.sp2 = z6)
    ## stop if differences are highly correlated if(ee$D.overlap.uncorrected>0.5 &
    ## correct.env==F){stop('Analysis stopped. In this situation, it is best to change
    ## 'correct.env=F' to 'correct.env=T'. Please change this and rerun to proceed.')}
    if (correct.env == T) {
        print("Just so you know, you input: 'correct.env=T'")
    }
    # plot Differences
    
    humboldt.plot.espace.diff(espace.diff = ee, correct.env = F, type = "environment")
    if (ee$e.uncor.sum != 0) {
        contour(z1$x, (sort((z1$y * -1))), z1$Z, add = T, levels = quantile(z1$Z[z1$Z > 0], c(0.1, 
            0.5, 0.75)), drawlabels = F, lty = c(1, 2, 3), lwd = c(1, 1, 1), col = "grey")
    }
    humboldt.plot.contrib(zz$pca.cal$co, zz$pca.cal$eig)
    
    dev.off()
    pdf(file = paste(inname, "_ENV2.pdf", sep = ""))  # create a pdf file named from the names of the 2 species
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
    AmaxF <- max(zz$scores.env1FULL[, 1])
    AminF <- min(zz$scores.env1FULL[, 1])
    BmaxF <- max(zz$scores.env2FULL[, 1])
    BminF <- min(zz$scores.env2FULL[, 1])
    CmaxF <- max(zz$scores.env1FULL[, 2])
    CminF <- min(zz$scores.env1FULL[, 2])
    DmaxF <- max(zz$scores.env2FULL[, 2])
    DminF <- min(zz$scores.env2FULL[, 2])
    Xmin1 = pmin(AminF, BminF)
    Xmax1 = pmax(AmaxF, BmaxF)
    Ymin1 = pmin(CminF, DminF)
    Ymax1 = pmax(CmaxF, DmaxF)
    
    AmaxFG <- max(zz$scores.env1FULL[, 3])
    AminFG <- min(zz$scores.env1FULL[, 3])
    BmaxFG <- max(zz$scores.env2FULL[, 3])
    BminFG <- min(zz$scores.env2FULL[, 3])
    CmaxFG <- max(zz$scores.env1FULL[, 4])
    CminFG <- min(zz$scores.env1FULL[, 4])
    DmaxFG <- max(zz$scores.env2FULL[, 4])
    DminFG <- min(zz$scores.env2FULL[, 4])
    
    ### plot scatter plots of environment before and after
    plot(zz$scores.env1FULL[1:2], col = "light grey", xlim = c(Xmin1, Xmax1), ylim = c(Ymin1, 
        Ymax1), main = "E-space Environment 1", xlab = "PC1", ylab = "PC2", pch = 20, cex.sub = 0.75, 
        cex = 0.3, sub = (paste("Removed:", zz$RedP1C, "sites of", nrow(zz$scores.env1FULL[1:2]), 
            "(", round((100 * zz$RedP1C/nrow(zz$scores.env1FULL[1:2])), 1), "%);", zz$RedP1sp, 
            "locs of", nrow(zz$scores.sp1FULL[1:2]), "(", round((100 * zz$RedP1sp/nrow(zz$scores.sp1FULL[1:2])), 
                1), "%)")))
    points(zz$scores.env1[1:2], pch = 19, cex = 0.3)
    points(zz$scores.sp1FULL[1:2], pch = 19, col = "blue", cex = 0.3)
    points(zz$scores.sp1[1:2], pch = 19, col = "red", cex = 0.3)
    rect(min(zz$scores.env2[, 1], zz$scores.env1[, 1]), min(zz$scores.env2[, 2], zz$scores.env1[, 
        2]), max(zz$scores.env2[, 1], zz$scores.env1[, 1]), max(zz$scores.env2[, 2], zz$scores.env1[, 
        2]), density = NULL, angle = 45, col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"))
    plot(zz$scores.env2FULL[1:2], col = "light grey", xlim = c(Xmin1, Xmax1), ylim = c(Ymin1, 
        Ymax1), main = "E-space Environment 2", xlab = "PC1", ylab = "PC2", pch = 20, cex.sub = 0.75, 
        cex = 0.3, sub = (paste("Removed:", zz$RedP2C, "sites of", nrow(zz$scores.env2FULL[1:2]), 
            "(", round((100 * zz$RedP2C/nrow(zz$scores.env2FULL[1:2])), 1), "%);", zz$RedP2sp, 
            "locs of", nrow(zz$scores.sp2FULL[1:2]), "(", round((100 * zz$RedP2sp/nrow(zz$scores.sp2FULL[1:2])), 
                1), "%)")))
    points(zz$scores.env2[1:2], pch = 19, cex = 0.3)
    points(zz$scores.sp2FULL[1:2], pch = 19, col = "blue", cex = 0.3)
    points(zz$scores.sp2[1:2], pch = 19, col = "red", cex = 0.3)
    rect(min(zz$scores.env2[, 1], zz$scores.env1[, 1]), min(zz$scores.env2[, 2], zz$scores.env1[, 
        2]), max(zz$scores.env2[, 1], zz$scores.env1[, 1]), max(zz$scores.env2[, 2], zz$scores.env1[, 
        2]), density = NULL, angle = 45, col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"))
    
    plot(zz$scores.env1FULL[3:4], col = "light grey", xlim = c(AminFG, AmaxFG), ylim = c(CminFG, 
        CmaxFG), main = "G-space Environment 1", xlab = "Longitude", ylab = "Latitude", pch = 20, 
        cex.sub = 0.75, cex = 0.3, sub = (paste("Removed:", zz$RedP1C, "sites of", nrow(zz$scores.env1FULL[3:4]), 
            "(", round((100 * zz$RedP1C/nrow(zz$scores.env1FULL[3:4])), 1), "%);", zz$RedP1sp, 
            "locs of", nrow(zz$scores.sp1FULL[3:4]), "(", round((100 * zz$RedP1sp/nrow(zz$scores.sp1FULL[3:4])), 
                1), "%)")))
    points(zz$scores.env1[3:4], pch = 19, cex = 0.3)
    points(zz$scores.sp1FULL[3:4], pch = 19, col = "blue", cex = 0.3)
    points(zz$scores.sp1[3:4], pch = 19, col = "red", cex = 0.3)
    plot(zz$scores.env2FULL[3:4], col = "light grey", xlim = c(BminFG, BmaxFG), ylim = c(DminFG, 
        DmaxFG), main = "G-space Environment 2", xlab = "Longitude", ylab = "Latitude", pch = 20, 
        cex.sub = 0.75, cex = 0.3, sub = (paste("Removed:", zz$RedP2C, "sites of", nrow(zz$scores.env2FULL[3:4]), 
            "(", round((100 * zz$RedP2C/nrow(zz$scores.env2FULL[3:4])), 1), "%);", zz$RedP2sp, 
            "locs of", nrow(zz$scores.sp2FULL[3:4]), "(", round((100 * zz$RedP2sp/nrow(zz$scores.sp2FULL[3:4])), 
                1), "%)")))
    points(zz$scores.env2[3:4], pch = 19, cex = 0.3)
    points(zz$scores.sp2FULL[3:4], pch = 19, col = "blue", cex = 0.3)
    points(zz$scores.sp2[3:4], pch = 19, col = "red", cex = 0.3)
    
    dev.off()
    ### plot of espace on geospace
    pdf(file = paste(inname, "_ENV3.pdf", sep = ""))  # create a pdf file named from the names of the 2 species
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
    par(mar = c(2.5, 4, 2.5, 2.5))
    options(warn = -1)
    n.env1 <- 1:nrow(zz$scores.env1)
    n.env2 <- (nrow(zz$scores.env1) + 1):(nrow(zz$scores.env1) + nrow(zz$scores.env2))
    env12 <- rbind(zz$scores.env1, zz$scores.env2)
    Z1 <- env12[, 1]
    Z2 <- env12[, 2]
    X <- env12[, 3]
    Y <- env12[, 4]
    df <- data.frame(X, Y, Z1, Z2)
    rangeMin <- min(df[, 3])
    rangeMax <- max(df[, 3])
    rangeEnv <- (rangeMax - rangeMin)/256
    rangeEnvL <- (rangeMax - rangeMin)
    rangeBreaks <- seq(rangeMin, rangeMax, rangeEnv)
    r.range <- c(rangeMin, rangeMax)
    div.r.range <- rangeEnvL/5
    df <- df[n.env1, ]
    r3 <- raster(res = env.reso, xmn = min(df[, 1]), xmx = max(df[, 1]), ymn = min(df[, 2]), 
        ymx = max(df[, 2]))
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
    cells <- cellFromXY(r3, df[, 1:2])
    r3[cells] <- df[, 3]
    plot(r3, col = cols, cex.axis = 0.6, legend = FALSE, breaks = rangeBreaks)
    plot(r3, legend.only = TRUE, col = cols, legend.width = 1, legend.shrink = 0.5, axis.args = list(at = seq(r.range[1], 
        r.range[2], div.r.range), labels = round(seq(r.range[1], r.range[2], div.r.range), 1), 
        cex.axis = 0.6), legend.args = list(text = "PC1", side = 4, font = 2, line = 2.5, cex = 0.8))
    
    df <- data.frame(X, Y, Z1, Z2)
    df <- df[n.env2, ]
    r4 <- raster(res = env.reso, xmn = min(df[, 1]), xmx = max(df[, 1]), ymn = min(df[, 2]), 
        ymx = max(df[, 2]))
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
    cells <- cellFromXY(r4, df[, 1:2])
    r4[cells] <- df[, 3]
    plot(r4, col = cols, cex.axis = 0.6, legend = FALSE, breaks = rangeBreaks)
    
    df <- data.frame(X, Y, Z1, Z2)
    rangeMin <- min(df[, 4])
    rangeMax <- max(df[, 4])
    rangeEnv <- (rangeMax - rangeMin)/256
    rangeEnvL <- (rangeMax - rangeMin)
    rangeBreaks <- seq(rangeMin, rangeMax, rangeEnv)
    r.range <- c(rangeMin, rangeMax)
    div.r.range <- rangeEnvL/5
    df <- df[n.env1, ]
    r1 <- raster(res = env.reso, xmn = min(df[, 1]), xmx = max(df[, 1]), ymn = min(df[, 2]), 
        ymx = max(df[, 2]))
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
    cells <- cellFromXY(r1, df[, 1:2])
    r1[cells] <- df[, 4]
    plot(r1, col = cols, cex.axis = 0.6, legend = FALSE, breaks = rangeBreaks)
    plot(r1, legend.only = TRUE, col = cols, legend.width = 1, legend.shrink = 0.5, axis.args = list(at = seq(r.range[1], 
        r.range[2], div.r.range), labels = round(seq(r.range[1], r.range[2], div.r.range), 1), 
        cex.axis = 0.6), legend.args = list(text = "PC2", side = 4, font = 2, line = 2.5, cex = 0.8))
    
    df <- data.frame(X, Y, Z1, Z2)
    df <- df[n.env2, ]
    r2 <- raster(res = env.reso, xmn = min(df[, 1]), xmx = max(df[, 1]), ymn = min(df[, 2]), 
        ymx = max(df[, 2]))
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
    cells <- cellFromXY(r2, df[, 1:2])
    r2[cells] <- df[, 4]
    plot(r2, col = cols, cex.axis = 0.6, legend = FALSE, breaks = rangeBreaks)
    mtext("E-space plotted in G-space", side = 3, line = -2, outer = TRUE)
    options(warn = 0)
    dev.off()
    
    ###############################################################
    # calculation of occurrence density and test of niche equivalency and background
    a <- humboldt.equivalency.test(z5, z6, rep = e.reps, correct.env = correct.env, nae = nacinZ, 
        kern.smooth = kern.smoothinZ, sample = e.sample.method, run.silent.equ = run.silent)  # test of niche equivalency and background according to Warren et al. 2008
    b <- humboldt.background.test(g2e = zz, correct.env = correct.env, env.reso = env.reso, sim.dir = 1, 
        rep = b.reps, kern.smooth = kern.smoothinZ, R = Rin, force.equal.sampl=b.force.equal.sample, run.silent.bak = run.silent)
    b2 <- humboldt.background.test(g2e = zz, correct.env = correct.env, env.reso = env.reso, 
        sim.dir = 2, rep = b.reps, kern.smooth = kern.smoothinZ, R = Rin, force.equal.sampl=b.force.equal.sample, run.silent.bak = run.silent)
    ################################################################ 
    pdf(file = paste(inname, "_SPP.pdf", sep = ""))  # create a pdf file named from the names of the 2 species
    layout(matrix(c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 5, 3, 3, 6, 7), 4, 4, byrow = TRUE))
    humboldt.plot.niche(z5, title = "E-space Species 1", name.axis1 = "PC1", name.axis2 = "PC2", correct.env = F)
    humboldt.plot.niche(z6, title = "E-space Species 2", name.axis1 ="PC1", name.axis2 = "PC2", 
        correct.env = F)
    # humboldt.plot.contrib(zz$pca.cal$co,zz$pca.cal$eig) plot difference
    humboldt.plot.espace.diff(espace.diff = ee, correct.env = F, type = "species")
    
    if (ee$s.uncor.sum != 0) {
        contour(z5$x, (sort((z5$y * -1))), z5$Z, add = T, levels = quantile(z5$Z[z5$Z > 0], c(0.1, 
            0.5, 0.75)), drawlabels = F, lty = c(1, 2, 3), lwd = c(1, 1, 1), col = "grey")
    }
    if (ee$s.uncor.sum == 0) {
        contour(z5$x, (sort((z5$y * -1))), z5$Z, add = T, levels = quantile(z5$Z[z5$Z > 0], c(0.1, 
            0.5, 0.75)), drawlabels = F, lty = c(1, 2, 3), lwd = c(1, 1, 1), col = "white")
    }
    # other plots
    DvalSp <- round(as.numeric(humboldt.niche.overlap(z5, z6, correct.env = correct.env, nae = "YES", 
        thresh.espace.z = threshinZ)[1]), 3)
    DvalSp2 <- round(as.numeric(humboldt.niche.overlap(z5, z6, correct.env = correct.env, nae = "NO", 
        thresh.espace.z = threshinZ)[1]), 3)
    if (DvalSp2 < DvalSp) {
        DvalSp2 <- DvalSp
    }
    DvalAdjSp <- round((DvalSp/DvalClim), 3)
    DvalAdjSp[DvalAdjSp > 1] <- 1
    print("Overlap:")
    print(DvalSp)
    # print('Overlap - adjusted:') print(DvalAdjSp)
    print("Overlap - analogous environments only:")
    print(DvalSp2)
    plot.new()
    text(0.5, 0.5, paste("Niche Overlap:", "\n", "D=", DvalSp, "\n", "Analog Env ONLY:", "\n", 
        "D=", DvalSp2), cex = 1)
    humboldt.plot.histrogram(a, "D", "Equivalency")
    humboldt.plot.density(b, "D", "Background 2->1")
    humboldt.plot.density(b2, "D", "Background 1->2")
    dev.off()
    
    print("FINISHED ANALYSES: now generating supplemental plots (if selected)")
    ##################################################################################################
    ############################# scatter plots and hist of datasets #################################
    ##################################################################################################
    if (p.scatter == T) {
        pdf(file = paste(inname, "_SCATTER.pdf", sep = ""))  # create a pdf file named from the names of the 2 species
        humboldt.plot.scatter(zz$scores.env1[1:2], xlab = "PC1", ylab = "PC2", main = "Environment 1 PCA")
        humboldt.plot.scatter(zz$scores.env2[1:2], xlab = "PC1", ylab = "PC2", main = "Environment 2 PCA")
        humboldt.plot.scatter(zz$scores.sp1[1:2], xlab = "PC1", ylab = "PC2", main = "Sp/pop 1 PCA")
        humboldt.plot.scatter(zz$scores.sp2[1:2], xlab = "PC1", ylab = "PC2", main = "Sp/pop 2 PCA")
        
        dev.off()
    }
    ##################################################################################################
    ############################## box and line plots of old/new datasets ############################
    ##################################################################################################

    if (p.boxplot == T) {
        AmaxPCA1 <- max(zz$scores.env1[, 1])
        AminPCA1 <- min(zz$scores.env1[, 1])
        AmeanPCA1 <- mean(zz$scores.env1[, 1])
        BmaxPCA1 <- max(zz$scores.env2[, 1])
        BminPCA1 <- min(zz$scores.env2[, 1])
        BmeanPCA1 <- mean(zz$scores.env2[, 1])
        XminPCA1 = pmin(AmaxPCA1, BmaxPCA1)
        XmaxPCA1 = pmax(AminPCA1, BminPCA1)
        AmaxPCA2 <- max(zz$scores.env1[, 2])
        AminPCA2 <- min(zz$scores.env1[, 2])
        AmeanPCA2 <- mean(zz$scores.env1[, 2])
        BmaxPCA2 <- max(zz$scores.env2[, 2])
        BminPCA2 <- min(zz$scores.env2[, 2])
        BmeanPCA2 <- mean(zz$scores.env2[, 2])
        XminPCA2 = pmin(AmaxPCA2, BmaxPCA2)
        XmaxPCA2 = pmax(AminPCA2, BminPCA2)
        CmaxPCA1 <- max(zz$scores.sp1[, 1])
        CminPCA1 <- min(zz$scores.sp1[, 1])
        CmeanPCA1 <- mean(zz$scores.sp1[, 1])
        DmaxPCA1 <- max(zz$scores.sp2[, 1])
        DminPCA1 <- min(zz$scores.sp2[, 1])
        DmeanPCA1 <- mean(zz$scores.sp2[, 1])
        CmaxPCA2 <- max(zz$scores.sp1[, 2])
        CminPCA2 <- min(zz$scores.sp1[, 2])
        CmeanPCA2 <- mean(zz$scores.sp1[, 2])
        DmaxPCA2 <- max(zz$scores.sp2[, 2])
        DminPCA2 <- min(zz$scores.sp2[, 2])
        DmeanPCA2 <- mean(zz$scores.sp2[, 2])
        EmaxPCA1 <- max(zz$scores.sp1FULL[, 1])
        EminPCA1 <- min(zz$scores.sp1FULL[, 1])
        EmeanPCA1 <- mean(zz$scores.sp1FULL[, 1])
        FmaxPCA1 <- max(zz$scores.sp2FULL[, 1])
        FminPCA1 <- min(zz$scores.sp2FULL[, 1])
        FmeanPCA1 <- mean(zz$scores.sp2FULL[, 1])
        EmaxPCA2 <- max(zz$scores.sp1FULL[, 2])
        EminPCA2 <- min(zz$scores.sp1FULL[, 2])
        EmeanPCA2 <- mean(zz$scores.sp1FULL[, 2])
        FmaxPCA2 <- max(zz$scores.sp2FULL[, 2])
        FminPCA2 <- min(zz$scores.sp2FULL[, 2])
        FmeanPCA2 <- mean(zz$scores.sp2FULL[, 2])
        
        
        
        pdf(file = paste(inname, "_BOXPLOT.pdf", sep = ""))
        layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))  # create a pdf file named from the names
        X1 = c(AminF, EminPCA1, EmeanPCA1, EmaxPCA1, AmaxF)
        Y1 = c(BminF, FminPCA1, FmeanPCA1, FmaxPCA1, BmaxF)
        XX1 = c(AminPCA1, CminPCA1, CmeanPCA1, CmaxPCA1, AmaxPCA1)
        YY1 = c(BminPCA1, DminPCA1, DmeanPCA1, DmaxPCA1, BmaxPCA1)
        df = data.frame(X1, Y1, XX1, YY1)
        
        
        bp = boxplot(df, horizontal = TRUE, axes = TRUE, staplewex = 1, lwd = 1, whisklty = 1, 
            medlwd = 0, medcol = c("olivedrab3", "skyblue2", "olivedrab3", "skyblue2"), col = c("olivedrab3", 
                "skyblue2", "olivedrab3", "skyblue2"), range = 0, xlab = "environmental space- PC1", 
            at = c(5, 4, 2, 1), names = c("sp/pop1", "sp/pop2", "overlap\nsp/pop1", "overlap\nsp/pop2"), 
            las = 1)
        abline(v = AmaxPCA1, lty = 2)  #XmaxPCA1
        abline(v = AminPCA1, lty = 2)
        aZ = data.frame(zz$scores.sp1FULL[, 1], "sp1F")
        bZ = data.frame(zz$scores.sp2FULL[, 1], "sp2F")
        aF = data.frame(zz$scores.sp1[, 1], "sp1")
        bF = data.frame(zz$scores.sp2[, 1], "sp2")
        
        names(aF)[1] <- "Val"
        names(aZ)[1] <- "Val"
        names(bF)[1] <- "Val"
        names(bZ)[1] <- "Val"
        names(aF)[2] <- "G"
        names(aZ)[2] <- "G"
        names(bF)[2] <- "G"
        names(bZ)[2] <- "G"
        df2 = rbind(aZ, bZ, aF, bF)
        stripchart(Val ~ G, vertical = FALSE, data = df2, method = "jitter", jitter = 0.3, add = TRUE, 
            pch = 20, col = "black", at = c(5, 4, 2, 1))
        
        X1 = c(CminF, EminPCA2, EmeanPCA2, EmaxPCA2, CmaxF)
        Y1 = c(DminF, FminPCA2, FmeanPCA2, FmaxPCA2, DmaxF)
        XX1 = c(AminPCA2, CminPCA2, CmeanPCA2, CmaxPCA2, AmaxPCA2)
        YY1 = c(BminPCA2, DminPCA2, DmeanPCA2, DmaxPCA2, BmaxPCA2)
        df3 = data.frame(X1, Y1, XX1, YY1)
        
        bp2 = boxplot(df3, horizontal = TRUE, axes = TRUE, staplewex = 1, lwd = 1, whisklty = 1, 
            medlwd = 0, medcol = c("olivedrab3", "skyblue2", "olivedrab3", "skyblue2"), col = c("olivedrab3", 
                "skyblue2", "olivedrab3", "skyblue2"), range = 0, xlab = "environmental space- PC2", 
            at = c(5, 4, 2, 1), names = c("sp/pop1", "sp/pop2", "overlap\nsp/pop1", "overlap\nsp/pop2"), 
            las = 1)
        abline(v = AminPCA2, lty = 2)
        abline(v = AmaxPCA2, lty = 2)
        aZ = data.frame(zz$scores.sp1FULL[, 2], "sp1F")
        bZ = data.frame(zz$scores.sp2FULL[, 2], "sp2F")
        aF = data.frame(zz$scores.sp1[, 2], "sp1")
        bF = data.frame(zz$scores.sp2[, 2], "sp2")
        
        names(aF)[1] <- "Val"
        names(aZ)[1] <- "Val"
        names(bF)[1] <- "Val"
        names(bZ)[1] <- "Val"
        names(aF)[2] <- "G"
        names(aZ)[2] <- "G"
        names(bF)[2] <- "G"
        names(bZ)[2] <- "G"
        df4 = rbind(aZ, bZ, aF, bF)
        stripchart(Val ~ G, vertical = FALSE, data = df4, method = "jitter", jitter = 0.3, add = TRUE, 
            pch = 20, col = "black", at = c(5, 4, 2, 1))
        dev.off()
    }
    ##################################################################################################
    ####################################### overlapPlots ############################################
    #################################################################################################
    if (p.overlap == T) {
        data.env.occ2 <- rbind(zz$occ.sp1, zz$occ.sp2, zz$env1, zz$env2)[e.var]
        data.env.occ2 <- data.env.occ2[1:2]
        nSP1 <- nrow(zz$occ.sp1[1:2])
        nSP2 <- nrow(zz$occ.sp2[1:2])
        nCM1 <- nrow(zz$env1[1:2])
        nCM2 <- nrow(zz$env2[1:2])
        
        pca.cal2 <- dudi.pca(data.env.occ2, center = T, scale = T, scannf = F, nf = 2)
        # pca.cal2$li$Axis2<-(pca.cal2$li$Axis2*-1)
        pca.cal2$li$Axis1 <- (pca.cal2$li$Axis1 * -1)  #use
        humboldt.plot.overlap(in.pca = pca.cal2, nsp1 = nSP1, nsp2 = nSP2, nenv1 = nCM1, nenv2 = nCM2, 
            pdf.out = T, pdfname = paste(inname, "-SP-OVERLAP.pdf", sep = ""))
    }
l$g2e <- zz
l$eqiv <- a
l$b21 <- b
l$b12 <- b2
invisible(l)
}



##################################################################################################
#######################Select top environmental variables for PCA 
##################################################################################################
#' Select top environmental variables for PCA 
#' @param env1 environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice
#' @param env2 environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn; with X1-Xn being any string label. If env1=env2, input the same file twice
#' @param sp1 occurrence sites for the species/population 1 at study area 1 (env1). Column names should be 'sp', 'x','y'
#' @param sp2 occurrence sites for the species/population 2 at study area 2 (env2). Column names should be 'sp', 'x','y'
#' @param rarefy.dist removes occurrences within a minimum distance (specified here) to each other (this function uses the humboldt.occ.rarefy function). Values need to be in km[recommended] or decimal degrees. See associated parameter rarefy.units. Note: rarefy.dist=0 will remove no occurrences 
#' @param rarefy.units the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees
#' @param env.reso the resolution of the input environmental data grid in decimal degrees
#' @param e.var Selection of variables to include in evaluation for each species
#' @param learning.rt1 value from 0.01 to 0.001 for building SDM, start with 0.01 and if prompted, change to 0.001. the default value is 0.01
#' @param learning.rt2 value from 0.01 to 0.001 for building SDM, start with 0.01 and if prompted, change to 0.001. The default value is 0.01
#' @param steps1 numbers of trees to add at each cycle for modelling sp1. Start with 50 and if you run into problems gradually decrease, stopping at 1. The default value is 50
#' @param steps2 numbers of trees to add at each cycle for modelling sp2. Start with 50 and if you run into problems gradually decrease, stopping at 1. The default value is 50
#' @param nvars.save if method="nvars",this variable is required. It is the number of the top variables to save per species. The kept variables are selected by their relative influence in predicting the species distribution, selecting for the highest contributing variables. Often the total variables retained is lower due to identical variables select among both species. The default value is 5.  This value will be ignored if method="estimate" or "contrib" 
#' @param method this determines how important environmental variables are selected.  There are three options: "estimate", "contrib", "nvars". If method="estimate", the boosted regression tree algorithm will choose the number of variables to include by systematically removing variables until average change in the model exceeds the original standard error of deviance explained.  This is the most computationally intensive method. If method="contrib", variables above a relative influence value will be kept. See associated parameter 'contrib.greater'. If method="nvars", a fixed number of user specified variables will be kept. The kept variables are selected by their relative influence, selecting for the highest contributing variables. See associated parameter 'nvars.save'
#' @param contrib.greater if method="contrib", this variable is required. The kept variables are selected for their relative influence in predicting the species' distribution.  Here users select variables equal to or above an input model contribution value. The default value for this method is 5 (= variables with 5 percent or higher contribution to model of either species are kept). This value will be ignored if method="estimate" or "nvars"
#' @param pa.ratio ratio of pseudoabsences to occurrence points, typically this is 4. The null value is 4
#' @return This function runs generalized boosted regression models (a machine learning SDM algorithm) to select top parameters for inclusion in PCA. This is important because you want the PC to reflect variables that are relevant to the species distribution. Alternatively you can run Maxent outside of R and manually curate the variables you include (also recommended). 
#' @import dismo
#' @import sp
#' @import raster
#' @import gbm
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##perform modeling to determin imporant variables
#' reduc.vars<- humboldt.top.env(env1=env1,env2=env2,sp1=occ.sp1,sp2=occ.sp2,rarefy.dist=40, rarefy.units="km", env.reso=0.0833338,learning.rt1=0.01,learning.rt2=0.01,e.var=(3:21),pa.ratio=4,steps1=50,steps2=50,method="contrib",contrib.greater=5)
#' 
#' ##use new variables for env1 and evn2, use as you normally would do for env1/env2 (input above)
#' 
#' ##for example, input into converted geographic space to espace
#' ##zz<-humboldt.g2e(env1=reduc.vars$env1, env2=reduc.vars$env2....

humboldt.top.env <- function(env1, env2, sp1, sp2, rarefy.dist = 0, rarefy.units="km", env.reso, learning.rt1 = 0.01, 
    learning.rt2 = 0.01, e.var, pa.ratio = 4, steps1 = 50, steps2 = 50, method = "contrib", nvars.save = 5, 
    contrib.greater = 5) # method=='estimate' method=='contrib' method=='nvars'
{
    l <- list()
    x2var <- e.var + 1
    occ.sp1in <- humboldt.occ.rarefy(in.pts = sp1, colxy = 2:3, rarefy.dist = rarefy.dist, rarefy.units=rarefy.units)
    occ.sp2in <- humboldt.occ.rarefy(in.pts = sp2, colxy = 2:3, rarefy.dist = rarefy.dist, rarefy.units=rarefy.units)
    # create sp occurrence dataset by adding environmental variables from the global
    # environmental datasets resolution should be the resolution of the environmental data grid
    occ.sp1 <- na.exclude(humboldt.sample.spp(dfsp = occ.sp1in, colspxy = 2:3, colspkept = NULL, 
        dfvar = env1, colvarxy = 1:2, colvar = "all", resolution = env.reso))
    occ.sp2 <- na.exclude(humboldt.sample.spp(dfsp = occ.sp2in, colspxy = 2:3, colspkept = NULL, 
        dfvar = env2, colvarxy = 1:2, colvar = "all", resolution = env.reso))
    ## proper number of PA
    nnsp1 <- nrow(occ.sp1)
    nnenv1 <- nrow(env1)
    nnpa1 <- nnsp1 * pa.ratio
    if (nnpa1 > nnenv1) {
        nnpa1 <- nnenv1
    }
    nnsp2 <- nrow(occ.sp2)
    nnenv2 <- nrow(env2)
    nnpa2 <- nnsp2 * pa.ratio
    if (nnpa2 > nnenv2) {
        nnpa2 <- nnenv2
    }
    
    ## sp1
    climSample <- env1[sample(nrow(env1), nnpa1), ]
    nClim <- nrow(climSample)
    ID <- rep(0, nClim)
    clim <- cbind(ID, climSample)
    nSp <- nrow(occ.sp1)
    ID <- rep(1, nSp)
    sp <- cbind(ID, occ.sp1)
    datamodel <- rbind(sp, clim)
    
    if (method == "nvars" || method == "contrib") {
        modelOutA <- gbm.step(data = datamodel, gbm.x = x2var, gbm.y = 1, family = "bernoulli", 
            tree.complexity = 5, learning.rate = learning.rt1, step.size = steps1, bag.fraction = 0.5, 
            plot.main = FALSE, plot.folds = FALSE)
    }
    if (method == "estimate") {
        modelOut1A <<- gbm.step(data = datamodel, gbm.x = x2var, gbm.y = 1, family = "bernoulli", 
            tree.complexity = 5, learning.rate = learning.rt1, step.size = steps1, bag.fraction = 0.5, 
            plot.main = FALSE, plot.folds = FALSE)
        mod.simpA <<- gbm.simplify(modelOut1A)
        min.y <- min(c(0, mod.simpA$deviance.summary$mean))
        n.vars1 <- match(min.y, c(0, mod.simpA$deviance.summary$mean)) - 1
        vars1.use <- mod.simpA$pred.list[[n.vars1]]
        names(env1)[vars1.use]
        modelOutA <- gbm.step(data = datamodel, gbm.x = mod.simpA$pred.list[[n.vars1]], gbm.y = 1, 
            family = "bernoulli", tree.complexity = 5, learning.rate = learning.rt1, step.size = steps1, 
            bag.fraction = 0.5, plot.main = FALSE, plot.folds = FALSE)
    }
    # varInc1<-modelOut$contributions[1:5,][,1]
    if (method == "nvars") {
        varInc1 <- toString(modelOutA$contributions[1:nvars.save, ][, 1])
    }
    if (method == "estimate") {
        varInc1 <- toString(modelOutA$contributions[, 1])
    }
    if (method == "contrib") {
        varInc1 <- toString(modelOutA$contributions$var[modelOutA$contributions$rel.inf > contrib.greater])
    }
    ## sp2
    climSample <- env2[sample(nrow(env2), nnpa2), ]
    nClim <- nrow(climSample)
    ID <- rep(0, nClim)
    clim <- cbind(ID, climSample)
    nSp <- nrow(occ.sp2)
    ID <- rep(1, nSp)
    sp <- cbind(ID, occ.sp2)
    datamodel <- rbind(sp, clim)
    
    if (method == "nvars" || method == "contrib") {
        modelOutB <- gbm.step(data = datamodel, gbm.x = x2var, gbm.y = 1, family = "bernoulli", 
            tree.complexity = 5, learning.rate = learning.rt2, step.size = steps2, bag.fraction = 0.5, 
            plot.main = FALSE, plot.folds = FALSE)
    }
    
    if (method == "estimate") {
        modelOut1B <- gbm.step(data = datamodel, gbm.x = x2var, gbm.y = 1, family = "bernoulli", 
            tree.complexity = 5, learning.rate = learning.rt2, step.size = steps2, bag.fraction = 0.5, 
            plot.main = FALSE, plot.folds = FALSE)
        mod.simpB <<- gbm.simplify(modelOut1B)
        min.y <- min(c(0, mod.simpB$deviance.summary$mean))
        n.vars2 <- match(min.y, c(0, mod.simpB$deviance.summary$mean)) - 1
        vars2.use <- mod.simpB$pred.list[[n.vars2]]  #-1
        names(env1)[vars2.use]
        modelOutB <- gbm.step(data = datamodel, gbm.x = mod.simpB$pred.list[[n.vars2]], gbm.y = 1, 
            family = "bernoulli", tree.complexity = 5, learning.rate = learning.rt2, step.size = steps2, 
            bag.fraction = 0.5, plot.main = FALSE, plot.folds = FALSE)
    }
    if (method == "nvars") {
        varInc2 <- toString(modelOutB$contributions[1:nvars.save, ][, 1])
    }
    if (method == "estimate") {
        varInc2 <- toString(modelOutB$contributions[, 1])
    }
    if (method == "contrib") {
        varInc2 <- toString(modelOutB$contributions$var[modelOutB$contributions$rel.inf > contrib.greater])
    }
    
    # select columns in model
    varInc1a <- gsub(" ", "", varInc1, fixed = TRUE)
    varInc1b <- unlist(strsplit(varInc1a[1], split = ","))
    varInc2a <- gsub(" ", "", varInc2, fixed = TRUE)
    varInc2b <- unlist(strsplit(varInc2a[1], split = ","))
    print(varInc1a)
    print(varInc2a)
    
    # begin partition environment
    env1sp1 <- subset(env1, select = c(varInc1b))
    env1sp2 <- subset(env1, select = c(varInc2b))
    env2sp1 <- subset(env2, select = c(varInc1b))
    env2sp2 <- subset(env2, select = c(varInc2b))
    # reconstruct original environmental variables with xy
    env1pre <- cbind(env1[1:2], env1sp1, env1sp2)
    env2pre <- cbind(env2[1:2], env2sp1, env2sp2)
    # store new environment columns and remove duplicate names
    l$env1 <- env1pre[, !duplicated(colnames(env1pre))]
    l$env2 <- env2pre[, !duplicated(colnames(env2pre))]
    return(l)
    
}

############################################################################################
#################### Espace correction of environment abundance test ######################
#################################################################################################
#' Espace correction of environmental abundance test
#' @param Z.env1 a grid of espace for environment 1 (output from humboldt.grid.espace)
#' @param Z.env2 a grid of espace for environment 2 (output from humboldt.grid.espace)
#' @param Z.sp1 a grid of espace for species/population 1 (output from humboldt.grid.espace)
#' @param Z.sp2 a grid of espace for species/population 2 (output from humboldt.grid.espace)
#' @return Determines if users should correct occurrence densities of each species by the prevalence of the environments in their range for equivalency, background, and overlap analyses (correct.env=T). Often datasets have high overlap among the differences between input environments and the differences between species distributions in environmental space. If ignored, equivalency tests are prone to type I errors and you may observed statistical differences that are entirely due to differences in the availability of environments vs. actual differences in occupied environmental space. When highly correlated, it is strongly encouraged to use 'correct.env=T' to correct espace observations by abundance of environments.\cr
#' \cr
#' Output: an environmental dataset with only 'important' variables for inclusion into a PCA (vs. all variables)
#' @seealso \code{humboldt.g2e, humboldt.grid.espace, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace socres for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run create a grid of Environmental Space Function
#' Z.sp1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' Z.sp2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
#' Z.env1<- humboldt.grid.espace(scores.env12,scores.env1,scores.env1,kern.smooth=1,R=100)
#' Z.env2<- humboldt.grid.espace(scores.env12,scores.env2,scores.env2,kern.smooth=1,R=100)
#' 
#' ee<- humboldt.espace.correction(Z.env1=Z.env1,Z.env2=Z.env2,Z.sp1=Z.sp1,Z.sp2=Z.sp2)
#'

humboldt.espace.correction <- function(Z.env1, Z.env2, Z.sp1, Z.sp2) {
    
    l <- list()
    ## new code
    ZspO <- humboldt.niche.overlap(Z.sp1, Z.sp2, correct.env = F, nae = "NO", thresh.espace.z = 0.001)
    ZenvO <- humboldt.niche.overlap(Z.env1, Z.env2, correct.env = F, nae = "NO", thresh.espace.z = 0.001)
    # if (ZspO$remaining.espace.per<20 & ZspO$D*10<ZspO$remaining.espace.per &
    # ZenvO$D*10<ZspO$remaining.espace.per)
    
    ### UNCORRECTED MEASUREMENTS measure differences between habitat environmental space
    Zdiff.C.uncor <- Z.env1$z.uncor - Z.env2$z.uncor
    Zdiff.S.uncor <- Z.sp1$z.uncor - Z.sp2$z.uncor
    Zdiff.C.cor <- Z.env1$z.cor - Z.env2$z.cor
    Zdiff.S.cor <- Z.sp1$z.cor - Z.sp2$z.cor
    
    # Z4[Z4<max(Z4)/1000]<-0
    Zdiff.C.uncor.sum <- sum(Zdiff.C.uncor)
    Zdiff.C.uncor.max <- max(Zdiff.C.uncor)
    Zdiff.C.uncor.min <- abs(min(Zdiff.C.uncor))
    Zdiff.C.uncor.Z <- 1/max(Zdiff.C.uncor.max, Zdiff.C.uncor.min)
    if (Zdiff.C.uncor.Z == Inf) {
        Zdiff.C.uncor.Z = 0
    }
    if (Zdiff.C.uncor.Z != 0) {
        Zdiff.C.uncor <- Zdiff.C.uncor * Zdiff.C.uncor.Z
        Zdiff.C.uncor <- Zdiff.C.uncor + 1
        Zdiff.C.uncor <- Zdiff.C.uncor/2
    }
    if (Zdiff.C.uncor.Z == 0) {
        Zdiff.C.uncor <- Zdiff.C.uncor
    }
    
    
    ## measure differences between species environmental space Z4[Z4<max(Z4)/1000]<-0
    Zdiff.S.uncor.sum <- sum(Zdiff.S.uncor)
    Zdiff.S.uncor.max <- max(Zdiff.S.uncor)
    Zdiff.S.uncor.min <- abs(min(Zdiff.S.uncor))
    Zdiff.S.uncor.Z <- 1/max(Zdiff.S.uncor.max, Zdiff.S.uncor.min)
    if (Zdiff.S.uncor.Z == Inf) {
        Zdiff.C.uncor.S = 0
    }
    if (Zdiff.S.uncor.Z != 0) {
        Zdiff.S.uncor <- Zdiff.S.uncor * Zdiff.S.uncor.Z
        Zdiff.S.uncor <- Zdiff.S.uncor + 1
        Zdiff.S.uncor <- Zdiff.S.uncor/2
    }
    if (Zdiff.S.uncor.Z == 0) {
        Zdiff.S.uncor <- Zdiff.S.uncor
    }
    
    ### CORRECTED MEASUREMENTS Z4[Z4<max(Z4)/1000]<-0
    Zdiff.C.cor.sum <- sum(Zdiff.C.cor)
    Zdiff.C.cor.max <- max(Zdiff.C.cor)
    Zdiff.C.cor.min <- abs(min(Zdiff.C.cor))
    Zdiff.C.cor.Z <- 1/max(Zdiff.C.cor.max, Zdiff.C.cor.min)
    if (Zdiff.C.cor.Z == Inf) {
        Zdiff.C.cor.Z = 0
    }
    if (Zdiff.C.cor.Z != 0) {
        Zdiff.C.cor <- Zdiff.C.cor * Zdiff.C.cor.Z
        Zdiff.C.cor <- Zdiff.C.cor + 1
        Zdiff.C.cor <- Zdiff.C.cor/2
    }
    if (Zdiff.C.cor.Z == 0) {
        Zdiff.C.cor <- Zdiff.C.cor
    }
    
    ## measure differences between species environmental space Z4[Z4<max(Z4)/1000]<-0
    Zdiff.S.cor.sum <- sum(Zdiff.S.cor)
    Zdiff.S.cor.max <- max(Zdiff.S.cor)
    Zdiff.S.cor.min <- abs(min(Zdiff.S.cor))
    Zdiff.S.cor.Z <- 1/max(Zdiff.S.cor.max, Zdiff.S.cor.min)
    if (Zdiff.S.cor.Z == Inf) {
        Zdiff.S.cor.Z = 0
    }
    if (Zdiff.S.cor.Z != 0) {
        Zdiff.S.cor.raw <- Zdiff.S.cor + 1
        Zdiff.S.cor.raw <- Zdiff.S.cor.raw/2
        # Zdiff.S.cor<-Zdiff.S.cor+0.5
        Zdiff.S.cor <- Zdiff.S.cor * Zdiff.S.cor.Z
        Zdiff.S.cor <- Zdiff.S.cor + 1
        Zdiff.S.cor <- Zdiff.S.cor/2
    }
    if (Zdiff.S.cor.Z == 0) {
        Zdiff.S.cor.raw <- Zdiff.S.cor + 0.5
        Zdiff.S.cor <- Zdiff.S.cor + 0.5
    }
    
    ### tests for difference - cor
    D.uncor = 0
    D.cor = 0
    if (Zdiff.S.cor.sum != 0 & Zdiff.C.cor.sum != 0) {
        D.cor <- cor(c(Zdiff.C.cor), c(Zdiff.S.cor))
    }
    
    ### tests for difference - uncor
    if (Zdiff.S.cor.sum != 0 & Zdiff.C.cor.sum != 0) {
        D.uncor <- cor(c(Zdiff.C.uncor), c(Zdiff.S.uncor))
    }
    
    ############################## 
    if (Zdiff.S.cor.sum != 0 & Zdiff.C.cor.sum != 0) {
        print("Uncorrected Environment Overlap of Differences Between Environments and Species:")
        print(D.uncor)
        print("Corrected Environment Overlap of Differences Between Environments and Species:")
        print(D.cor)
        print("Analogous climate space percentage")
        print(ZspO$remaining.espace.per)
        
    }
    # if (D.uncor>=0.5 & Zdiff.S.uncor.sum!=0 & Zdiff.C.uncor.sum!=0 &
    # ZspO$remaining.espace.per<20 & ZspO$D*100<30){print('Your datasets have high overlap among
    # the differences between input environments and the differences between species
    # distributions in environmental space.  Normally, when highly correlated, it is strongly
    # encouraged to use 'correct.env=T' to correct espace observations by abundance of
    # environments.  However in this case, there exists very little similarity in analogous
    # espace among environments and this overlap is expected due to this.  Further, in shared
    # environments, the two species share little overlap. In this sitation, we recommend using
    # 'correct.env=F' for analyes.')}
    
    # if (D.uncor>=0.5 & Zdiff.S.uncor.sum!=0 & Zdiff.C.uncor.sum!=0 &
    # ZspO$remaining.espace.per<20 & ZspO$D*100>30){print('Your datasets have high overlap among
    # the differences between input environments and the differences between species
    # distributions in environmental space. If ignored, equivalency tests are prone to type I
    # errors and you may observed statistical differences that are entirely due to differences in
    # the availability of environments vs. actual differences in occupied environmental space. In
    # this case, though there exists very little similarity in analogous espace among
    # environments, which can cause these correlations. And beacause in shared environments, the
    # two species share considerable overlap, it is encouraged to use 'correct.env=T' to correct
    # espace observations by abundance of environments.')} & ZspO$remaining.espace.per>20
    
    if (D.uncor >= 0.5 & Zdiff.S.uncor.sum != 0 & Zdiff.C.uncor.sum != 0) {
        print("Your datasets have high overlap among the differences between input environments and the differences between species distributions in environmental space. This is despite the fact that considserable e-space is shared between the two environments. If ignored, equivalency tests are prone to type I errors and you may observed statistical differences that are entirely due to differences in the availability of environments vs. actual differences in occupied environmental space. When highly correlated and their exists considerable e-space overlap, it is strongly encouraged to use 'correct.env=T' to correct espace observations by abundance of environments.")
    }
    if (Zdiff.C.uncor.sum == 0 & Zdiff.S.uncor.sum > 0) {
        print("Input environments are identical, use either 'correct.env=T' or 'correct.env=F'.")
    }
    if (Zdiff.C.uncor.sum == 0 & Zdiff.S.uncor.sum == 0) {
        print("Input environments and observation points are identical, no need for analysis- the observed niches are identical.")
    }
    if (Zdiff.C.uncor.sum > 0 & Zdiff.S.uncor.sum == 0) {
        print("Input observation points are identical, no need for analysis- the observed niches are identical.")
    }
    l$D.overlap.uncorrected <- D.uncor
    l$D.overlap.corrected <- D.cor
    l$env1x <- Z.env1$x
    l$env1y <- Z.env1$y
    l$env2x <- Z.env2$x
    l$env2y <- Z.env2$y
    l$s.uncor <- Zdiff.S.uncor
    l$s.cor <- Zdiff.S.cor
    l$e.uncor <- Zdiff.C.uncor
    l$e.cor <- Zdiff.C.cor
    l$s.nonzero <- Zdiff.S.uncor.Z
    l$e.nonzero <- Zdiff.C.uncor.Z
    l$s.uncor.sum <- Zdiff.S.uncor.sum
    l$s.cor.sum <- Zdiff.S.cor.sum
    l$e.uncor.sum <- Zdiff.C.uncor.sum
    l$e.cor.sum <- Zdiff.C.cor.sum
    
    return(l)
}

############################################################################################
#################### Plot differences in espace between environments and species' distributions ####################
#################################################################################################
#' Plot differences in espace between environments and species' distributions
#' @param espace.diff output from humboldt.espace.correction
#' @param correct.env if correct.env=T, this function displays the espace differences between occurrence densities of both species or differences between both environments, which are corrected by the prevalence of the corresponding environments in their range. If correct.env=F, this function displays the espace differences between occurrence densities of both species or differences between both environments without correcting the prevalence of the corresponding environments in their range.
#' @param type if type="species", plots will show differences in espace between the two input species.  if type="environments", plots will show differences in espace among environments.  
#' @return Function plots differences observed in species or environmental espace. Used to determine if users should correct occurrence densities of each species by the prevalence of the environments in their range for equivalency, background and overlap analyses (correct.env=T). Often datasets have high overlap among the differences between input environments and the differences between species distributions in environmental space. If ignored, equivalency tests are prone to type I errors and you may observed statistical differences that are entirely due to differences in the availability of environments vs. actual differences in occupied environmental space. When highly correlated, it is strongly encouraged to use 'correct.env=T' to correct espace observations by abundance of environments"
#' @seealso \code{humboldt.g2e, humboldt.grid.espace, humboldt.equivalency.test, humboldt.background.test, humboldt.niche.overlap, humboldt.plot.niche,humboldt.doitall} which use or depend on outputs of this function 
#' @importFrom scales alpha
#' @export
#' @examples
#' library(humboldt)
#'
#' ##load environmental variables for all sites of the study area 1 (env1). Column names should be x,y,X1,X2,...,Xn)
#' env1<-read.delim("env1.txt",h=T,sep="\t")
#'
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' env2<-read.delim("env2.txt",h=T,sep="\t") 
#'
#' ## remove NAs and make sure all variables are imported as numbers
#' env1<-humboldt.scrub.env(env1)
#' env2<-humboldt.scrub.env(env2)
#'
#' ##load occurrence sites for the species at study area 1 (env1). Column names should be 'sp', 'x','y'
#' occ.sp1<-na.exclude(read.delim("sp1.txt",h=T,sep="\t"))
#'
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be 'sp', 'x','y'. 
#' occ.sp2<-na.exclude(read.delim("sp2.txt",h=T,sep="\t"))
#' 
#' ##convert geographic space to espace
#' zz<-humboldt.g2e(env1=env1, env2=env2, sp1=occ.sp1, sp2=occ.sp2, reduce.env = 2, reductype = "PCA", non.analogous.environments = "NO", env.trim.mcp = T, e.var=c(3:21),  col.env = e.var, mcp.buffer.sp1 = 200, mcp.buffer.sp2 = 200, env.pca = T, rarefy.dist = 50, rarefy.units="km", env.reso=0.41666669, kern.smooth = 1, PROJ = F,  R = 100, run.silent = F)
#' 
#' ##store espace socres for sp1 and environments 1,2 and both environments combined output from humboldt.g2e
#' scores.env1<-zz$scores.env1[1:2]
#' scores.env2<-zz$scores.env2[1:2]
#' scores.env12<- rbind(zz$scores.env1[1:2],zz$scores.env2[1:2])
#' scores.sp1<-zz$scores.sp1[1:2]
#' scores.sp2<-zz$scores.sp2[1:2]
#' 
#' ## run create a grid of Environmental Space Function
#' z.sp1<- humboldt.grid.espace(scores.env12,scores.env1,scores.sp1,kern.smooth=1,R=100)
#' z.sp2<- humboldt.grid.espace(scores.env12,scores.env2,scores.sp2,kern.smooth=1,R=100)
#' z.env1<- humboldt.grid.espace(scores.env12,scores.env1,scores.env1,kern.smooth=1,R=100)
#' z.env2<- humboldt.grid.espace(scores.env12,scores.env2,scores.env2,kern.smooth=1,R=100)
#' 
#' ee<- humboldt.espace.correction(Z.env1=z1,Z.env2=z2,Z.sp1=z3,Z.sp2=z4)
#' 
#' ## plot differences between species' espaces
#' humboldt.plot.espace.diff<-espace.diff=ee, correct.env=F, type="species")
#' 
#' ##plot contour lines of environment 1 if env1 and env2 are not identical
#' if(ee$s.uncor.sum!=0){
#' contour(z.env1$x,(sort((z.env1$y*-1))),z.env1$Z,add=T,levels=quantile(z.env1$Z[z.env1$Z>0],c(0.1,0.5,0.75)),drawlabels=F,lty=c(1,2,3), lwd=c(1,1,1), col="grey")}
#' 
#' ## plot differences between environments' espaces
#' humboldt.plot.espace.diff<-espace.diff=ee, correct.env=F, type="environments")
#'
#' ##plot contour lines of environmental 1 if env1 and env2 are not identical
#' if(ee$e.uncor.sum!=0){
#' contour(z.env1$x,(sort((z.env1$y*-1))),z.env1$Z,add=T,levels=quantile(z.env1$Z[z.env1$Z>0],c(0.1,0.5,0.75)),drawlabels=F,lty=c(1,2,3), lwd=c(1,1,1), col="grey")}

humboldt.plot.espace.diff <- function(espace.diff, correct.env = F, type = "species") {
    if (correct.env == T) {
        if (type == "species") {
            image(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$s.cor, 
                col = colorRampPalette(c("#000099", "#FFFFFF", "#FF3100"))(256), zlim = c(0, 
                  1), xlab = "PC1", ylab = "PC2")
            if (espace.diff$s.cor.sum != 0) {
                contour(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$s.cor, 
                  add = T, levels = pretty(espace.diff$s.cor[espace.diff$s.cor > 0], 10), drawlabels = F, 
                  col = alpha("black", 0.75))
                title(main = "Spp. Differences: E-space 1 vs. E-space 2", sub = paste("Red=S1>S2, Blue=S1<S2, Cor. Uncorrected=", 
                  round(espace.diff$D.overlap.uncorrected, 3)))
            }
        }
        if (type == "environment") {
            image(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$e.cor, 
                col = colorRampPalette(c("#000099", "#FFFFFF", "#FF3100"))(256), zlim = c(0, 
                  1), xlab = "PC1", ylab = "PC2")
            if (espace.diff$e.cor.sum != 0) {
                contour(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$e.cor, 
                  add = T, levels = pretty(espace.diff$e.cor[espace.diff$e.cor > 0], 10), drawlabels = F, 
                  col = alpha("black", 0.75))
                title(main = "Env. Differences: E-space 1 vs. E-space 2", sub = paste("Red=E1>E2, Blue=E1<E2, Cor. Uncorrected=", 
                  round(espace.diff$D.overlap.uncorrected, 3)))
            }
        }
    }
    if (correct.env == F) {
        if (type == "species") {
            image(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$s.uncor, 
                col = colorRampPalette(c("#000099", "#FFFFFF", "#FF3100"))(256), zlim = c(0, 
                  1), xlab = "PC1", ylab = "PC2")
            if (espace.diff$s.uncor.sum != 0) {
                contour(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$s.uncor, 
                  add = T, levels = pretty(espace.diff$s.uncor[espace.diff$s.uncor > 0], 10), 
                  drawlabels = F, col = alpha("black", 0.75))
            }
            title(main = "Spp. Differences: E-space 1 vs. E-space 2", sub = paste("Red=S1>S2, Blue=S1<S2, Cor. Uncorrected=", 
                round(espace.diff$D.overlap.uncorrected, 3)))
        }
        if (type == "environment") {
            image(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$e.uncor, 
                col = colorRampPalette(c("#000099", "#FFFFFF", "#FF3100"))(256), zlim = c(0, 
                  1), xlab = "PC1", ylab = "PC2")
            
            if (espace.diff$e.uncor.sum != 0) {
                contour(espace.diff$env1x, (sort((espace.diff$env1y * -1))), z = espace.diff$e.uncor, 
                  add = T, levels = pretty(espace.diff$e.uncor[espace.diff$e.uncor > 0], 10), 
                  drawlabels = F, col = alpha("black", 0.75))
                title(main = "Env. Differences: E-space 1 vs. E-space 2", sub = paste("Red=E1>E2, Blue=E1<E2, Cor. Uncorrected=", 
                  round(espace.diff$D.overlap.uncorrected, 3)))
            }
        }
    }
}
