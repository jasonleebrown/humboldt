#' Example Data: Conium maculatum's invasive distribution in Americas
#' 
#' @docType data
#'
#' @usage data(sp2)
#'
#' @format An data frame with 482 rows and 3 variables:
#' \describe{
#'   \item{sp}{species name}
#'   \item{x}{longitude}
#'   \item{y}{latittude}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' 
#' library(humboldt)
#' 
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
#' occ.sp1<-sp1
#' 
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
#' data(sp2)
#' occ.sp2<-sp2

"sp2"

#' Example Data: Conium maculatum's native distribution
#' 
#' @docType data
#'
#' @usage data(sp1)
#'
#' @format An data frame with 4975 rows and 3 variables:
#' \describe{
#'   \item{sp}{species name}
#'   \item{x}{longitude}
#'   \item{y}{latittude}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' 
#' library(humboldt)
#' 
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
#' occ.sp1<-sp1
#' 
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
#' data(sp2)
#' occ.sp2<-sp2

"sp1"

#' Example Data: Bioclimatic data for Europe and surrounding areas
#' 
#' @docType data
#'
#' @usage data(env2)
#'
#' @format An data frame with 23758 rows and 21 variables:
#' \describe{
#'   \item{x}{longitude}
#'   \item{y}{latittude}
#'   \item{var1}{any continious enviromental variable}
#'   \item{var2}{any continious enviromental variable}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' 
#' library(humboldt)
#' 
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
#' occ.sp1<-sp1
#' 
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
#' data(sp2)
#' occ.sp2<-sp2

"env2"

#' Example Data: Bioclimatic data for the Americas and surrounding areas
#' 
#' @docType data
#'
#' @usage data(env1)
#'
#' @format An data frame with 14095 rows and 21 variables:
#' \describe{
#'   \item{x}{longitude}
#'   \item{y}{latittude}
#'   \item{var1}{any continious enviromental variable}
#'   \item{var2}{any continious enviromental variable}
#' } 
#'
#' @keywords datasets
#'
#' @examples
#' 
#' library(humboldt)
#' 
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
#' occ.sp1<-sp1
#' 
#' ##load occurrence sites for the species at study area 2 (env2). Column names should be sp,x,y
#' data(sp2)
#' occ.sp2<-sp2

"env1"

