# print(system("hostname", intern=TRUE))
# 
# tmp.install.packages <- function(pack, dependencies=TRUE, ...) {
#   path <- tempdir()
#   ## Add 'path' to .libPaths, and be sure that it is not
#   ## at the first position, otherwise any other package during
#   ## this session would be installed into 'path'
#   firstpath <- .libPaths()[1]
#   .libPaths(c(firstpath, path))
#   install.packages(pack, dependencies=dependencies, lib=path, ...)
# }
# list.of.packages <- c("BMhyb")
# if((list.of.packages %in% installed.packages()[,"Package"])) {
# 	if(packageVersion("BMhyb")!="1.3.3") {
# 		tmp.install.packages("BMhyb_1.3.3.tar.gz", repos = NULL, type="source")
# 	}
# }
# if(!(list.of.packages %in% installed.packages()[,"Package"])) {
# 	tmp.install.packages("BMhyb_1.3.3.tar.gz", repos = NULL, type="source")
# }
#library(BMhyb)

rm(list=ls())
setwd("~/Dropbox/CollabJhwuengOMeara/ReReSubmission/TestModifiedModel")
source("bmhyb.r")
source("GeneratePossibleSets.R")

#load("~/Dropbox/CollabJhwuengOMeara/ReReSubmission/Simulations/RunToLookAtSurfaceWithoutSearchingRun2_clusterl/ParamCombo273.Rep8.RData")


#names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#GetMeansModified <- function(x, phy, flow, actual.params) {
#params must be named vector
#SimulateTipData <- function(phy, flow, params) {

#SimulateNetwork <- function(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE) {

#RunFromSet<-function(x, id, rep.id, possible.sets) {
#	print("starting")
#	print(c(id, rep.id))
  


for(Index in 1:dim(possible.sets)[1]){
  print(Index)
  x<-possible.sets[Index,]
  params<-rep(NA,5)
	names(params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	for (i in sequence(length(params))) {
		params[i]<-x[which(names(x)==names(params)[i])]
	  }
  actual.params<-rep(TRUE, 5)
  names(actual.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	params<-unlist(params)
	
	network<-SimulateNetwork(ntax.nonhybrid=as.numeric(x['ntax.nonhybrid']), ntax.hybrid=as.numeric(x['ntax.hybrid']), flow.proportion=as.numeric(x['flow.proportion']), origin.type=as.character(x['origin.type']), birth = 1, death = 0.5, sample.f = 0.5, tree.height = as.numeric(x['tree.height']), allow.ghost=FALSE)
	try(tips<-SimulateTipData(network$phy, network$flow, params))
  V.modified <- GetVModified(params, network$phy, network$flow, actual.params)
  means.modified <- GetMeansModified(params, network$phy, network$flow, actual.params)
  kappa.val <- kappa(V.modified, exact=TRUE)
  
  BMhyb(data=tips, phy=network$phy, flow=network$flow, opt.method="Nelder-Mead", models=c(1,2,3,4), verbose=FALSE, get.se=FALSE, plot.se=FALSE, store.sims=FALSE, precision=2, auto.adjust=FALSE, likelihood.precision=0.001, allow.extrapolation=FALSE, n.points=5000, measurement.error=NULL, do.kappa.check=FALSE, number.of.proportions=101, number.of.proportions.adaptive=101, allow.restart=TRUE) 
  
  
  
# 	try(likelihood.extrap.no.kappa.no <- CalculateLikelihood(params, tips, network$phy, network$flow, actual.params, allow.extrapolation=FALSE, do.kappa.check=FALSE))
# 	try(likelihood.extrap.yes.kappa.no <- CalculateLikelihood(params, tips, network$phy, network$flow, actual.params, allow.extrapolation=TRUE, do.kappa.check=FALSE))
# 	try(likelihood.extrap.no.kappa.yes<- CalculateLikelihood(params, tips, network$phy, network$flow, actual.params, allow.extrapolation=FALSE, do.kappa.check=TRUE))
# 	try(likelihood.extrap.yes.kappa.yes<- CalculateLikelihood(params, tips, network$phy, network$flow, actual.params, allow.extrapolation=TRUE, do.kappa.check=TRUE))
# 
#   likelihood.extrap.no.kappa.no
#   likelihood.extrap.yes.kappa.no
#   likelihood.extrap.no.kappa.yes
#   likelihood.extrap.yes.kappa.yes
  
  	#save(list=ls(), file=paste("ParamCombo", id, ".Rep", rep.id, ".RData", sep=""))
}


rm(list=ls())
setwd("~/Dropbox/CollabJhwuengOMeara/ReReSubmission/Simulations/simsSig06_Expect100")
load("sim.tree.est.1.Rsave")
bd.est
plot(bd.est$phy)
yule.est
plot(yule.est$phy)
