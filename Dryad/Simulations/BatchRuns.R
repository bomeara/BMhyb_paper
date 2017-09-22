rm(list=ls())
library(BMhyb)
library(parallel)
library(R.utils)

DoRun <- function(id, identifier, sets) {
	session.info.variable <- sessionInfo()
	x <- sets[id,]
	params<-rep(NA,5)
	hostname <- system("hostname -s", intern=TRUE)
	names(params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	for (i in sequence(length(params))) {
		params[i]<-x[which(names(x)==names(params)[i])]
	}
	actual.params<-rep(TRUE, 5)
	names(actual.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	params<-unlist(params)
	result1 <- result2 <- result3 <- result4 <- NULL
	tips <- NULL
	attempts=1
	while(is.null(tips) & attempts < 2000) {
		try(network<-SimulateNetwork(ntax.nonhybrid=as.numeric(x['ntax.nonhybrid']), ntax.hybrid=as.numeric(x['ntax.hybrid']), flow.proportion=as.numeric(x['flow.proportion']), origin.type=as.character(x['origin.type']), birth = 1, death = 1, sample.f = 0.5, tree.height = as.numeric(x['tree.height']), allow.ghost=FALSE))
		try(tips<-SimulateTipData(network$phy, network$flow, params))

		attempts <- attempts+1
	}
	if(!is.null(tips)) {
		#try(result1 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=1, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE))
		#try(result2 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=2, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE))
		#try(result3 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=3, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE))
		#try(result4 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=4, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE))
		#try(result.best.with.measurement.error <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=which.min(c(result1$AICc, result2$AICc, result3$AICc, result4$AICc)), get.se=FALSE, measurement.error=NULL, attempt.deletion.fix=FALSE))
		result.grid <- NULL
		result.grid <- try(BMhybGrid(tips, network$phy, network$flow, plot.se=FALSE, store.sims=TRUE, measurement.error=0, attempt.deletion.fix=FALSE, n.points=10000))
		save(list=ls(), file=paste(hostname, "_Try9_ParamCombo", id, ".PID", identifier, "FirstGrid.RData", sep=""), compress=TRUE)
		pdf(file=paste(hostname, "_Try9_ParamCombo", id, ".PID", identifier, "SUCCESS.pdf", sep=""), width=20, height=8)
		try(PlotAICRegion(result.grid$sims, true.params=params))
		dev.off()
		starting.names <- colnames(result.grid$sims[which.min(result.grid$sims$AICc)[1],1:5])
		starting.values <- as.numeric(result.grid$sims[which.min(result.grid$sims$AICc)[1],1:5])
		names(starting.values) <- starting.names

			R.utils::withTimeout(try(result1 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=1, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE, starting.values=starting.values)), timeout=3600*4, onTimeout="warning")
			R.utils::withTimeout(try(result2 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=2, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE, starting.values=starting.values)), timeout=3600*4, onTimeout="warning")
			R.utils::withTimeout(try(result3 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=3, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE, starting.values=starting.values)), timeout=3600*4, onTimeout="warning")
			R.utils::withTimeout(try(result4 <- BMhyb(tips, network$phy, network$flow, plot.se=FALSE, models=4, get.se=FALSE, measurement.error=0, attempt.deletion.fix=FALSE, starting.values=starting.values)), timeout=3600*4, onTimeout="warning")


				result.grid.from.starting <- NULL
				R.utils::withTimeout(result.grid.from.starting <- try(BMhybGrid(tips, network$phy, network$flow, plot.se=FALSE, store.sims=TRUE, measurement.error=0, attempt.deletion.fix=FALSE, n.points=10000, starting.values=starting.values)), timeout=3600*4, onTimeout="warning")
				save(list=ls(), file=paste(hostname, "_Try9_ParamCombo", id, ".PID", identifier, "OPT_and_GRID.RData", sep=""), compress=TRUE)

				#save(list=ls(), file=paste(hostname, "_Try9_ParamCombo", id, ".PID", identifier, "GRID_only.RData", sep=""), compress=TRUE)
#	try(system(paste0('git add ', paste(hostname, "_ParamCombo", id, ".Rep", rep.id, ".RData", sep=""))))
#	try(system("git pull"))
#	try(system(paste0('git commit -m"',paste("ParamCombo", id, ".Rep", rep.id, " from ", hostname, " added", sep=""), '" -a')))
	}
}


source("GeneratePossibleSets.R")

#try(system("git lfs install"))
#try(system("git pull"))
#try(system('git lfs track "*.RData"'))
#try(system('git add .gitattributes'))
#try(system('git commit -m"adding lfs" -a'))

set.indices <- sample.int(nrow(possible.sets), nrow(possible.sets))

for (i in sequence(length(set.indices))) {
    DoRun(id=set.indices[i], identifier=Sys.getpid(), sets=possible.sets)
}

#for (sim.rep in sequence(nreps)) {
#	empty <- mclapply(set.indices, DoRun, rep.id=sim.rep, sets=possible.sets, mc.cores=24)
#}
