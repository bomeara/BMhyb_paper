rm(list=ls())
library(BMhyb)
setwd("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Simulations/RunWithFewerFlowParamsv1.3.3_cluster_l")
source("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Simulations/RunWithFewerFlowParamsv1.3.3_cluster_l/GeneratePossibleSets.R")
counts.of.results <- rep(0, dim(possible.sets)[1])
possible.sets.renamed <- possible.sets
colnames(possible.sets.renamed) <- paste(colnames(possible.sets.renamed), ".true", sep="")

TestUltrametric <- function(x) {
	ultrametric=TRUE
	for (i in sequence(dim(x)[1])) {
		for (j in sequence(dim(x)[2])) {
			if(x[i,j] < min(c(x[i,], x[,j]))) {
				ultrametric=FALSE
			}
			if (i!=j) {
				if(x[i,i] < x[i,j]) {

					ultrametric=FALSE
				}
			}
		}
	}	
	return(ultrametric)
}

all.results <- data.frame()
average.results <- data.frame()
true.vs.estimated.vcv <- list()
for(rep.id in sequence(nreps)) {
	id.converter <- sample.int(dim(possible.sets)[1])
	for (j in sequence(dim(possible.sets)[1])) {
		id <- id.converter[j] #this is done to reorder which params are run first
		file <- paste("ParamCombo", id, ".Rep", rep.id, ".RData", sep="")
		try(load(file), silent=TRUE)
		if(any(grepl("model.1", ls())) & any(grepl("model.2", ls())) & any(grepl("model.3", ls())) & any(grepl("model.4", ls()))) {
			model.2$Model=2
			model.3$Model=3
			model.4$Model=4
			local.results <- rbind(model.1, model.2, model.3, model.4)
			local.results$combo = as.numeric(id)
			local.results$rep.id = as.numeric(rep.id)
			local.results <- local.results[,-(which(grepl("lower", colnames(local.results))))]
			local.results <- local.results[,-(which(grepl("upper", colnames(local.results))))]
			local.results$deltaAICc <- local.results$AICc - min(local.results$AICc)
			rel.lik <- exp(-0.5* local.results$deltaAICc)
			local.results$AkaikeWeight <- rel.lik / sum(rel.lik)
			local.results$flag.as.weird=FALSE
			local.results$convergence.possible.fail=FALSE
			
			if(max(local.results$deltaAICc>200)) {
				local.results$flag.as.weird=TRUE
				print("possible weirdness")
				print(local.results[,1:13])	
			}
			if(max(local.results$NegLogL)>1e100) {
				local.results$convergence.possible.fail=TRUE
				local.results$flag.as.weird=TRUE
				print("possible lack of convergence")
				print(local.results[,1:13])	
			}
			local.average <- local.results[-(2:4),]
			for (i in 1:10) {
				local.average[1,i] <- weighted.mean(local.results[,i], w=local.results$AkaikeWeight)	
			}

			
			local.results <- cbind(local.results, possible.sets.renamed[local.results$combo[1],])
			local.results$origin.type.true <- as.character(local.results$origin.type.true)
			local.average <- cbind(data.frame(local.average, stringsAsFactors=FALSE), local.results[1,c((length(local.average)+1):dim(local.results)[2])])
			counts.of.results[as.numeric(id)] <- counts.of.results[as.numeric(id)]+1
			counts.per.condition <- cbind(counts.of.results, possible.sets)
			
		  	free.parameters<-rep(TRUE, 5)
	  		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
			estimated.params <- c(local.average$sigma.sq, local.average$mu, local.average$bt, local.average$vh, local.average$SE)
			names(estimated.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
			estimated.VCV <- GetVModified(estimated.params, network$phy, network$flow, actual.params= free.parameters)
			estimated.means <- GetMeansModified(estimated.params, network$phy, network$flow, actual.params= free.parameters)
			estimated.likelihood <- CalculateLikelihood(estimated.params, tips, network$phy, network$flow, actual.params= free.parameters)

			true.params <- c(local.average$sigma.sq.true, local.average$mu.true, local.average$bt.true, local.average$vh.true, local.average$SE.true)
			names(true.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
			true.VCV <- GetVModified(true.params, network$phy, network$flow, actual.params= free.parameters)
			true.means <- GetMeansModified(true.params, network$phy, network$flow, actual.params= free.parameters)
			true.likelihood <- CalculateLikelihood(true.params, tips, network$phy, network$flow, actual.params= free.parameters)
			
	
			true.min.eigen <- min(eigen(true.VCV)$values)
			estimated.min.eigen <- min(eigen(estimated.VCV)$values)
			
			local.average$true.min.eigen = true.min.eigen
			local.average$estimated.min.eigen = estimated.min.eigen
			local.results$true.min.eigen = true.min.eigen
			local.results$estimated.min.eigen = estimated.min.eigen
						local.average$true.is.ultrametric = TestUltrametric(true.VCV)
			local.average$estimated.is.ultrametric = TestUltrametric(estimated.VCV)
			local.results$true.is.ultrametric = TestUltrametric(true.VCV)
			local.results$estimated.is.ultrametric = TestUltrametric(estimated.VCV)

						all.results <- rbind(all.results, local.results)
			average.results <- rbind(average.results, local.average)

			save(all.results, counts.of.results, counts.per.condition , file="/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/FiguresAndTables/SimulationResults.RData")
			print(paste("loaded", sum(counts.of.results), "results"))

			rm(model.1, model.2, model.3, model.4)
		}
	}
}

best.results <- subset(all.results, deltaAICc==0)
save(all.results, best.results, average.results, counts.of.results, counts.per.condition, possible.sets.renamed, possible.sets, file="/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/FiguresAndTables/SimulationResults.RData")

print(table(all.results$flag.as.weird)/dim(all.results)[1])
print(table(all.results$convergence.possible.fail)/dim(all.results)[1])