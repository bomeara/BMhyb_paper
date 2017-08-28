rm(list=ls())
library(BMhyb)
library(ape)
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/Simulations")
#system("rm PDF*pdf")
source("GeneratePossibleSets.R")
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
best.results <- data.frame()
true.vs.estimated.vcv <- list()
current.file.count=0
all.files <- c(system(paste0("ls -1 *Try9*OPT_and_GRID.RData"),intern=TRUE), system(paste0("ls -1 *Try9*GRID_only.RData"),intern=TRUE))
for (j in sequence(dim(possible.sets)[1])) {
	id <- j
	#files <- system(paste0("ls -1 *Try9*ParamCombo",id,".P*SUCCESS.RData"),intern=TRUE, ignore.stderr=TRUE)
	files <- all.files[grepl(paste0("ParamCombo",id,".P"), all.files)]
	for (file.index in sequence(length(files))) {
		rm(result1)
		rm(result2)
		rm(result3)
		rm(result4)
		rm(result.grid)
		rm(result.grid.from.starting)
		print(files[file.index])

		try(load(files[file.index]))
		if(!is.null(result.grid) & class(result.grid)!="try-error") {
			if(nrow(result.grid$results)>0) {
			    #print(result4)
				file.split <- unlist(strsplit(strsplit(files[file.index],"_")[[1]],"\\."))
				current.file.count <- current.file.count + 1
				if(current.file.count%%5==0) {
					print(paste0("Done ", current.file.count, " of ", length(all.files), " files to do"))
				}
				local.results <- result.grid$results
				run.type <- "none"
				if(grepl("GRID_only", files[file.index])) {
					local.results <- result.grid.from.starting$results
					run.type <- "grid"
				} else {
					try(local.results <- rbind(result1, result2, result3, result4, result.grid.from.starting$results))
					run.type <- c(rep("optim",4), rep("grid", 4))
				}
				if(!is.null(local.results)) {
					local.results$run.type <- run.type
					#try(local.results <- rbind(local.results, result.best.with.measurement.error))
					local.results$combo = as.numeric(id)
					local.results$computer = file.split[1]
					local.results$pid = as.numeric(gsub("GRID","",gsub("PID", "", gsub("SUCCESS", "", file.split[4]))))
					#print(local.results)
					local.results <- local.results[,-(which(grepl("lower", colnames(local.results))))]
					local.results <- local.results[,-(which(grepl("upper", colnames(local.results))))]
					local.results$deltaAICc <- local.results$AICc - min(local.results$AICc)
					rel.lik <- exp(-0.5* local.results$deltaAICc)
					local.results$AkaikeWeight <- rel.lik / sum(rel.lik)
					local.results$flag.as.weird=FALSE
					local.results$convergence.possible.fail=FALSE
					if(max(local.results$deltaAICc>200)) {
						local.results$flag.as.weird=TRUE
						#print("possible weirdness")
						#print(local.results[,1:13])
					}
					if(max(local.results$NegLogL)>1e100) {
						local.results$convergence.possible.fail=TRUE
						local.results$flag.as.weird=TRUE
						#print("possible lack of convergence")
						#print(local.results[,1:13])
					}
					local.average <- local.results[-(2:4),]
					for (i in 1:10) {
						local.average[1,i] <- weighted.mean(local.results[,i], w=local.results$AkaikeWeight)
					}
					try(local.results <- cbind(local.results, possible.sets.renamed[local.results$combo[1],]))
					local.results$origin.type.true <- as.character(local.results$origin.type.true)
					local.average <- cbind(data.frame(local.average, stringsAsFactors=FALSE), local.results[1,c((length(local.average)+1):dim(local.results)[2])])
					counts.of.results[as.numeric(id)] <- counts.of.results[as.numeric(id)]+1
					try(counts.per.condition <- cbind(counts.of.results, possible.sets))

					free.parameters<-rep(TRUE, 5)
					names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
					estimated.params <- c(local.average$sigma.sq, local.average$mu, local.average$bt, local.average$vh, local.average$SE)
					names(estimated.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")

					best.params <- c(local.results$sigma.sq[which.min(local.results$deltaAICc)], local.results$mu[which.min(local.results$deltaAICc)], local.results$bt[which.min(local.results$deltaAICc)], local.results$vh[which.min(local.results$deltaAICc)], local.results$SE[which.min(local.results$deltaAICc)])
					names(best.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")

					estimated.VCV <- GetVModified(estimated.params, network$phy, network$flow, actual.params= free.parameters)
					estimated.means <- GetMeansModified(estimated.params, network$phy, network$flow, actual.params= free.parameters)
					estimated.likelihood <- NA
					try(estimated.likelihood <- CalculateLikelihood((estimated.params), tips, network$phy, network$flow, actual.params= free.parameters, lower.b = c(0, -Inf, 0.000001, 0, 0)[actual.params], upper.b=c(10,Inf,100,100,100)[actual.params]))

					best.VCV <- GetVModified(best.params, network$phy, network$flow, actual.params= free.parameters)
					best.means <- GetMeansModified(best.params, network$phy, network$flow, actual.params= free.parameters)
					best.likelihood <- NA
					try(best.likelihood <- CalculateLikelihood((best.params), tips, network$phy, network$flow, actual.params= free.parameters, lower.b= c(0, -Inf, 0.000001, 0, 0)[actual.params], upper.b=c(10,Inf,100,100,100)[actual.params]))


					true.params <- c(local.average$sigma.sq.true, local.average$mu.true, local.average$bt.true, local.average$vh.true, local.average$SE.true)
					names(true.params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
					true.VCV <- GetVModified(true.params, network$phy, network$flow, actual.params= free.parameters)
					true.means <- GetMeansModified(true.params, network$phy, network$flow, actual.params= free.parameters)
					true.likelihood <- NA
					try(true.likelihood <- CalculateLikelihood((true.params), tips, network$phy, network$flow, actual.params= free.parameters, lower.b = c(0, -Inf, 0.000001, 0, 0)[actual.params], upper.b=c(10,Inf,100,100,100)[actual.params]))


					true.min.eigen <- min(eigen(true.VCV)$values)
					estimated.min.eigen <- min(eigen(estimated.VCV)$values)
					best.min.eigen <- min(eigen(best.VCV)$values)
					phy.only.min.eigen <- min(eigen(ape::vcv(network$phy))$values)


					local.average$true.min.eigen = true.min.eigen
					local.average$estimated.min.eigen = estimated.min.eigen
					local.results$true.min.eigen = true.min.eigen
					local.results$estimated.min.eigen = estimated.min.eigen
					local.average$true.is.ultrametric = TestUltrametric(true.VCV)
					local.average$estimated.is.ultrametric = TestUltrametric(estimated.VCV)
					local.results$true.is.ultrametric = TestUltrametric(true.VCV)
					local.results$estimated.is.ultrametric = TestUltrametric(estimated.VCV)
					local.results$true.kappa = kappa(true.VCV)
					local.results$estimated.kappa = kappa(estimated.VCV)
					local.results$true.det = det(true.VCV)
					local.results$estimated.det = det(estimated.VCV)

					local.average$best.min.eigen = best.min.eigen
					local.results$best.min.eigen = best.min.eigen
					local.average$best.is.ultrametric = TestUltrametric(best.VCV)
					local.results$best.is.ultrametric = TestUltrametric(best.VCV)
					local.results$best.kappa = kappa(best.VCV)
					local.results$best.det = det(best.VCV)
					local.results$phy.only.min.eigen <- phy.only.min.eigen
					#print(c(true.min.eigen , phy.only.min.eigen))

					local.results$true.likelihood <- true.likelihood
					local.results$best.likelihood <- best.likelihood
					local.results$estimated.likelihood <- estimated.likelihood
					pdf(file=paste0("PDF_",ifelse(true.min.eigen < 0, "BAD", "GOOD"),"_",gsub("RData", "pdf", files[file.index])))
					PlotNetwork(network$phy, network$flow)
					dev.off()


					try(all.results <- rbind(all.results, local.results))
					try(average.results <- rbind(average.results, local.average))
					save(all.results, average.results, best.results, counts.of.results, counts.per.condition , file="SimulationResults.RData")
				}
			}
		}
	}
}


best.results <- subset(all.results, deltaAICc==0)
save(all.results, best.results, average.results, counts.of.results, counts.per.condition, possible.sets.renamed, possible.sets, file="SimulationResults.RData")

#print(table(all.results$flag.as.weird)/dim(all.results)[1])
#print(table(all.results$convergence.possible.fail)/dim(all.results)[1])
