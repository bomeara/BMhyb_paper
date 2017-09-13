library(ape)
load("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidAug2017/CichlidResults_AllCombined.RData")
cichlid.phy <- read.tree("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidAug2017/cichlid.phy")
results.cichlid <- all.results
load("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/NicotianaAug2017/NicotianaResults_AllCombined.RData")
nicotiana.phy <- read.tree("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/NicotianaAug2017/nicotiana.phy")
results.nicotiana<- all.results



MakePrettyTable <- function(results, model.names=c("vH free", "beta free", "vH and beta fixed", "vH and beta free"), params = c("sigma.sq", "mu", "bt", "vh", "SE", "sigma.sq_times_treeheight", "SE2"), round.factor=1, phy) {
	results$sigma.sq_times_treeheight <- unlist(results$sigma.sq) * max(branching.times(phy))
	results$sigma.sq_times_treeheight.lower <- unlist(results$sigma.sq.lower) * max(branching.times(phy))
	results$sigma.sq_times_treeheight.upper <- unlist(results$sigma.sq.upper) * max(branching.times(phy))
	results$SE2 <- (unlist(results$SE))^2
	results$SE2.lower <- (unlist(results$SE.lower))^2
	results$SE2.upper <- (unlist(results$SE.upper))^2

	results.with.avg <- rbind(results, rep(NA, dim(results)[2]))
	model.avg.index <- dim(results.with.avg)[1]
	weights <- results$AkaikeWeight

	for (i in sequence(dim(results)[2])) {
		try(results.with.avg[model.avg.index, i] <- weighted.mean(unlist(results[,i]), weights))
	}
	if(!any(grepl("error", colnames(results.with.avg)))) {
		results.with.avg$error <- "unspecified"
	}

	pretty <- results.with.avg[,c("error", "NegLogL", "K", "deltaAICc", "AkaikeWeight")]
	pretty <- cbind(data.frame(ModelName=model.names[unlist(results.with.avg$Model)], stringsAsFactors=FALSE), pretty)
	pretty$ModelName[nrow(pretty)] <- "average"
	pretty.orig.width <- dim(pretty)[2]
	pretty <- cbind(pretty, data.frame(matrix(nrow=dim(pretty)[1], ncol=length(params))))
	for (i in sequence(length(params))) {
		for (j in sequence(dim(pretty)[1])) {
			summary <- paste(round(unlist(results.with.avg[j,params[i]]), round.factor), " (", round(unlist(results.with.avg[j,paste(params[i], ".lower", sep="")]), round.factor), ", ", round(unlist(results.with.avg[j,paste(params[i], ".upper", sep="")]), round.factor), ")", sep="")
			pretty[j, i+ pretty.orig.width] <- summary
		}
		colnames(pretty)[i+ pretty.orig.width] <- params[i]
	}
	pretty$NegLogL <- round(unlist(pretty$NegLogL), round.factor)
	pretty$K <- round(unlist(pretty$K), 0)
	pretty$deltaAICc <- round(unlist(pretty$deltaAICc), round.factor)
	pretty$AkaikeWeight <- round(unlist(pretty$AkaikeWeight), round.factor)


	pretty$deltaAICc[which(pretty$ModelName=="average")] <- 1e1000
	pretty<- pretty[order(pretty$deltaAICc),]
	pretty$deltaAICc[which(pretty$ModelName=="average")] <- NA
	pretty$NegLogL[which(pretty$ModelName=="average")] <- NA
	pretty$K[which(pretty$ModelName=="average")] <- NA
	pretty$AkaikeWeight[which(pretty$ModelName=="average")] <- NA
	pretty$SE[which(pretty$error=="Empirical")] <- "Empirical"

	return(pretty)
}
cichlid.pretty <- MakePrettyTable(results.cichlid, phy=cichlid.phy)
nicotiana.pretty <- MakePrettyTable(results.nicotiana, phy=nicotiana.phy)
output.table <- rbind(cichlid.pretty , nicotiana.pretty)
output.table <- cbind(Dataset=c(rep("Cichlid", nrow(cichlid.pretty)), rep("Nicotiana", nrow(nicotiana.pretty ))), output.table)
write.csv(output.table, file="/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/FiguresAndTablesAug2017/EmpiricalResults.csv")
