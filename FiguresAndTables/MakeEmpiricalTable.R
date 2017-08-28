library(ape)
load("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Cichlid/CichlidResults.RSave")
cichlid.phy <- phy
load("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Nicotiana/Nicotiana.RSave")
nicotiana.phy <- phy

MakePrettyTable <- function(results, model.names=c("vH free", "beta free", "vH and beta fixed", "vH and beta free", "average"), params = c("sigma.sq", "mu", "bt", "vh", "SE", "sigma.sq_times_treeheight", "SE2"), round.factor=3, phy) {
	results$sigma.sq_times_treeheight <- results$sigma.sq * max(branching.times(phy))
	results$sigma.sq_times_treeheight.lower <- results$sigma.sq.lower * max(branching.times(phy))
	results$sigma.sq_times_treeheight.upper <- results$sigma.sq.upper * max(branching.times(phy))
	results$SE2 <- (results$SE)^2
	results$SE2.lower <- (results$SE.lower)^2
	results$SE2.upper <- (results$SE.upper)^2

	results.with.avg <- rbind(results, rep(NA, dim(results)[2]))
	model.avg.index <- dim(results.with.avg)[1]
	weights <- results$AkaikeWeight

	for (i in sequence(dim(results)[2])) {
		results.with.avg[model.avg.index, i] <- weighted.mean(results[,i], weights)	
	}
	
	pretty <- results.with.avg[,c("NegLogL", "K", "deltaAICc", "AkaikeWeight")]
	pretty <- cbind(data.frame(Model=model.names), pretty)
	pretty.orig.width <- dim(pretty)[2]
	pretty <- cbind(pretty, data.frame(matrix(nrow=dim(pretty)[1], ncol=length(params))))
	for (i in sequence(length(params))) {
		for (j in sequence(dim(pretty)[1])) {
			summary <- paste(round(results.with.avg[j,params[i]], round.factor), " (", round(results.with.avg[j,paste(params[i], ".lower", sep="")], round.factor), ", ", round(results.with.avg[j,paste(params[i], ".upper", sep="")], round.factor), ")", sep="")
			pretty[j, i+ pretty.orig.width] <- summary
		}	
		colnames(pretty)[i+ pretty.orig.width] <- params[i]
	}
	pretty<- pretty[order(pretty$K),]
	return(pretty)	
}

output.table <- rbind(MakePrettyTable(results.cichlid$results, phy=cichlid.phy), MakePrettyTable(results.nicotiana$results, phy=nicotiana.phy))
output.table <- cbind(Dataset=c(rep("Cichlid", 5), rep("Nicotiana", 5)), output.table)
write.csv(output.table, file="/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/FiguresAndTables/EmpiricalResults.csv")