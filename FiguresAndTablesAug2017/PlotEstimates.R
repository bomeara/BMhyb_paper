setwd('/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/FiguresAndTablesAug2017')
library(BMhyb)
rm(list=ls())

load("../Simulations/SimulationResults.RData")
source("../Simulations/GeneratePossibleSets.R")

params.to.plot.main.fig <- c("bt", "vh")
params.to.plot.supp.fig <- c("bt", "vh", "mu", "SE", "sigma.sq")

logs.to.use.main.fig <-    c("y", "")
logs.to.use.supp.fig <-    c("y", "", "", "", "y")
hybrid.labels <- c("i", "v", "x")
nonhybrid.labels <- c("30", "100")

only.ok.runs <- FALSE
#origins.chosen <- "individual"
origins.chosen <- NULL

model1.only <- subset(all.results, all.results$Model==1)

print("Replicates per combo:")
print(t(t(table(model1.only$combo))))

if(only.ok.runs) {
	average.results <- subset(average.results, !flag.as.weird & !convergence.possible.fail)
	best.results <- subset(best.results, !flag.as.weird & !convergence.possible.fail)
	all.results <- subset(all.results, !flag.as.weird & !convergence.possible.fail)
}

if(!is.null(origins.chosen)) {
	average.results <- subset(average.results, origin.type.true==origins.chosen)
	best.results <- subset(best.results, origin.type.true==origins.chosen)
	all.results <- subset(all.results, origin.type.true==origins.chosen)
}


for (figure.index in sequence(2)) {
	params.to.plot <- c()
	pdf.name <- NA
	if(figure.index==1) {
		params.to.plot <- params.to.plot.main.fig
		pdf.name <- "FigureSimulationResults.pdf"
		logs.to.use <- logs.to.use.main.fig
		#summaries.by.estimator <- summaries.by.estimator.main
	} else {
		params.to.plot <- params.to.plot.supp.fig
		pdf.name <- "SuppFigureSimulationResults.pdf"
		logs.to.use <- logs.to.use.supp.fig
		#summaries.by.estimator <- summaries.by.estimator.supp
	}
	pdf(file=pdf.name, width=5*length(params.to.plot), height=5)
	par(mfcol=c(1, length(params.to.plot)))
	for (param.index in sequence(length(params.to.plot))) {

		possible.values.true <- sort(unique(average.results[,paste(params.to.plot[param.index],".true", sep="")]))
		values.ntip.nonhybrids <- sort(unique(average.results$ntax.nonhybrid.true))
		values.ntip.hybrids <- sort(unique(average.results$ntax.hybrid.true))
		values.true <- average.results[,paste(params.to.plot[param.index],".true", sep="")]
		values.estimated <- average.results[,paste(params.to.plot[param.index], sep="")]
		combos.of.parameters <- expand.grid(nonhybrids=values.ntip.nonhybrids, hybrids=values.ntip.hybrids, true.param=possible.values.true)
		combos.of.ntax <- expand.grid(nonhybrids=values.ntip.nonhybrids, hybrids=values.ntip.hybrids)
		y.range <- range(c(quantile(values.true, c(0.01, 0.99)), values.estimated))
		if(params.to.plot[param.index]=="vh") {
			y.range <- c(-5, 20)
		}
		if(params.to.plot[param.index]=="SE") {
			y.range <- c(-.2, 1.5)
		}
		
		if(params.to.plot[param.index]=="mu") {
		    y.range <- c(0, 2)
		}
		if(params.to.plot[param.index]=="bt") {
		    y.range <- c(1e-7, 1e3)
		}
		plot(x=c(0, dim(combos.of.parameters)[1]+3+length(values.ntip.nonhybrids)), y=y.range, type="n", xaxt="n", xlab="", ylab=params.to.plot[param.index], bty="n", log=logs.to.use[param.index], main=paste(params.to.plot[param.index]))
		base.position = 0
		for (possible.values.index in sequence(length(possible.values.true))) {
			true.value <- possible.values.true[possible.values.index]
			base.position <- base.position+1
			base.start = base.position
			base.end = base.start
			for (nonhybrids.index in sequence(length(values.ntip.nonhybrids))) {
				nonhybrids = values.ntip.nonhybrids[nonhybrids.index]
				if(nonhybrids.index>1) {
					base.position <- base.position+1
				}
				max.depth=Inf
				for (hybrids.index in sequence(length(values.ntip.hybrids))) {
					hybrids = values.ntip.hybrids[hybrids.index]
					true.column <- which(colnames(average.results)==paste(params.to.plot[param.index], ".true", sep=""))
					est.column <- which(colnames(average.results)==params.to.plot[param.index])
					local.df <- average.results[which(average.results[,true.column]==true.value),]
					local.df <- subset(local.df, ntax.nonhybrid.true== nonhybrids & ntax.hybrid.true==hybrids)
					values.y <- local.df[,est.column]
					quantile.results <- quantile(values.y, c(0.025, 0.975))

					points(base.position, median(values.y), pch=20)
					arrows(x0=base.position, y0=max(quantile.results[1], min(y.range)), x1=base.position, y1=min(quantile.results[2], max(y.range)), length=0, col="black")
					base.end=base.position
					max.depth = min(min(quantile.results), max.depth)
					text(base.position, max(quantile.results), hybrid.labels[hybrids.index], pos=3, col="gray")
					base.position <- base.position+1
				}
				text(base.position-2, max.depth, nonhybrid.labels[nonhybrids.index], pos=1, col="gray")

			}
			lines(c(base.start, base.end), c(true.value, true.value), col="gray")


		}
	}
	dev.off()
	system(paste("open ", pdf.name))
}

combo.counts <- rep(0, nrow(possible.sets))
for (i in sequence(length(combo.counts))) {
    combo.counts[i] <- nrow(subset(all.results, all.results$combo==i))   
}
possible.sets.and.results <- possible.sets
possible.sets.and.results$finished.sims <- combo.counts