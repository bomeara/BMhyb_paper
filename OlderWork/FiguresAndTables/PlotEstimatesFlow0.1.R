setwd('~/Dropbox/CollabJhwuengOMeara/Simulations/CollabTony/')
load("AllSummarizedByRep.RData")
flow.to.plot <- 0.1
params.to.plot.main.fig <- c("bt", "vh")
params.to.plot.supp.fig <- c("bt", "vh", "mu", "SE", "sigma.sq")

logs.to.use.main.fig <-    c("y", "")
logs.to.use.supp.fig <-    c("y", "", "", "", "y")
hybrid.labels <- c("i", "v", "x")
nonhybrid.labels <- c("30", "100")

summaries$Nnonhybrid <- summaries$Ntax - summaries$Nhybrid
summaries$flow.prop0vs.1 <- 0
summaries$flow.prop0vs.1[which(summaries$flow.prop==0.1)] <- 1
summaries$flow.prop0vs.5 <- 0
summaries$flow.prop0vs.5[which(summaries$flow.prop==0.5)] <- 1
summaries <- subset(summaries, NumberModels==4 &flow.prop == flow.to.plot) #only do sims that completed; only do correct flow rate

bad.rows <- which(is.na(summaries), TRUE)[,1]
summaries <- summaries[-bad.rows,]



summaries.avg <-summaries[,c("tree.type", "Nnonhybrid", "Nhybrid", "flow.prop0vs.1", "flow.prop0vs.5", "bt", "vh", "mu", "SE", "sigma.sq", "bt.true", "vh.true", "mu.true", "SE.true", "sigma.sq.true")]
summaries.avg$estimator <- 0 #use avg model
summaries.avg$vh.fixed <- 0
summaries.vh.fixed <-summaries[,c("tree.type", "Nnonhybrid", "Nhybrid", "flow.prop0vs.1", "flow.prop0vs.5", "bt.vhfixed", "vh.vhfixed", "mu.vhfixed", "SE.vhfixed", "sigma.sq.vhfixed", "bt.true", "vh.true", "mu.true", "SE.true", "sigma.sq.true")]
summaries.vh.fixed$estimator <- 0
summaries.vh.fixed$vh.fixed <- 1


colnames(summaries.vh.fixed) <- gsub(".vhfixed", "", colnames(summaries.vh.fixed))

summaries.by.estimator <- list(summaries.avg, summaries.vh.fixed)
names(summaries.by.estimator) <- c("model average", "models with vH fixed at 0")
summaries.by.estimator.main <- summaries.by.estimator[-2]
summaries.by.estimator.supp <- summaries.by.estimator
for (figure.index in sequence(2)) {
	params.to.plot <- c()
	pdf.name <- NA
	if(figure.index==1) {
		params.to.plot <- params.to.plot.main.fig
		pdf.name <- "/Users/bomeara/Dropbox/CollabJhwuengOMeara/PCM_Hybridization/Resubmission/NewFiguresTables/SimulationSummaries/FigureSimulationResultsFlow0.1.pdf"
		logs.to.use <- logs.to.use.main.fig
		summaries.by.estimator <- summaries.by.estimator.main
	} else {
		params.to.plot <- params.to.plot.supp.fig
		pdf.name <- "/Users/bomeara/Dropbox/CollabJhwuengOMeara/PCM_Hybridization/Resubmission/NewFiguresTables/SimulationSummaries/SuppFigureSimulationResultsFlow0.1.pdf"
		logs.to.use <- logs.to.use.supp.fig
		summaries.by.estimator <- summaries.by.estimator.supp
	}
	pdf(file=pdf.name, width=5*length(params.to.plot), height=5*length(summaries.by.estimator))
	par(mfcol=c(length(summaries.by.estimator), length(params.to.plot)))
	for (param.index in sequence(length(params.to.plot))) {
		for (estimator.index in sequence(length(summaries.by.estimator))) {
			#pdf(file=paste("Figure_", params.to.plot[param.index], "_estimator_", names(summaries.by.estimator)[estimator.index], ".pdf", sep=""), width=5, height=5)
			local.summary <- summaries.by.estimator[[estimator.index]]
			values.true.x <- sort(unique(local.summary[,paste(params.to.plot[param.index],".true", sep="")]))
			values.ntip.nonhybrids <- sort(unique(local.summary$Nnonhybrid))
			values.ntip.hybrids <- sort(unique(local.summary$Nhybrid))
			observed.column.index <- which(colnames(local.summary)==paste(params.to.plot[param.index], sep=""))
			
			#plot(x=c(0,(3*length(values.true.x)+1)), y=range(local.summary[,observed.column.index], na.rm=TRUE), type="n", xaxt="n", xlab="", ylab=params.to.plot[param.index], bty="n", log=logs.to.use[param.index])
			y.range <- range(c(local.summary[,observed.column.index], values.true.x))
			if(params.to.plot[param.index]=="vh") {
				y.range <- c(-5, 75)	
			}
			if(params.to.plot[param.index]=="SE") {
				y.range <- c(-.2, 1.5)	
			}
			#y.range <- c(-500, 500)
			plot(x=c(1.5,(3*length(values.true.x)+1)), y=y.range, type="n", xaxt="n", xlab="", ylab=params.to.plot[param.index], bty="n", log=logs.to.use[param.index], main=paste(params.to.plot[param.index], " using ",names(summaries.by.estimator)[estimator.index], sep=""))
			for (values.true.index in sequence(length(values.true.x))) {
				base.position.set <- 2.7*values.true.index -1 + .9*seq(from=0, to=2, length.out=(length(values.ntip.nonhybrids) - 1) + length(values.ntip.nonhybrids)*length(values.ntip.hybrids))
				base.position.index = 0
				for (values.ntip.nonhybrid.index in sequence(length(values.ntip.nonhybrids))) {
					#base.position <- base.position + .4*(values.ntip.nonhybrid.index - 0.5)
					if(values.ntip.nonhybrid.index>1) {
						base.position.index<-base.position.index+1
					}
					for (values.ntip.hybrid.index in sequence(length(values.ntip.hybrids))) {
						base.position.index <- base.position.index+1
						#base.position <- base.position + .1*(values.ntip.hybrid.index - 1)
						#base.position.set <- append(base.position.set, base.position)
						base.position <- base.position.set[base.position.index]
						very.local.summary <- subset(local.summary, Nnonhybrid == values.ntip.nonhybrids[values.ntip.nonhybrid.index] & Nhybrid == values.ntip.hybrids[values.ntip.hybrid.index])
						true.column.index <- which(colnames(very.local.summary)==paste(params.to.plot[param.index],".true", sep=""))
						observed.column.index <- which(colnames(very.local.summary)==paste(params.to.plot[param.index], sep=""))
						very.local.summary <- very.local.summary[which(very.local.summary[, true.column.index] == values.true.x[values.true.index]),]
						values.y <- very.local.summary[, observed.column.index]
						#text(base.position, median(values.y), paste(values.ntip.nonhybrids[values.ntip.nonhybrid.index], values.ntip.hybrids[values.ntip.hybrid.index], values.true.x[values.true.index], sep="_"), cex=0.3)
						points(base.position, median(values.y), pch=20)
						quantile.results <- quantile(values.y, c(0.025, 0.975))
						print(quantile.results)
						print(median(values.y))
						arrows(x0=base.position, y0=max(quantile.results[1], min(y.range)), x1=base.position, y1=min(quantile.results[2], max(y.range)), length=0, col="black")
						text(base.position, min(quantile.results[2], max(y.range)), hybrid.labels[values.ntip.hybrid.index], pos=3, cex=0.8, col="darkgray")
						if(values.ntip.hybrid.index==2) {
							text(base.position, min(quantile.results[1], max(y.range)), nonhybrid.labels[values.ntip.nonhybrid.index], adj=c(.5,2), cex=0.8, col="darkgray")		
						}
						#lines(rep(base.position,2), median(values.y)+c(sd1, -sd1), col="red")
						
						print(paste(base.position, paste(values.ntip.nonhybrids[values.ntip.nonhybrid.index], values.ntip.hybrids[values.ntip.hybrid.index], values.true.x[values.true.index], sep="_"), median(values.y), length(values.y), min(values.y), max(values.y)))
					}
				}
				lines(c(min(base.position.set), max(base.position.set)), rep(values.true.x[values.true.index], 2), col="gray")	
			}
			print("done")
			#dev.off()
		}
	}
	dev.off()
}