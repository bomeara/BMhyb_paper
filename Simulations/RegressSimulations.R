library(MuMIn)
rm(list=ls())
load("SimulationResults.RData")
average.results$bt.deviation <- abs(average.results$bt.true - average.results$bt) / ifelse(average.results$bt.true==0, 1, average.results$bt.true)
average.results$sigma.sq.deviation <- abs(average.results$sigma.sq.true - average.results$sigma.sq) / ifelse(average.results$sigma.sq.true==0, 1, average.results$sigma.sq.true)
average.results$vh.deviation <- abs(average.results$vh.true - average.results$vh) / ifelse(average.results$vh.true==0, 1, average.results$vh.true)
average.results$mu.deviation <- abs(average.results$mu.true - average.results$mu) / ifelse(average.results$mu.true==0, 1, average.results$mu.true)

options(na.action = "na.fail")
bt.model <- lm(bt.deviation ~ ntax.nonhybrid.true + ntax.hybrid.true + flow.proportion.true + origin.type.true + tree.height.true, data=average.results)
bt.dredge <- dredge(bt.model)

quantile(average.results$bt.deviation, c(0.025, .5, 0.975))

quantile(average.results$vh.deviation, c(0.025, .5, 0.975))

quantile(average.results$sigma.sq.deviation, c(0.025, .5, 0.975))

quantile(average.results$mu.deviation, c(0.025, .5, 0.975))



quantile(average.results$bt.deviation[which(average.results$ntax.nonhybrid.true==30 & average.results$ntax.hybrid.true==1)], c(0.025, .5, 0.975))

quantile(average.results$bt.deviation[which(average.results$ntax.nonhybrid.true==30 & average.results$ntax.hybrid.true==10)], c(0.025, .5, 0.975))

quantile(average.results$bt.deviation[which(average.results$ntax.nonhybrid.true==100 & average.results$ntax.hybrid.true==1)], c(0.025, .5, 0.975))

quantile(average.results$bt.deviation[which(average.results$ntax.nonhybrid.true==100 & average.results$ntax.hybrid.true==10)], c(0.025, .5, 0.975))

quantile(average.results$bt.deviation[which(average.results$ntax.nonhybrid.true==30 & average.results$ntax.hybrid.true==1)], c(0.025, .5, 0.975))

quantile(average.results$bt.deviation[which(average.results$ntax.nonhybrid.true==30 & average.results$ntax.hybrid.true==10)], c(0.025, .5, 0.975))

quantile(average.results$bt[which(average.results$ntax.nonhybrid.true==100 & average.results$ntax.hybrid.true==10 & average.results$bt.true==1)], c(0.025, .5, 0.975))
