rm(list=ls())
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidDraft4_July2018/")
library(BMhyb)
library(geiger)
library(ape)
#library(optimx)
library(corpcor)
library(phylobase)
library(phytools)
options(warn=1)
phy <- read.tree("Cichlid.phy")
phy <- root(phy, "Julidochromis_ornatus", resolve.root=TRUE)
data.raw <- read.csv("Cichlid.csv", stringsAsFactors=FALSE)
data.vector.raw <- as.numeric(data.raw[,2])
names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
error.vector.raw <- as.numeric(data.raw[,3])
names(error.vector.raw) <- gsub(" ", "_", data.raw[,1])
flow<-data.frame(donor=c("Lamprologus_callipterus", "Lamprologus_meleagris", "Lamprologus_ocellatus", "Neolamprologus_wauthioni", "Lamprologus_speciosus"), recipient=rep("Lamprologus_meleagris", 5), gamma=rep(0.5,5), time.from.root.donor=rep(nodeheight(phy, node=findMRCA(phy, tips=c("Lamprologus_callipterus", "Lamprologus_meleagris", "Lamprologus_ocellatus", "Neolamprologus_wauthioni", "Lamprologus_speciosus"), type="node")), 5), time.from.root.recipient=rep(max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_meleagris"))], 5), stringsAsFactors=FALSE)
flow<-rbind(flow, data.frame(donor=c("Lamprologus_callipterus"), recipient="Neolamprologus_wauthioni", gamma=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_wauthioni"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_wauthioni"))], stringsAsFactors=FALSE))
flow<-rbind(flow, data.frame(donor=c("Neolamprologus_brevis"), recipient="Lamprologus_speciosus", gamma=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_speciosus"))],time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_speciosus"))], stringsAsFactors=FALSE))
flow<-rbind(flow, data.frame(donor=c("Lamprologus_callipterus"), recipient="Neolamprologus_fasciatus", gamma=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_fasciatus"))],time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_fasciatus"))], stringsAsFactors=FALSE))

ingroup<-phy$tip.label[-(match(c("Julidochromis_ornatus", "Telmatochromis_vittatus", "Variabilichromis_moorii"), phy$tip.label))]

#repeat ingroup as with L. melaeagris
flow<-rbind(flow, data.frame(donor=ingroup, recipient=rep("Neolamprologus_multifasciatus", length(ingroup)), gamma=rep(0.5, length(ingroup)), time.from.root.donor=rep(nodeheight(phy, node=findMRCA(phy, ingroup, type="node")), length(ingroup)), time.from.root.recipient=rep(max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_multifasciatus"))], length(ingroup)), stringsAsFactors=FALSE))
pdf(file="CichlidNetwork.pdf")
PlotNetwork(ladderize(phy,right=FALSE), flow, head.length=0.1, arrow.width=3, try.rotations=TRUE)
dev.off()
try(system("echo 'Making cichlid tree done' | terminal-notifier -sound default"))

options(error = utils::recover)
results.cichlid.no.measurement.error <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=0, badval.if.not.positive.definite=FALSE, attempt.deletion.fix=FALSE)
save(results.cichlid.no.measurement.error, phy, flow, data.vector.raw, file="CichlidResults_zeroerror.RData")
try(system("echo 'Measurement error cichlid done' | terminal-notifier -sound default"))
try(system("rm -rf WithFixedZeroError"))
system("mkdir WithFixedZeroError")
system("mv Model*.pdf WithFixedZeroError")
system("mv CichlidResults_zeroerror.RData WithFixedZeroError")

try(system("rm -rf FixedActualError"))
try(system("mkdir FixedActualError"))
results.cichlid.known_nonzero.measurement.error <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=error.vector.raw, badval.if.not.positive.definite=FALSE, attempt.deletion.fix=FALSE)
save(results.cichlid.known_nonzero.measurement.error, phy, flow, data.vector.raw, error.vector.raw, file="CichlidResults_fixednonzeroerror.RData")
try(system("echo 'Measurement error cichlid done' | terminal-notifier -sound default"))
system("mv Model*.pdf FixedActualError")
system("mv CichlidResults_fixednonzeroerror.RData FixedActualError")

try(system("rm -rf InferredError"))
try(system("mkdir InferredError"))

results.cichlid.inferrederror <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=NULL, badval.if.not.positive.definite=FALSE, attempt.deletion.fix=FALSE)
try(system("echo 'Estimating SE cichlid done' | terminal-notifier -sound default"))

save(results.cichlid.inferrederror, phy, flow, data.vector.raw, file="CichlidResults_inferrederror.RData")

try(system("echo 'Inferred Measurement error cichlid done' | terminal-notifier -sound default"))
system("mv Model*.pdf InferredError")
system("mv CichlidResults_inferrederror.RData InferredError")

results.cichlid.no.measurement.error$results$error <- "Zero"
results.cichlid.known_nonzero.measurement.error$results$error <- "Empirical"
results.cichlid.inferrederror$results$error <- "Inferred"
results.cichlid.no.measurement.error$sims$error <- "Zero"
results.cichlid.known_nonzero.measurement.error$sims$error <- "Empirical"
results.cichlid.inferrederror$sims$error <- "Inferred"
all.results <- rbind(results.cichlid.no.measurement.error$results, results.cichlid.known_nonzero.measurement.error$results, results.cichlid.inferrederror$results)
all.results$deltaAICc <- all.results$AICc - min(all.results$AICc)
all.results$AkaikeWeight <- AkaikeWeight(all.results$deltaAICc)
all.sims <- rbind(results.cichlid.no.measurement.error$sims, results.cichlid.known_nonzero.measurement.error$sims, results.cichlid.inferrederror$sims)
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidDraft4_July2018/")
save(all.results, all.sims, file="CichlidResults_AllCombined.RData")
save(list=ls(), file="CichlidResults_DumpOfEverything.RData")

#phy2<-drop.tip(phy, "Nicotiana_alata")
#pruned2<-treedata(phy2, data.log)
# data2.log <- pruned2$data[,1]
# names(data2.log) <-rownames(pruned2$data)
# flow2<-data.frame(donor=c("Nicotiana_sylvestris", "Nicotiana_paniculata"), recipient=c("Nicotiana_tabacum", "Nicotiana_rustica"), gamma=c(0.5, 0.5), time.from.root=c(0.042861 - 0.000524, 0.042861 - 0.001893), stringsAsFactors=FALSE)

# results <- BMhyd(data2.log, pruned2$phy, flow2, models=c(1))
