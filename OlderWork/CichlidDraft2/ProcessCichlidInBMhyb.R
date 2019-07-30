rm(list=ls())
library(BMhyb)
library(geiger)
library(ape)
#library(optimx)
library(corpcor)
library(phylobase)
library(phytools)
phy <- read.tree("Cichlid.phy")
phy <- root(phy, "Julidochromis_ornatus", resolve.root=TRUE)
data.raw <- read.csv("Cichlid.csv", stringsAsFactors=FALSE)
data.vector.raw <- as.numeric(data.raw[,2])
names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
flow<-data.frame(donor=c("Lamprologus_callipterus", "Lamprologus_meleagris", "Lamprologus_ocellatus", "Neolamprologus_wauthioni", "Lamprologus_speciosus"), recipient=rep("Lamprologus_meleagris", 5), m=rep(0.5,5), time.from.root.donor=rep(nodeheight(phy, node=findMRCA(phy, tips=c("Lamprologus_callipterus", "Lamprologus_meleagris", "Lamprologus_ocellatus", "Neolamprologus_wauthioni", "Lamprologus_speciosus"), type="node")), 5), time.from.root.recipient=rep(max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_meleagris"))], 5), stringsAsFactors=FALSE)
flow<-rbind(flow, data.frame(donor=c("Lamprologus_callipterus"), recipient="Neolamprologus_wauthioni", m=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_wauthioni"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_wauthioni"))], stringsAsFactors=FALSE))
flow<-rbind(flow, data.frame(donor=c("Neolamprologus_brevis"), recipient="Lamprologus_speciosus", m=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_speciosus"))],time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_speciosus"))], stringsAsFactors=FALSE))
flow<-rbind(flow, data.frame(donor=c("Lamprologus_callipterus"), recipient="Neolamprologus_fasciatus", m=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_fasciatus"))],time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_fasciatus"))], stringsAsFactors=FALSE))

ingroup<-phy$tip.label[-(match(c("Julidochromis_ornatus", "Telmatochromis_vittatus", "Variabilichromis_moorii"), phy$tip.label))]

#repeat ingroup as with L. melaeagris
flow<-rbind(flow, data.frame(donor=ingroup, recipient=rep("Neolamprologus_multifasciatus", length(ingroup)), m=rep(0.5, length(ingroup)), time.from.root.donor=rep(nodeheight(phy, node=findMRCA(phy, ingroup, type="node")), length(ingroup)), time.from.root.recipient=rep(max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_multifasciatus"))], length(ingroup)), stringsAsFactors=FALSE))
pdf(file="CichlidNetwork.pdf")
PlotNetwork(phy, flow)
dev.off()
system("echo 'Making cichlid tree done' | terminal-notifier -sound default")

options(error = utils::recover)
results.cichlid.no.measurement.error <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=0)
save(results.cichlid.no.measurement.error, phy, flow, data.vector.raw, file="CichlidResults.RData")
system("echo 'Measurement error cichlid done' | terminal-notifier -sound default")
system("mv *.pdf WithFixedZeroError")
results.cichlid <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE)
system("echo 'Estimating SE cichlid done' | terminal-notifier -sound default")

save(results.cichlid, results.cichlid.no.measurement.error, phy, flow, data.vector.raw, file="CichlidResults.RData")
#phy2<-drop.tip(phy, "Nicotiana_alata")
#pruned2<-treedata(phy2, data.log)
# data2.log <- pruned2$data[,1]
# names(data2.log) <-rownames(pruned2$data)
# flow2<-data.frame(donor=c("Nicotiana_sylvestris", "Nicotiana_paniculata"), recipient=c("Nicotiana_tabacum", "Nicotiana_rustica"), m=c(0.5, 0.5), time.from.root=c(0.042861 - 0.000524, 0.042861 - 0.001893), stringsAsFactors=FALSE)

# results <- BMhyd(data2.log, pruned2$phy, flow2, models=c(1))
