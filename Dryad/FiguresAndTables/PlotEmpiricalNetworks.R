setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidAug2017")
library(BMhyb)
library(geiger)
library(ape)
#library(optimx)
library(corpcor)
library(phylobase)
library(phytools)

load("CichlidResults_DumpOfEverything.RData")
# phy <- read.tree("Cichlid.phy")
# data.raw <- read.csv("Cichlid.csv", stringsAsFactors=FALSE)
# data.vector.raw <- as.numeric(data.raw[,2])
# names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
# data.log <- log(data.vector.raw)
# flow<-data.frame(donor=c("Lamprologus_callipterus", "Lamprologus_meleagris", "Lamprologus_ocellatus", "Neolamprologus_wauthioni", "Lamprologus_speciosus"), recipient=rep("Lamprologus_meleagris", 5), m=rep(0.5,5), time.from.root.donor=rep(nodeheight(phy, node=findMRCA(phy, tips=c("Lamprologus_callipterus", "Lamprologus_meleagris", "Lamprologus_ocellatus", "Neolamprologus_wauthioni", "Lamprologus_speciosus"), type="node")), 5), time.from.root.recipient=rep(max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_meleagris"))], 5), stringsAsFactors=FALSE)
# flow<-rbind(flow, data.frame(donor=c("Lamprologus_callipterus"), recipient="Neolamprologus_wauthioni", m=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_wauthioni"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_wauthioni"))], stringsAsFactors=FALSE))
# flow<-rbind(flow, data.frame(donor=c("Neolamprologus_brevis"), recipient="Lamprologus_speciosus", m=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_speciosus"))],time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Lamprologus_speciosus"))], stringsAsFactors=FALSE))
# flow<-rbind(flow, data.frame(donor=c("Lamprologus_callipterus"), recipient="Neolamprologus_fasciatus", m=0.5, time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_fasciatus"))],time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_fasciatus"))], stringsAsFactors=FALSE))
# 
# ingroup<-phy$tip.label[-(match(c("Julidochromis_ornatus", "Telmatochromis_vittatus", "Variabilichromis_moorii"), phy$tip.label))]
# 
# #repeat ingroup as with L. melaeagris
# flow<-rbind(flow, data.frame(donor=ingroup, recipient=rep("Neolamprologus_multifasciatus", length(ingroup)), m=rep(0.5, length(ingroup)), time.from.root.donor=rep(nodeheight(phy, node=findMRCA(phy, ingroup, type="node")), length(ingroup)), time.from.root.recipient=rep(max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Neolamprologus_multifasciatus"))], length(ingroup)), stringsAsFactors=FALSE))
# phy$tip.label <- gsub("_", " ", phy$tip.label)
# flow$donor <- gsub("_", " ", flow$donor)
# flow$recipient <- gsub("_", " ", flow$recipient)
# 



phy.cichlid <- phy
flow.cichlid <- flow


setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/NicotianaAug2017")

load("NicotianaResults_DumpOfEverything.RData")

# phy <- read.tree("Nicotiana.phy")
# data.raw <- read.csv("Nicotiana.csv", stringsAsFactors=FALSE)
# data.vector.raw <- as.numeric(data.raw[,2])
# names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
# data.log <- log(data.vector.raw)
# flow<-data.frame(donor=c("Nicotiana_sylvestris"), recipient=c("Nicotiana_tabacum"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], stringsAsFactors=FALSE)
# 
# 
# flow<-rbind(flow, data.frame(donor=c("Nicotiana_undulata"), recipient=c("Nicotiana_rustica"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))
# flow<-rbind(flow,data.frame(donor=c("Nicotiana_arentsii"), recipient=c("Nicotiana_rustica"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))
# 
# 
# flow<-rbind(flow, data.frame(donor=c("Nicotiana_wigandioides"), recipient=c("Nicotiana_arentsii"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))],stringsAsFactors=FALSE))
# 
# phy$tip.label <- gsub("_", " ",gsub("Nicotiana", "N.", phy$tip.label))
# flow$donor <- gsub("_", " ",gsub("Nicotiana", "N.", flow$donor))
# flow$recipient <- gsub("_", " ",gsub("Nicotiana", "N.", flow$recipient))

phy.nicotiana <- phy
flow.nicotiana <- flow

setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/FiguresAndTablesAug2017")

pdf(file="EmpiricalNetworks.pdf", width=10, height=5)
par(mfcol=c(1,2))
cex.factor = 0.4
padding.factor = 2
length.factor = 0.07
PlotNetwork(phy.cichlid, flow.cichlid, name.padding=1.5, cex=0.5, head.length=length.factor, col.donor="black", main="Cichlid", try.rotations = TRUE)
PlotNetwork(phy.nicotiana, flow.nicotiana, name.padding=1.3, cex=0.35, head.length=length.factor, col.donor="black", main="Nicotiana")
dev.off()
