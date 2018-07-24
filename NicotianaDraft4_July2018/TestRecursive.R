rm(list=ls())
library(BMhyb)
library(geiger)
library(ape)
#library(optimx)
library(corpcor)
library(phylobase)
phy <- read.tree("nicotiana.phy")
data.raw <- read.csv("nicotiana.csv", stringsAsFactors=FALSE)
data.vector.raw <- as.numeric(data.raw[,2])
names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
flow<-data.frame(donor=c("Nicotiana_sylvestris"), recipient=c("Nicotiana_tabacum"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], stringsAsFactors=FALSE)


flow<-rbind(flow, data.frame(donor=c("Nicotiana_undulata"), recipient=c("Nicotiana_rustica"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))
flow<-rbind(flow,data.frame(donor=c("Nicotiana_arentsii"), recipient=c("Nicotiana_rustica"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))


flow<-rbind(flow, data.frame(donor=c("Nicotiana_wigandioides"), recipient=c("Nicotiana_arentsii"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))],stringsAsFactors=FALSE))
options(error = utils::recover)
n.points=200000
test <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE,  n.points=n.points, models=2)
