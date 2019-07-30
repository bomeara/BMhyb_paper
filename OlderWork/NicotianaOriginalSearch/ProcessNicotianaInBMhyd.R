rm(list=ls())
library(BMhyb)
library(geiger)
library(ape)
#library(optimx)
library(corpcor)
library(phylobase)
phy <- read.tree("Nicotiana.phy")
data.raw <- read.csv("Nicotiana.csv", stringsAsFactors=FALSE)
data.vector.raw <- as.numeric(data.raw[,2])
names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
data.log <- log(data.vector.raw)
flow<-data.frame(donor=c("Nicotiana_sylvestris"), recipient=c("Nicotiana_tabacum"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], stringsAsFactors=FALSE)


flow<-rbind(flow, data.frame(donor=c("Nicotiana_undulata"), recipient=c("Nicotiana_rustica"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))
flow<-rbind(flow,data.frame(donor=c("Nicotiana_arentsii"), recipient=c("Nicotiana_rustica"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))


flow<-rbind(flow, data.frame(donor=c("Nicotiana_wigandioides"), recipient=c("Nicotiana_arentsii"), m=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))],stringsAsFactors=FALSE))

pdf(file="NicotianaNetwork.pdf")
PlotNetwork(phy, flow)
dev.off()

results.nicotiana <- BMhyb(data.log, phy, flow, models=c(1:4), store.sims=TRUE, get.se=TRUE)
save(results.nicotiana, data.log, phy, flow, file="Nicotiana.RSave")
#library(subselect)
#phy2<-drop.tip(phy, trim.matrix(vcv(phy))$names.discarded)
#pruned2<-treedata(phy2, data.log)
# data2.log <- pruned2$data[,1]
# names(data2.log) <-rownames(pruned2$data)
# flow2<-data.frame(donor=c("Nicotiana_sylvestris", "Nicotiana_paniculata"), recipient=c("Nicotiana_tabacum", "Nicotiana_rustica"), m=c(0.5, 0.5), time.from.root=c(0.042861 - 0.000524, 0.042861 - 0.001893), stringsAsFactors=FALSE)

## results2 <- BMhyd(data2.log, pruned2$phy, flow, models=c(2), get.se=FALSE)
 #x<-c(2.141756e+00  ,1.681108e+01,  1.273650e+01 , 0, 1.428251e+00 )
#		free.parameters<-rep(TRUE, 5)
 #		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")

# CalculateLikelihood(x, data.log, phy, flow, actual.params=free.parameters, proportion.mix.with.diag=1)
