rm(list=ls())
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/NicotianaDraft4_July2018")
library(BMhyb)
library(geiger)
library(ape)
#library(optimx)
library(corpcor)
library(phylobase)
#options(warn=2)
options(warn=-1)
phy <- read.tree("nicotiana.phy")
data.raw <- read.csv("nicotiana.csv", stringsAsFactors=FALSE)
data.vector.raw <- as.numeric(data.raw[,2])
names(data.vector.raw) <- gsub(" ", "_", data.raw[,1])
flow<-data.frame(donor=c("Nicotiana_sylvestris"), recipient=c("Nicotiana_tabacum"), gamma=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_tabacum"))], stringsAsFactors=FALSE)


flow<-rbind(flow, data.frame(donor=c("Nicotiana_undulata"), recipient=c("Nicotiana_rustica"), gamma=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))
flow<-rbind(flow,data.frame(donor=c("Nicotiana_arentsii"), recipient=c("Nicotiana_rustica"), gamma=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_rustica"))], stringsAsFactors=FALSE))


flow<-rbind(flow, data.frame(donor=c("Nicotiana_wigandioides"), recipient=c("Nicotiana_arentsii"), gamma=c(0.5), time.from.root.donor=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))], time.from.root.recipient=max(branching.times(phy)) - phy$edge.length[which(phy$edge[,2]==which(phy$tip.label=="Nicotiana_arentsii"))],stringsAsFactors=FALSE))

pdf(file="NicotianaNetwork.pdf")
PlotNetwork(phy, flow)
dev.off()

#try(system("echo 'Making nicotiana tree done' | terminal-notifier -sound default"))

options(error = utils::recover)
n.points=50000
results.nicotiana.no.measurement.error <- NULL
#results.nicotiana.no.measurement.error3 <- BMhybGrid(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=0, n.points=n.points, attempt.deletion.fix=FALSE, get.se=TRUE, plot.se=TRUE, models=3)

try(results.nicotiana.no.measurement.error <- BMhybGrid(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=0, n.points=n.points, attempt.deletion.fix=FALSE, get.se=TRUE, plot.se=TRUE))
save(results.nicotiana.no.measurement.error, phy, flow, data.vector.raw, file="NicotianaResults.RData")

results.nicotiana.no.measurement.error.optimize <- NULL
try(results.nicotiana.no.measurement.error.optimize <- BMhyb(data.vector.raw, phy, flow, store.sims=TRUE, measurement.error=0, n.points=n.points, attempt.deletion.fix=FALSE, get.se=TRUE, plot.se=FALSE))
save(results.nicotiana.no.measurement.error, results.nicotiana.no.measurement.error.optimize, phy, flow, data.vector.raw, file="NicotianaResults.RData")

system("echo 'Zero measurement error nicotiana done' | terminal-notifier -sound default")
try(system("rm -rf WithFixedZeroError"))
system("mkdir WithFixedZeroError")
system("mv Model*.pdf WithFixedZeroError")
results.nicotiana.inferred.measurement.error <- NULL
try(results.nicotiana.inferred.measurement.error <- BMhybGrid(data.vector.raw, phy, flow, store.sims=TRUE, n.points=n.points, measurement.error=NULL, attempt.deletion.fix=FALSE, get.se=TRUE, plot.se=TRUE))
system("echo 'Estimating SE nicotiana done' | terminal-notifier -sound default")
save(results.nicotiana.inferred.measurement.error , results.nicotiana.no.measurement.error.optimize, results.nicotiana.no.measurement.error, phy, flow, data.vector.raw, file="NicotianaResults.RData")
try(system("rm -rf InferredError"))
try(system("mkdir InferredError"))
system("mv Model*.pdf InferredError")
all.results <- rbind(results.nicotiana.no.measurement.error$results, results.nicotiana.inferred.measurement.error$results)
all.results$deltaAICc <- all.results$AICc - min(all.results$AICc)
all.results$AkaikeWeight <- AkaikeWeight(all.results$deltaAICc)
all.sims <- rbind(results.nicotiana.no.measurement.error$sims, results.nicotiana.inferred.measurement.error$sims)
save(all.results, all.sims, file="NicotianaResults_AllCombined.RData")
save(list=ls(), file="NicotianaResults_DumpOfEverything.RData")
