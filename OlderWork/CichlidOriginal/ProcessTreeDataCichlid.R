library(ape)
library(rfishbase)
phy<-read.nexus("CichlidClockResolved.tre")
tip.names<-rep(NA, Ntip(phy))
for (i in sequence(Ntip(phy))) {
	tip.names[i] <- gsub("\\d", "", strsplit(phy$tip.label[i], "'")[[1]][2])
}
phy$tip.label <- tip.names

#traits<-read.csv("KomoriEtAlTable2.csv", stringsAsFactors=FALSE) #get traits from Rfishbase instead
fish.data.local<-loadCache() #Ask Brian this
fish.data.cichlids <- fish.data.local[findSpecies(phy$tip.label, fish.data.local)]
#updateCache()
traits <- getSize(fish.data.cichlids)#getSize is a function?

library(utilitree)
library(rjson)
library(geiger)

#need to rescale tree to use Kolbmuller et al calibration. They said: Within the lamprologines, the inferred hybridization partners are only distantly related, as the split between the mitochondrial lineages including N. brevis/calliurus and the lineage containing L. callipterus and N. fasciatus is esti- mated as 3.80 (± 1.33) – 4.58 (± 1.60) million years before presence
#We take the mean of this range for our point estimate
patristic.distance.of.node <- cophenetic(phy)["Lamprologus_callipterus", "Neolamprologus_calliurus"]
phy$edge.length <- phy$edge.length*mean(c(3.80, 4.58)) / (0.5*patristic.distance.of.node) #divided by two as it's distance down and back up


#tree.names <- phy$tip.label
#tree.names.resolved <- resolveNames(tree.names)
#trait.names.resolved <- resolveNames(names(traits))
#names(traits.focal)<-trait.names.resolved
#phy$tip.label <- tree.names.resolved

pruned<-treedata(drop.tip(phy, tip=which(duplicated(phy$tip.label))), traits, sort=TRUE)
pruned.phy <- pruned$phy
pruned.phy <- drop.tip(pruned.phy, tip=which(duplicated(pruned.phy$tip.label)))
pruned.traits <- pruned$data

write.tree(pruned.phy, file="Cichlid.phy")
write.csv(pruned.traits[rev(sequence(dim(pruned.traits)[1])),], "Cichlid.csv")
