library(ape)
phy<-read.nexus("ChaseClockResolved.tre")[[1]] #only want one of the trees
traits<-read.csv("KomoriEtAlTable2.csv", stringsAsFactors=FALSE)
library(utilitree)
library(rjson)
library(geiger)

#need to rescale tree to Nicotiana stem group age of 15.3 MY, from Clarkson et al. 2005, Long-term genome diploidization in allopolyploid Nicotiana sectionRepandae (Solanaceae)
patristic.distance.of.node <- cophenetic(phy)["Symonanthus_bancroftii", "Nicotiana_acaulis"]
phy$edge.length <- phy$edge.length*15.3 / (0.5*patristic.distance.of.node) #divided by two as it's distance down and back up


tree.names <- phy$tip.label
tree.names.resolved <- resolveNames(tree.names)
trait.names <- paste("Nicotiana ", traits$Species, sep="")
trait.names.resolved <- resolveNames(trait.names)
traits.focal<-traits$Mannitol
names(traits.focal)<-trait.names.resolved
phy$tip.label <- tree.names.resolved
pruned<-treedata(phy, traits.focal, sort=TRUE)
pruned.phy <- pruned$phy
pruned.traits <- pruned$data
write.tree(pruned.phy, file="Nicotiana.phy")
write.csv(pruned.traits[rev(sequence(dim(pruned.traits)[1])),], "Nicotiana.csv")