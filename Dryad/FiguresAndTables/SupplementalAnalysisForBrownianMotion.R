library(BMhyb)
library(TreeSim)
library(ape)
library(geiger)
is.good <- TRUE
tries <- 0
while(is.good) {
    phy.original <- sim.bd.taxa(n=round(runif(1, 30, 100)), numbsim = 1, lambda=.1, mu=0.05, frac=0.5, complete=TRUE)[[1]]
    phy <- phy.original
     phy$edge.length[unique(sample.int(Nnode(phy), size=round(runif(1, 1, Ntip(phy)))))] <- 1e-16
    is.good <- IsPositiveDefinite(vcv(phy))
    tries <- tries+1
}
plot.phylo(phy, show.tip.label=FALSE)
axisPhylo()
traits <- geiger::sim.char(phy, par=1, model="BM")

#This should return an error
phytools::ratebytree(c(phy, phy.original), list(traits[,,1]))

#note that Geiger won't return an error when analyzing likelihood, as it adopts Ho and Ané's approach that doesn't require matrix inversion. This is clever, but the three point condition does not hold for networks.