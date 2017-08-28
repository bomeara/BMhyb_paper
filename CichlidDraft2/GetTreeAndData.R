rm(list=ls())
library(ape)
library(geiger)
library(rfishbase)
#they had one as AY7403786 but it should be AY740386 (species matches)
cichlid.accessions<-c("DQ055016", "AY682545", "EF191082", "EF191108", "EF191105", "EF191121", "AF398229", "EF191107", "AF398226", "DQ055023", "EF191083", "EF191117", "EF191096", "AY740388", "DQ055024", "EF191122", "EF191119", "EF191120", "EF191103", "EF191104", "EF191089", "EF191090", "EF191091", "DQ055030", "EF191099", "EF191100", "EF191116", "EF191118", "EF191124", "EF191125", "EF191123", "EF191126", "EF191085", "EF191084", "EF191087", "EF191088", "AY740386", "DQ055019", "EF191093", "DQ055027", "EF191097", "EF191098", "EF191113", "EF191114", "EF191115", "EF191109", "EF191110", "EF191111", "EF191112", "EF191086", "DQ055032", "EF191101", "EF191102", "DQ055036", "DQ055037", "DQ055057", "DQ055034", "DQ055040", "EF191092", "DQ055018", "DQ055041", "DQ055051", "DQ055025", "DQ055038", "DQ055052", "EF191094", "EF191095") #From Koblmüller et al. 2007
bad<-c()
system("rm Cichlid.fasta")
seq<-read.GenBank(cichlid.accessions)
seq2<-seq
names(seq2) <- attr(seq, "species")
write.dna(seq2, file="Cichlid.fasta", format="fasta", nbcol=-1, append=FALSE, colsep = "")
system("mafft --auto Cichlid.fasta > Cichlid.aln.fasta")
seq3 <- read.FASTA("Cichlid.aln.fasta")
names(seq3) <- gsub("'", "", names(seq3))
seq4 <- as.character(seq3)
seq4<-seq3[!duplicated(names(seq3))]
write.dna(seq4, file="Cichlid.for.raxml.txt")
system("raxml -s Cichlid.for.raxml.txt -g Cichlid.Constraint.phy -f d -m GTRGAMMA -p 12345 -n Cichlid.GTRGAMMA")
constraint.tree <- read.tree("Cichlid.Constraint.phy")
raxml.tree <- read.tree("RAxML_bestTree.Cichlid.RAxML")
raxml.tree <- reorder(root(raxml.tree, outgroup="Julidochromis_ornatus"))
write.tree(raxml.tree, file="treePL.input.tre")

#Note: there were numerical issues with a fixed age of 4.19. So used age of 419, then correct brlen in R
treePL.correction = 100
cat(paste("treefile = treePL.input.tre\nsmooth = 100\nnumsites = ", dim(as.matrix(seq3))[2], "\nmrca = FROMPAPER Lamprologus_callipterus Neolamprologus_calliurus\nmin = FROMPAPER ", treePL.correction*round(mean(c(3.80, 4.58)),3), "\nmax = FROMPAPER ", treePL.correction*round(mean(c(3.80, 4.58)),3),"\noutfile = Cichlid.chronogram.tre\nthorough\nprime",sep=""), file="batch.treePL")
system("treePL batch.treePL")

#manually got the optimal parameters from first run. Ran again.

cat(paste("treefile = treePL.input.tre
smooth = 100
numsites = 1047
mrca = FROMPAPER Lamprologus_callipterus Neolamprologus_calliurus
min = FROMPAPER 419
max = FROMPAPER 419
outfile = Cichlid.chronogram.tre
opt = 5
moredetail
optad = 2
optcvad = 2
moredetailcvad
randomcv
"), file="batch.treePL")
system("treePL batch.treePL")



phy <- read.tree("Cichlid.chronogram.tre")
phy$edge.length <- phy$edge.length/treePL.correction

#need to rescale tree to use Kolbmuller et al calibration. They said: Within the lamprologines, the inferred hybridization partners are only distantly related, as the split between the mitochondrial lineages including N. brevis/calliurus and the lineage containing L. callipterus and N. fasciatus is esti- mated as 3.80 (± 1.33) – 4.58 (± 1.60) million years before presence
#We take the mean of this range for our point estimate
#patristic.distance.of.node <- cophenetic(phy)["Lamprologus_callipterus", "Neolamprologus_calliurus"]
#phy$edge.length <- phy$edge.length*mean(c(3.80, 4.58)) / (0.5*patristic.distance.of.node) #divided by two as it's distance down and back up


fish.data<- species(gsub("_", " ",(phy$tip.label)))
traits <- log(fish.data$Length)
names(traits) <- gsub(" ", "_", fish.data$sciname)
pruned<-treedata(phy, traits, sort=TRUE)



write.tree(pruned$phy, file="Cichlid.phy")
write.csv(pruned$data, "Cichlid.csv")
