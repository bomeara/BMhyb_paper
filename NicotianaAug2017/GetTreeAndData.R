rm(list=ls())
library(ape)
library(geiger)
library(taxize)
setwd("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Nicotiana")
accessions <- read.delim(file="Accessions.txt", header=FALSE)
seq <- read.GenBank(accessions[,1])
seq2<-seq
names(seq2) <- attr(seq, "species")
write.dna(seq2, file="nicotiana.fasta", format="fasta", nbcol=-1, append=FALSE, colsep = "")
system("mafft --auto nicotiana.fasta > nicotiana.aln.fasta")
seq3 <- read.FASTA("nicotiana.aln.fasta")
names(seq3) <- gsub("'", "", names(seq3))
seq4 <- as.character(seq3)
seq4<-seq3[!duplicated(names(seq3))]
write.dna(seq4, file="nicotiana.for.raxml.txt")
system("raxml -s nicotiana.for.raxml.txt -f d -m GTRGAMMA -p 12345 -n nicotiana.GTRGAMMA")
raxml.tree <- read.tree("RAxML_bestTree.nicotiana.GTRGAMMA")
raxml.tree <- reorder(root(raxml.tree, outgroup="Petunia_axillaris", resolve.root=TRUE))
write.tree(raxml.tree, file="treePL.input.tre")
treePL.correction = 100
#The genera Nicotiana and Symonanthus (tribe Anthocercideae) are sister groups (Clarkson et al., 2004) that split approx. 15.3 Myr ago 
cat(paste("treefile = treePL.input.tre\nsmooth = 100\nnumsites = ", dim(as.matrix(seq3))[2], "\nmrca = FROMPAPER Symonanthus_bancroftii Nicotiana_undulata\nmin = FROMPAPER ", treePL.correction*15.3, "\nmax = FROMPAPER ", treePL.correction*15.3,"\noutfile = nicotiana.chronogram.tre\nthorough\nprime",sep=""), file="batch.treePL")
system("treePL batch.treePL")
system("rm nicotiana.chronogram.tre")
#manually got the optimal parameters from first run. Ran again.


cat(paste("treefile = treePL.input.tre\nsmooth = 100\nnumsites = ", dim(as.matrix(seq3))[2], "\nmrca = FROMPAPER Symonanthus_bancroftii Nicotiana_undulata\nmin = FROMPAPER ", treePL.correction*15.3, "\nmax = FROMPAPER ", treePL.correction*15.3,"\noutfile = nicotiana.chronogram.tre
opt = 5
optad = 2
moredetailad
optcvad = 2
moredetailcvad
randomcv
",sep=""), file="batch.treePL")

system("treePL batch.treePL")



phy <- read.tree("nicotiana.chronogram.tre")
phy$edge.length <- phy$edge.length/treePL.correction

new.names <- rep(NA, length(phy$tip.label))
for (i in sequence(length(phy$tip.label))) {
	new.names[i] <- gnr_resolve(phy$tip.label[i], best_match_only=TRUE, data_source_ids= gnr_datasources()$id[which(gnr_datasources()$title=="NCBI")])$matched_name

}

phy$tip.label <- gsub(" ", "_", new.names)

traits<-read.csv("KomoriEtAlTable2.csv", stringsAsFactors=FALSE)
trait.names <- paste("Nicotiana ", traits$Species, sep="")
trait.names.resolved <- rep(NA, length(trait.names))
for (i in sequence(length(trait.names))) {
	trait.names.resolved[i] <- gsub(" ", "_",gnr_resolve(trait.names[i], best_match_only=TRUE, data_source_ids= gnr_datasources()$id[which(gnr_datasources()$title=="NCBI")])$matched_name)

}

traits.focal <- log(traits$Mannitol)
names(traits.focal) <- trait.names.resolved


pruned<-treedata(phy, traits.focal, sort=TRUE)



write.tree(pruned$phy, file="nicotiana.phy")
write.csv(pruned$data, "nicotiana.csv")
