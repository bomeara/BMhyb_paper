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
system("/usr/local/Cellar/raxml/8.2.11_1/bin/raxmlHPC-PTHREADS -s Cichlid.for.raxml.txt -g Cichlid.Constraint.phy -f d -m GTRGAMMA -p 12345 -n Cichlid.GTRGAMMA")
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


## Add validate names function

fishbase.names <- rfishbase::validate_names(gsub("_", " ", phy$tip.label))
fishbase.names <- unique(fishbase.names[!is.na(fishbase.names)])

#Pulled data from https://fms02.lsa.umich.edu/fmi/webd/ummz_fish for each genus. Saved as <genus>.tab

files.from.umich <- system("ls -1 *.tab", intern=TRUE)
umich.data <- read.delim(files.from.umich[1], header=FALSE, stringsAsFactors = FALSE)
for (i in 2:length(files.from.umich)) {
    umich.data <- rbind(umich.data, read.delim(files.from.umich[2], header=FALSE, stringsAsFactors = FALSE))   
}
umich.clean <- data.frame(genus=umich.data$V9, species=umich.data$V18, length=umich.data$V17, stringsAsFactors = FALSE)
umich.clean <- umich.clean[which(nchar(umich.clean$length)>0),]
umich.clean$length <- gsub(" mm SL", "", umich.clean$length)
umich.clean$length <- gsub(" mmS L", "", umich.clean$length) #human typos
umich.clean$min <- NA
umich.clean$max <- NA
for (i in sequence(nrow(umich.clean))) {
    if(grepl("-", umich.clean$length[i])) {
        lengths <- strsplit(umich.clean$length[i], "-")[[1]]   
        umich.clean$min[i] <- lengths[1]
        umich.clean$max[i] <- lengths[2]
    } else {
        if(grepl(",", umich.clean$length[i])) {
          lengths <- strsplit(umich.clean$length[i], ", ")[[1]]   
           umich.clean$min[i] <- lengths[1]
           umich.clean$max[i] <- lengths[2]
        } else {
            umich.clean$max[i] <- umich.clean$length[i]
            umich.clean$min[i] <- umich.clean$length[i]
        }
    }
}
umich.clean$min <- as.numeric(umich.clean$min)/10 #divided by 10 b/c fishbase is in cm, not mm. Note that this is also SL, whereas fishbase is TL, but still rough idea
umich.clean$max <- as.numeric(umich.clean$max)/10
umich.clean$range <- umich.clean$max - umich.clean$min
umich.clean$guess.at.sd <- sqrt((umich.clean$range))/2
umich.no.zero.range <- umich.clean[which(umich.clean$range>0),]
umich.for.prediction <- umich.no.zero.range[which(umich.no.zero.range$min>0.5*umich.no.zero.range$max),] #advice from borstein to toss those where there are probably subadults included
umich.for.prediction$guess <- (umich.for.prediction$guess.at.sd)
umich.for.prediction$log.max <- log(umich.for.prediction$max)
regression.result <- lm(guess ~ log.max, data=umich.for.prediction)
plot(umich.for.prediction$log.max, umich.for.prediction$guess)
abline(regression.result)


fish.data<- species(fishbase.names)
traits <- log(fish.data$Length)
names(traits) <- gsub(" ", "_", fish.data$sciname)
pruned<-treedata(phy, traits, sort=TRUE)

final.data <- pruned$data
data.to.fit <- (data.frame(log.max=final.data[,1]))
final.se <- (predict(regression.result, data.to.fit))
final.data <- cbind(final.data, final.se)




write.tree(pruned$phy, file="Cichlid.phy")
write.csv(final.data, "Cichlid.csv")
