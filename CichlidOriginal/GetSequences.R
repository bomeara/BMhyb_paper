library(ape)
#they had one as AY7403786 but it should be AY740386 (species matches)
cichlid.accessions<-c("DQ055016", "AY682545", "EF191082", "EF191108", "EF191105", "EF191121", "AF398229", "EF191107", "AF398226", "DQ055023", "EF191083", "EF191117", "EF191096", "AY740388", "DQ055024", "EF191122", "EF191119", "EF191120", "EF191103", "EF191104", "EF191089", "EF191090", "EF191091", "DQ055030", "EF191099", "EF191100", "EF191116", "EF191118", "EF191124", "EF191125", "EF191123", "EF191126", "EF191085", "EF191084", "EF191087", "EF191088", "AY740386", "DQ055019", "EF191093", "DQ055027", "EF191097", "EF191098", "EF191113", "EF191114", "EF191115", "EF191109", "EF191110", "EF191111", "EF191112", "EF191086", "DQ055032", "EF191101", "EF191102", "DQ055036", "DQ055037", "DQ055057", "DQ055034", "DQ055040", "EF191092", "DQ055018", "DQ055041", "DQ055051", "DQ055025", "DQ055038", "DQ055052", "EF191094", "EF191095") #From Koblmüller et al. 2007
bad<-c()
system("rm Cichlid.fasta")
seq<-read.GenBank(cichlid.accessions)
seq2<-seq
names(seq2) <- attr(seq, "species")
write.dna(seq2, file="Cichlid.fasta", format="fasta", nbcol=-1, append=FALSE, colsep = "")
system("mafft --auto Cichlid.fasta > Cichlid.aln.fasta")
