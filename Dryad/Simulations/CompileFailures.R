rm(list=ls())
library(BMhyb)
library(ape)
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/Simulations")
source("GeneratePossibleSets.R")
counts.of.results <- rep(0, dim(possible.sets)[1])
possible.sets.renamed <- possible.sets
colnames(possible.sets.renamed) <- paste(colnames(possible.sets.renamed), ".true", sep="")


all.results <- data.frame()
current.file.count=0
all.files <- system(paste0("ls -1 *Try4*FAILURE.RData"),intern=TRUE)
for (j in sequence(dim(possible.sets)[1])) {
	id <- j
	files <- system(paste0("ls -1 *Try4*ParamCombo",id,".P*FAILURE.RData"),intern=TRUE, ignore.stderr=TRUE)
	for (file.index in sequence(length(files))) {
		rm(x)


		try(load(files[file.index]))
		if(!is.null(x)) {
			if(length(x)>0) {
			    #print(result4)
				file.split <- unlist(strsplit(strsplit(files[file.index],"_")[[1]],"\\."))
				current.file.count <- current.file.count + 1
				if(current.file.count%%5==0) {
					print(paste0("Done ", current.file.count, " of ", length(all.files), " files to do"))
				}
				try(all.results <- rbind(all.results, data.frame(x, stringsAsFactors=FALSE)))

				save(all.results , file="FailureResults.RData")
			}
		}
	}
}
