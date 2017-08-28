
library(BMhyb)
source("GeneratePossibleSets.R")

CreateRun<-function(x, id, rep.id, possible.sets) {
	R.run.file=paste("ParamCombo", id, ".Rep", rep.id, ".R", sep="")
	batch.run.file=paste("ParamCombo", id, ".Rep", rep.id, ".batch", sep="")
	file.root = paste("ParamCombo", id, ".Rep", rep.id, sep="")

	cat("source('ComputeIndividualLikelihood.R')\n", file=R.run.file, append=FALSE)
	cat(paste("RunFromSet(possible.sets[",id, ",], id=as.character(", id, "), rep.id=", rep.id, ", possible.sets=possible.sets)\n", sep=""),  file=R.run.file, append=TRUE)
	cat(paste("#! /bin/sh\n\nR CMD BATCH ", R.run.file, sep=""), file=batch.run.file, append=TRUE)
    system(paste("chmod u+x ", batch.run.file, sep=""))
    cat(paste("executable=", batch.run.file, "\nuniverse=vanilla\noutput=results.output.", file.root, "\nerror=error.", file.root, "\ntransfer_input_files=GeneratePossibleSets.R,BMhyb_1.3.3.tar.gz,", R.run.file, ",", batch.run.file, ",ComputeIndividualLikelihood.R", "\nlog=log.",file.root,"\nnotification=never\nshould_transfer_files=YES\nwhen_to_transfer_output=ON_EXIT\nqueue 1\n", sep=""), file=paste("myjob.",file.root,".submit", sep=""))
    system(paste("/condor/condor-installed/bin/condor_submit myjob.",file.root,".submit", sep=""))
    Sys.sleep(2)
}



for(rep.id in sequence(nreps)) {
	id.converter <- sample.int(dim(possible.sets)[1])
	for (j in sequence(dim(possible.sets)[1])) {
		id <- id.converter[j] #this is done to reorder which params are run first
		try(CreateRun(possible.sets[id,], id=as.character(id), rep.id=rep.id, possible.sets=possible.sets))
	}
}
