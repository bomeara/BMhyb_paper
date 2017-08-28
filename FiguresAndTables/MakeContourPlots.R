#load("Users/bomeara/Dropbox/CollabJhwuengOMeara/CichlidExampleRevision/CichlidResults.RSave")
library(BMhyd)
load("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Cichlid/CichlidResults.RSave")


sims.cichlid <- data.frame(results.cichlid$sims[[4]])

load("/Users/bomeara/Dropbox/CollabJhwuengOMeara/ReReSubmission/Nicotiana/Nicotiana.RSave")

sims.nicotiana <- data.frame(results.nicotiana$sims[[4]])

all.params <- colnames(sims.nicotiana)[-1]

#colnames(sims) <- c("neglnL", colnames(sims)[-length(colnames(sims))]) #to temporarily fix data coming in wrong order

pdf(file="AllContours.pdf", width=15, height=15)
par(mfrow=c(5,5))

for (x.param in sequence(length(all.params))) {
	for (y.param in sequence(length(all.params))) {
		sims <- NA
		title.text <- NA
		if (x.param != y.param) {
			if(x.param<y.param) {
				sims <- sims.cichlid
				title.text <- "Cichlid"
			} else {
				sims <- sims.nicotiana	
				title.text <- "Nicotiana"
			}
			ContourFromAdaptiveSampling(sims, params.of.interest=c(all.params[x.param], all.params[y.param]))
			title(main=title.text)
		} else {
			plot(c(0,1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", type="n")	
		}
	}
}

#ContourFromAdaptiveSampling(sims)
dev.off()
system("open AllContours.pdf")

#For this, want to probably use Nicotiana data (fits better): take best model (or most complex), look at contour.
#Probably do on a large simulated dataset, too.


