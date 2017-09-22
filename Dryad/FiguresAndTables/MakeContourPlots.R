rm(list=ls())
library(BMhyb)
library(ape)


CalculateAllWeights <- function(x, phy) {
    x$AICc <- AICc(Ntip(phy),k=x$K, x$NegLogL)
    x$deltaAICc <- x$AICc - min(x$AICc)
    x$AkaikeWeight <- AkaikeWeight(x$deltaAICc)
    return(x)
}


load("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidAug2017/CichlidResults_AllCombined.RData")
cichlid.phy <- read.tree("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/CichlidAug2017/cichlid.phy")
nicotiana.phy <- read.tree("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/NicotianaAug2017/nicotiana.phy")


sims.cichlid <- data.frame()
SE.names <- c("Zero", "Empirical", "Inferred")
for (i in 1:12) { #hack to deal with sims that don't align
    local.df <- all.sims[[i]]
    K <- ncol(local.df) - 1
    local.df$K <- K
    if(!any(grepl("bt", names(local.df)))) {
        local.df$bt <- 1   
    }
    if(!any(grepl("vh", names(local.df)))) {
        local.df$vh <- 0
    }
    if(!any(grepl("SE", names(local.df)))) {
        local.df$SE <- SE.names[1+i%%3]
    }
    if(i==1) {
        sims.cichlid <- local.df   
    } else {
        sims.cichlid <- rbind(sims.cichlid, local.df)
    }
}

names(sims.cichlid)[grepl("neglnL", names(sims.cichlid))] <- "NegLogL"

sims.cichlid <- CalculateAllWeights(sims.cichlid, cichlid.phy)


load("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/NicotianaAug2017/NicotianaResults_AllCombined.RData")

sims.nicotiana <- all.sims

sims.nicotiana <- CalculateAllWeights(sims.nicotiana, nicotiana.phy)

all.params <- c("sigma.sq", "mu", "vh", "bt")

#colnames(sims) <- c("neglnL", colnames(sims)[-length(colnames(sims))]) #to temporarily fix data coming in wrong order

setwd("/Users/bomeara/Documents/MyDocuments/GitClones/BMhyb_paper/FiguresAndTablesAug2017")
pdf(file="AllContours.pdf", width=15, height=15)
par(mfrow=c(length(all.params), length(all.params)))

sims.cichlid$neglnL <- sims.cichlid$NegLogL
sims.nicotiana$neglnL <- sims.nicotiana$NegLogL


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


