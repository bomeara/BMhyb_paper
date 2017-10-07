rm(list=ls())
library(ape)
library(BMhyb) # Version 1.5.1, published 2017-09-11 on the CRAN
# (see session infos at the botom of the script)

## Parameters
sigma2 <- 1
gamma <- 0.5

## Function to create the network presented in the article
create_paper_network <- function(gamma, t1, t2, t3){
  phy <- read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
  network <- list(phy = phy,
                  flow = data.frame(donor = "X",
                                    recipient = "R",
                                    gamma = gamma,
                                    time.from.root.donor = t1,
                                    time.from.root.recipient = t1 + t2))
  network$flow$donor <- as.character(network$flow$donor)
  network$flow$recipient <- as.character(network$flow$recipient)
  return(network)
}



###############################################################################
## Network with several tips descendants from on hybrid
###############################################################################
# This illustrate an other potential problem, when several tips are descending
# from a single hybridization event.

# (Going through the code, I noticed that "V[recipient1, recipient2]" does not
# seem to be modified, which is coherent with what I find here.)
#######################################################
##CASE 2: (same donor) flow came in "after" speciation
#######################################################
## Underlying tree
t1 <- 0.3; t2 <- 0.4; t3 <- 0.3;
phy <- read.tree(text = paste0("((X:", t2+t3, ",R1:", t2+t3, "):", t1, ",R2:", t1 + t2 + t3, ");"))
plot(phy)

## Network 
don_recp <- expand.grid(c("X"), c("R1", "R2"))
network <- list(phy = phy,
                flow = data.frame(donor = don_recp[,1],
                                  recipient = don_recp[,2],
                                  gamma = rep(gamma, 2),
                                  time.from.root.donor = rep(t1+t2, 2),
                                  time.from.root.recipient = rep(t1+t2, 2)))
network$flow$donor <- as.character(network$flow$donor)
network$flow$recipient <- as.character(network$flow$recipient)
## Plot
PlotNetwork(network$phy, network$flow)
axis(1, at = c(0, t1, t1+t2, t1+t2+t3), labels = c("0", "t1", "t1+t2", "t1+t2+t3"))
## VCV network
GetVModified(x, network$phy, network$flow, actual.params)
#      R    Y    X
# R 0.85 0.70 0.15
# Y 0.70 0.85 0.15
# X 0.15 0.15 1.00
#V.modified[recipient.index, donor.index] <- (1 - gamma) * V.original[recipient.index, donor.index] + (gamma) * (flow$time.from.root.recipient[flow.index]) * sigma.sq


## VCV underlying tree
vcv.phylo(network$phy, model="Brownian")
#     R   Y X
# R 1.0 0.7 0
# Y 0.7 1.0 0
# X 0.0 0.0 1
bt <- 1
vh <- 0
sigma.sq <- 1
mu <- 0
SE <- 0
phy<-network$phy
flow<-network$flow
phy
flow
times.original <- vcv.phylo(phy, model = "Brownian")
times.original
sigma.sq<-1
V.original <- sigma.sq * times.original
V.modified <- V.original

for (flow.index in sequence(dim(flow)[1])) {
  #flow.index<-1
  recipient.index <- which(rownames(V.modified) == flow$recipient[flow.index])
  recipient.index
  gamma <- flow$gamma[flow.index]
  recipient.index
  gamma
  rownames(V.modified)
  # if (length(recipient.index) != 1) {
  #   stop(paste("Tried to find ", flow$recipient[flow.index], 
  #              " but instead found ", paste(rownames(V.modified)[recipient.index], 
  #                                           sep = " ", collapse = " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree", 
  #              sep = ""))
  # }
  donor.index <- which(rownames(V.modified) == flow$donor[flow.index])
  donor.index
  # if (length(donor.index) != 1) {
  #   stop(paste("Tried to find ", flow$donor[flow.index], 
  #              " but instead found ", paste(rownames(V.modified)[donor.index], 
  #                                           sep = " ", collapse = " "), "; make sure the taxon names in the flow dataframe donor match that of your tree", 
  #              sep = ""))
  # }
  V.modified
  
  V.modified[recipient.index, donor.index] <- (1 - gamma) * V.original[recipient.index, donor.index] + (gamma) * (flow$time.from.root.recipient[flow.index]) * sigma.sq
  V.modified
  V.modified[donor.index, recipient.index] <- V.modified[recipient.index, donor.index]
  V.modified
  recipient.index
  V.modified[recipient.index, recipient.index] <- (V.original[recipient.index, recipient.index] - sigma.sq * flow$time.from.root.recipient[flow.index]) + (gamma^2 + (1 - gamma)^2) * (flow$time.from.root.recipient[flow.index]) * sigma.sq + 2 * gamma * (1 - gamma) * V.original[recipient.index, donor.index] + vh

  
  
  if(dim(flow)[1]>1){
    flow.1.index<-flow.index
    recipient.1.index <- which(rownames(V.modified) == flow$recipient[flow.1.index])
    donor.1.index<-donor.index
    for (flow.2.index in sequence(dim(flow)[1])) {
     #flow.2.index<-2
       recipient.2.index<-which(rownames(V.modified)==flow$recipient[flow.2.index])
      recipient.2.index
       donor.2.index<-which(rownames(V.modified)==flow$donor[flow.2.index])
      # case 1: (same donor) flow came in "before" speciation
      if(flow.1.index!=flow.2.index && flow$donor[flow.1.index]==flow$donor[flow.2.index] && V.original[recipient.1.index,recipient.2.index]> max(flow$time.from.root.recipient[flow.1.index],flow$time.from.root.recipient[flow.2.index] ) ){ 
         V.modified[recipient.1.index, recipient.2.index]<-(V.original[recipient.1.index,recipient.2.index] - sigma.sq*flow$time.from.root.recipient[flow.1.index]) + (1-gamma)^2*sigma.sq*flow$time.from.root.recipient[flow.1.index] + gamma^2*sigma.sq*flow$time.from.root.recipient[flow.2.index]  + gamma*(1-gamma)*sigma.sq*V.original[recipient.1.index,donor.1.index] + + gamma*(1-gamma)*sigma.sq*V.original[recipient.2.index,donor.1.index]
        }#end case 1
      # case 2: (same donor) flow came in "after" speciation 
      if(flow.1.index!=flow.2.index && flow$donor[flow.1.index]==flow$donor[flow.2.index] && V.original[recipient.1.index,recipient.2.index] < flow$time.from.root.recipient[flow.2.index] ){ 
        V.modified[recipient.1.index, recipient.2.index]<- (1-gamma)^2*sigma.sq*flow$time.from.root.recipient[flow.1.index] + gamma^2*sigma.sq*flow$time.from.root.recipient[flow.2.index]+gamma*(1-gamma)*sigma.sq*V.original[recipient.1.index,donor.1.index]  + gamma*(1-gamma)*sigma.sq*V.original[recipient.2.index,donor.1.index]                     
      }#end case 2
      #cases 3 different donors
     if(flow.1.index!=flow.2.index && flow$donor[flow.1.index]!=flow$donor[flow.2.index]){
       V.modified[recipient.1.index, recipient.2.index]<- gamma^2*V.original[recipient.1.index,donor.1.index] + (1-gamma)^2*V.original[recipient.2.index,donor.2.index]+gamma*(1-gamma)*V.original[recipient.1.index,donor.2.index] + (1-gamma)*gamma*V.original[recipient.2.index,donor.1.index]
      }#end case 3
     # V.modified[recipient.2.index, recipient.1.index]<- V.modified[recipient.1.index, recipient.2.index]
      }#end for flow.index.2
    }#end if(dim(flow)[1]>1)  

}#end for loop
print(V.modified)
#         
  
  
#


