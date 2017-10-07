rm(list=ls())

library(ape)
library(BMhyb)

#modified version for GetVModified function after Paul's finding Oct. 06 2017
GetVModified<-function (x, phy, flow, actual.params, measurement.error = NULL)
{
    bt <- 1
    vh <- 0
    sigma.sq <- x[1]
    mu <- x[2]
    SE <- 0
    if (is.null(measurement.error)) {
        SE <- x[length(x)]
    }
    bt.location <- which(names(actual.params) == "bt")
    if (length(bt.location) == 1) {
        bt <- x[bt.location]
    }
    vh.location <- which(names(actual.params) == "vh")
    if (length(vh.location) == 1) {
        vh <- x[vh.location]
    }
    times.original <- vcv.phylo(phy, model = "Brownian")
    V.original <- sigma.sq * times.original
    V.modified <- V.original
    for (flow.index in sequence(dim(flow)[1])) {
        recipient.index <- which(rownames(V.modified) == flow$recipient[flow.index])
        gamma <- flow$gamma[flow.index]
        if (length(recipient.index) != 1) {
            stop(paste("Tried to find ", flow$recipient[flow.index],
                " but instead found ", paste(rownames(V.modified)[recipient.index],
                  sep = " ", collapse = " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree",
                sep = ""))
        }
        donor.index <- which(rownames(V.modified) == flow$donor[flow.index])
        if (length(donor.index) != 1) {
            stop(paste("Tried to find ", flow$donor[flow.index],
                " but instead found ", paste(rownames(V.modified)[donor.index],
                  sep = " ", collapse = " "), "; make sure the taxon names in the flow dataframe donor match that of your tree",
                sep = ""))
        }
        V.modified[recipient.index, donor.index] <- (1 - gamma) *
            V.original[recipient.index, donor.index] + (gamma) *
            (flow$time.from.root.recipient[flow.index]) * sigma.sq
        V.modified[donor.index, recipient.index] <- V.modified[recipient.index,
            donor.index]
        V.modified[recipient.index, recipient.index] <- (V.original[recipient.index,
            recipient.index] - sigma.sq * flow$time.from.root.recipient[flow.index]) +
            (gamma^2 + (1 - gamma)^2) * (flow$time.from.root.recipient[flow.index]) *
                sigma.sq + 2 * gamma * (1 - gamma) * V.original[recipient.index,
            donor.index] + vh


     if(dim(flow)[1]>1){
     flow.1.index<-flow.index
     recipient.1.index <- which(rownames(V.modified) == flow$recipient[flow.1.index])
     donor.1.index<-donor.index
     for (flow.2.index in sequence(dim(flow)[1])) {
       recipient.2.index<-which(rownames(V.modified)==flow$recipient[flow.2.index])
       donor.2.index<-which(rownames(V.modified)==flow$donor[flow.2.index])
       # case 1: (same donor) flow came in "before" speciation
       if(flow.1.index!=flow.2.index && flow$donor[flow.1.index]==flow$donor[flow.2.index] && V.original[recipient.1.index,recipient.2.index]> max(flow$time.from.root.recipient[flow.1.index],flow$time.from.root.recipient[flow.2.index] ) ){
          V.modified[recipient.1.index, recipient.2.index]<-(V.original[recipient.1.index,recipient.2.index] - sigma.sq*flow$time.from.root.recipient[flow.1.index]) + (1-gamma)^2*sigma.sq*flow$time.from.root.recipient[flow.1.index] + gamma^2*sigma.sq*flow$time.from.root.recipient[flow.2.index] + gamma*(1-gamma)*sigma.sq*(V.original[recipient.1.index,donor.index]+V.original[recipient.2.index,donor.index])
          }#end case 1
       # case 2: (same donor) flow came in "after" speciation
       if(flow.1.index!=flow.2.index && flow$donor[flow.1.index]==flow$donor[flow.2.index] && V.original[recipient.1.index,recipient.2.index] < min(flow$time.from.root.recipient[flow.1.index], flow$time.from.root.recipient[flow.2.index]) ){
         V.modified[recipient.1.index, recipient.2.index]<- (1-gamma)^2*sigma.sq*flow$time.from.root.recipient[flow.1.index] + gamma^2*sigma.sq*flow$time.from.root.recipient[flow.2.index]+gamma*(1-gamma)*sigma.sq*(V.original[recipient.1.index,donor.index]+V.original[recipient.2.index,donor.index])
         }#end case 2
       }#end for flow.index.2
     }#end if(dim(flow)[1]>1)

    }
    if (is.null(measurement.error)) {
        diag(V.modified) <- diag(V.modified) + SE^2
    }
    else {
        diag(V.modified) <- diag(V.modified) + measurement.error^2
    }
    return(V.modified)
}


## Parameters
sigma2 <- 1
x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
actual.params <- c("sigma.sq", "mu", "bt", "vh", "SE")

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
##  Network with several tips descendants from on hybrid
###############################################################################
# This illustrate an other potential problem, when several tips are descending
# from a single hybridization event.

# (Going through the code, I noticed that "V[recipient1, recipient2]" does not
# seem to be modified, which is coherent with what I find here.)

#######################################################
##CASE 1: (same donor) flow came in "before" speciation
#######################################################
## Underlying tree
t1 <- 0.3; t2 <- 0.4; t3 <- 0.3;
phy <- read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
plot(phy)
## Network
don_recp <- expand.grid(c("X"), c("Y", "R"))
network <- list(phy = phy,
                flow = data.frame(donor = don_recp[,1],
                                  recipient = don_recp[,2],
                                  gamma = rep(gamma, 2),
                                  time.from.root.donor = rep(t1, 2),
                                  time.from.root.recipient = rep(t1, 2)))
network$flow$donor <- as.character(network$flow$donor)
network$flow$recipient <- as.character(network$flow$recipient)
## Plot
PlotNetwork(network$phy, network$flow)
axis(1, at = c(0, t1, t1+t2, t1+t2+t3), labels = c("0", "t1", "t1+t2", "t1+t2+t3"))

## VCV network
GetVModified(x, network$phy, network$flow, actual.params)
#   R    Y    X
#R 0.85 0.55 0.15
#Y 0.55 0.85 0.15
#X 0.15 0.15 1.00
# Cov[Y,R] = sigma^2 * [(gamma^2 + (1-gamma)^2)*t1 + t2] = 0.55
# `GetVModified` gives 0.55 same as the forumala.
# The ancestral hybridization event is taken into account.



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
#    X    R1    R2
#X  1.00 0.500 0.350
#R1 0.50 0.800 0.425
#R2 0.35 0.425 0.650
