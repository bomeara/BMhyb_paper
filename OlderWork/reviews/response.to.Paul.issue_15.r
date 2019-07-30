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


#issue 15 Variance between hybrid descendants

#Hi again, @bomeara and @djhwueng
#This might be related to #13 and #14.
#I tried a network with several hybridization events:

gamma1 <- 0.5; gamma2 <- 0.5;
## Underlying tree
t1 <- 0.2; t2 <- 0.2; t3 <- 0.2; t4 <- 0.2; t5 <- 0.2;
phy <- read.tree(text = paste0("(((R:",t4+t5,",Y:",t4+t5,"):",t3,",X:",t3+t4+t5,"):",t1+t2,",Z:",t1+t2+t3+t4+t5,");"))
plot(phy)
## Network
don_recp <- rbind(expand.grid(c("Z"), c("Y", "R", "X")),
                  expand.grid(c("X"), c("R")))
network <- list(phy = phy,
                flow = data.frame(donor = don_recp[,1],
                                  recipient = don_recp[,2],
                                  gamma = c(rep(gamma1, 3), gamma2),
                                  time.from.root.donor = c(rep(t1, 3), t1+t2+t3+t4),
                                  time.from.root.recipient = c(rep(t1, 3), t1+t2+t3+t4)))
network$flow$donor <- as.character(network$flow$donor)
network$flow$recipient <- as.character(network$flow$recipient)
## Plot
PlotNetwork(network$phy, network$flow)
axis(1, at = c(0, t1, t1+t2, t1+t2+t3, t1+t2+t3+t4, t1+t2+t3+t4+t5),
     labels = c("0", "t1", "t1+t2", "t1+t2+t3", "t1+t2+t3+t4", "t1+t2+t3+t4+t5"))
#This gives the folowing variance matrix:

sigma2 = 1
x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
actual.params <- c("sigma.sq", "mu", "bt", "vh", "SE")

GetVModified(x, network$phy, network$flow, actual.params)
#   R   Y   X   Z
#R 0.8 0.6 0.6 0.1
#Y 0.6 0.9 0.4 0.1
#X 0.6 0.4 0.9 0.1
#Z 0.1 0.1 0.1 1.0
#I think that the variance of R is not coherent with the model of trait evolution. If my computations are correct, we should have:

########
#Var[R] = (gamma2^2 + (1-gamma2)^2)*((gamma1^2 + (1-gamma1)^2)*t1+t2+t3+t4) + 2*gamma2*(1-gamma2)*((gamma1^2 + (1-gamma1)^2)*t1+t2) + t5     = 0.7 \neq 0.8
########
#(Note that the covariances between R and Y and X might also have problems, see #14).

#Browsing through the code, this might be linked with the fact that a new hybridization "erases" an older one in your algorithm. Indeed, all the computations are made using V.original, that do not take ancestral hybrids into account. Here, if there were only one hybridization (the second one), then we would have:

#Var[R] = (gamma2^2 + (1-gamma2)^2)*(t1+t2+t3+t4) + 2*gamma2*(1-gamma2)*(t1+t2) + t5 = 0.8 which is the result given by GetVModified.

#I think this is a seperate problem from the two other ones, hence the new issue. Again, I'm sorry if I mis-used your functions or made mistakes, please correct me if I did.

#Thanks !
