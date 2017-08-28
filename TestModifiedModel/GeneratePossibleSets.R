non.hybrid.vector <- c(30, 100)
hybrid.vector <- c(1, 5, 10)
flow.vector <-c(.5)
origins.vector <- c("clade", "individual")
sigmasq.vector <- c(0.01) 
tree.height.vector <- c(50)
mu.vector <- c(1)
bt.vector <- c(0.5, 1, 2)
vh.vector <- c(0, .1*max(sigmasq.vector)*max(tree.height.vector), max(sigmasq.vector)*max(tree.height.vector)) #zero, as much variance as you get from 10% of evo history, or as much variance as you get from the whole tree
SE.vector <- c(0, 0.5*sqrt(max(sigmasq.vector)*max(tree.height.vector)))
nreps <- 25
possible.sets <- expand.grid(ntax.nonhybrid=non.hybrid.vector, ntax.hybrid=hybrid.vector, flow.proportion=flow.vector, origin.type=origins.vector, sigma.sq=sigmasq.vector, mu=mu.vector, bt=bt.vector, vh=vh.vector, SE=SE.vector, tree.height=tree.height.vector)

dim(possible.sets)
