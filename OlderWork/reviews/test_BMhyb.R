rm(list=ls())
library(ape)
library(BMhyb) # Version 1.5.1, published 2017-09-11 on the CRAN
# (see session infos at the botom of the script)

## Parameters
sigma2 <- 1
x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
actual.params <- c("sigma.sq", "mu", "bt", "vh", "SE")

gamma <- 0.5


###############################################################################
## Computations of the VCV matrix for the network presented in the article
###############################################################################
# I illustrate here how the `time.from.root.donor` parameter might be not
# correctly taken into account.
# I compare the results of the function `GetVModified` to the results of the
# actual formulas of the paper applied to this simple network.

# (Actually, in a quick browse through the function code, I could find no reference
# of the `time.from.root.donor` parameter, that hence does not appear to be taken
# into account into the computations.)

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

## Function to compute the VCV matrix using the formulas of the article,
## for the network presented in it.
vcv_formula_paper <- function(gamma, t1, t2, t3, sigma2) {
  vcv_formula <- matrix(nrow = 3, ncol = 3)
  colnames(vcv_formula) <- c("R", "Y", "X")
  rownames(vcv_formula) <- c("R", "Y", "X")
  vcv_formula[1, 1] <- (gamma^2 + (1-gamma)^2) * (t1 + t2) + t3
  vcv_formula[1, 2] <- vcv_formula[2, 1] <- (1-gamma) * (t1 + t2)
  vcv_formula[1, 3] <- vcv_formula[3, 1] <- gamma * t1
  vcv_formula[2, 2] <- t1 + t2 + t3
  vcv_formula[2, 3] <- vcv_formula[3, 2] <- 0
  vcv_formula[3, 3] <- t1 + t2 + t3
  vcv_formula <- sigma2 * vcv_formula
  return(vcv_formula)
}

###############################################################################
## Network 1

t1 <- 0.3; t2 <- 0.4; t3 <- 0.3; # unit height
network <- create_paper_network(gamma, t1, t2, t3)
## Plot
PlotNetwork(network$phy, network$flow)
axis(1, at = c(0, t1, t1+t2, t1+t2+t3), labels = c("0", "t1", "t1+t2", "t1+t2+t3"))

## VCV
vcv_BMhyb <- GetVModified(x, network$phy, network$flow, actual.params)
vcv_BMhyb
#      R   Y    X
# R 0.65 0.7 0.35
# Y 0.70 1.0 0.00
# X 0.35 0.0 1.00
#V.modified[recipient.index, donor.index] <- (1 - gamma) * V.original[recipient.index, donor.index] + (gamma) * (flow$time.from.root.recipient[flow.index]) * sigma.sq

## VCV formula
vcv_formula <- vcv_formula_paper(gamma, t1, t2, t3, sigma2)
#      R    Y    X
# R 0.65 0.35 0.15
# Y 0.35 1.00 0.00
# X 0.15 0.00 1.00

## Discrepency
vcv_BMhyb - vcv_formula

# Cov[X, R] = sigma^2 * gamma * t1 = 1*0.5*0.3 = 0.15 in this case
# `GetVModified` gives 0.35 for it.
# Cov[Y, R] = sigma^2 * (1-gamma) * (t1 + t2) = 0.35
# `GetVModified` gives 0.7 for it.

###############################################################################
## Network 2: same network, but change donor time

t1 <- 0.7; t2 <- 0; t3 <- 0.3; # same underlying tree, just donor time is different
network2 <- create_paper_network(gamma, t1, t2, t3)
## Plot
PlotNetwork(network2$phy, network2$flow)
axis(1, at = c(0, t1, t1+t2, t1+t2+t3), labels = c("0", "t1", "t1+t2", "t1+t2+t3"))

## VCV function
vcv_BMhyb2 <- GetVModified(x, network2$phy, network2$flow, actual.params) 
# unchanged compared to previous one.
#      R   Y    X
# R 0.65 0.7 0.35
# Y 0.70 1.0 0.00
# X 0.35 0.0 1.00
#V.modified[recipient.index, donor.index] <- (1 - gamma) * V.original[recipient.index, donor.index] + (gamma) * (flow$time.from.root.recipient[flow.index]) * sigma.sq


## VCV formula
vcv_formula2 <- vcv_formula_paper(gamma, t1, t2, t3, sigma2)
#      R    Y    X
# R 0.65 0.35 0.35
# Y 0.35 1.00 0.00
# X 0.35 0.00 1.00

## Discrepency
vcv_BMhyb2 - vcv_formula2

# Cov[X, R] = sigma^2 * gamma * t1 = 0.5*0.7 = 0.35 in this case
# `GetVModified` gives 0.35 for it.
# Cov[Y, R] = sigma^2 * (1-gamma) * (t1 + t2) = 0.35
# `GetVModified` gives 0.7 for it.

###############################################################################
## Network with several tips descendants from on hybrid
###############################################################################
# This illustrate an other potential problem, when several tips are descending
# from a single hybridization event.

# (Going through the code, I noticed that "V[recipient1, recipient2]" does not
# seem to be modified, which is coherent with what I find here.)

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

# Cov[Y,R] = sigma^2 * [(gamma^2 + (1-gamma)^2)*t1 + t2] = 0.55
# `GetVModified` gives 0.7 (same as a tree).
# The ancestral hybridization event is not correctly taken into account.

###############################################################################
###############################################################################

# As a comparison, here is the Julia PhyloNetworks code producing the results
# for these three examples.
#
# Ex 1
# julia> using PhyloNetworks
# julia> net = readTopology("((Y:0.3,(R:0.3)#H1:0::0.5):0.7,(#H1:0.4::0.5,X:0.7):0.3);")
# julia> sharedPathMatrix(truenet)[:Tips]
# 3×3 Array{Float64,2}:
#   1.0   0.35  0.0 
#   0.35  0.65  0.15
#   0.0   0.15  1.0 
#
# Ex 2
# julia> net = readTopology("((Y:0.3,(R:0.3)#H1:0::0.5):0.7,(#H1:0::0.5,X:0.3):0.7);")
# julia> sharedPathMatrix(truenet)[:Tips]
# 3×3 Array{Float64,2}:
#   1.0   0.35  0.0 
#   0.35  0.65  0.35
#   0.0   0.35  1.0 
#
# Ex 3
# julia> net = readTopology("(((Y:0.3,R:0.3):0.4)#H1:0.3::0.5,(#H1:0::0.5,X:0.7):0.3);")
# julia> sharedPathMatrix(truenet)[:Tips]
# 3×3 Array{Float64,2}:
#   0.85  0.55  0.15
#   0.55  0.85  0.15
#   0.15  0.15  1.0 

###############################################################################
###############################################################################
# > sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.18.so
# 
# locale:
# [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
# [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] BMhyb_1.5.1 ape_4.1    
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_0.12.12            subplex_1.4-1           msm_1.6.4               mvtnorm_1.0-6          
# [5] lattice_0.20-35         tidyr_0.7.1             corpcor_1.6.9           prettyunits_1.0.2      
# [9] assertthat_0.2.0        digest_0.6.12           foreach_1.4.3           R6_2.2.2               
# [13] plyr_1.8.4              phytools_0.6-20         coda_0.19-1             httr_1.3.1             
# [17] ggplot2_2.2.1           progress_1.1.2          rlang_0.1.2.9000        uuid_0.1-2             
# [21] lazyeval_0.2.0          curl_2.8.1              data.table_1.10.4       taxize_0.9.0           
# [25] phangorn_2.2.0          Matrix_1.2-11           RNeXML_2.0.7            combinat_0.0-8         
# [29] splines_3.4.1           stringr_1.2.0           igraph_1.1.2            munsell_0.4.3          
# [33] compiler_3.4.1          numDeriv_2016.8-1       geiger_2.0.6            pkgconfig_2.0.1        
# [37] mnormt_1.5-5            tibble_1.3.4            gridExtra_2.2.1         TreeSim_2.3            
# [41] expm_0.999-2            quadprog_1.5-5          codetools_0.2-15        XML_3.98-1.9           
# [45] reshape_0.8.7           viridisLite_0.2.0       dplyr_0.7.2             MASS_7.3-47            
# [49] crul_0.3.8              grid_3.4.1              nlme_3.1-131            jsonlite_1.5           
# [53] gtable_0.2.0            magrittr_1.5            scales_0.5.0            stringi_1.1.5          
# [57] reshape2_1.4.2          viridis_0.4.0           bindrcpp_0.2            scatterplot3d_0.3-40   
# [61] phylobase_0.8.4         xml2_1.1.1              fastmatch_1.1-0         deSolve_1.20           
# [65] iterators_1.0.8         tools_3.4.1             rncl_0.8.2              ade4_1.7-8             
# [69] bold_0.5.0              glue_1.1.1              purrr_0.2.3             maps_3.2.0             
# [73] plotrix_3.6-6           parallel_3.4.1          survival_2.41-3         colorspace_1.3-2       
# [77] bindr_0.1               animation_2.5           clusterGeneration_1.3.4