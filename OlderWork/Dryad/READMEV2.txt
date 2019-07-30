This data package contains all relevant R and RData files that generate the Fig. 3, Fig. 4, Fig. 5, Fig. S3, and Table 1 for the submission manuscript: “Trait Evolution on Phylogenetic Networks” by Dwueng-Chwuan Jhwueng and Brian C. O’Meara.

The raw code for the R package BMhyd is included in this data package as well. It is also on CRAN, but the CRAN code will continue to evolve, while this gives a static snapshot of the code used for the analyses. See the following for detailed description of data. 

——————————————
Figure 3
——————————————
R file “test.mean.r” generates the RData file “test.mean.RData” that mainly contains the data array “sim.mean.array” used for testing the parameter of the mean using the rooted phylogenetic network of 19 species of cichlids.  

R file “test.bt.r” generates the RData file “test.bt.RData” that mainly contains the data array “sim.bt.array” used for testing the parameter of the hybrid vigor using the rooted phylogenetic network of 19 species of cichlids.  

R file “test.sigma_sq.” generates the RData file “test.sigma_sq.RData” that mainly contains the data array “sim.sigma.array” used for testing the parameter of the variance using the rooted phylogenetic network of 19 species of cichlids.  
  

R file “test.v.H.r” generates the RData file “test.v.H.RData” that mainly contains the data array “sim.sigma.h.array” used for testing the parameter of the hybrid variation using the rooted phylogenetic network of 19 species of cichlids.  


R file “test.h.r” generates the RData file “test.h.RData” that mainly contains the data array “sim.h.array” used for testing the parameter of the hybrid proportion using the rooted phylogenetic network of 19 species of cichlids.  
  

R file “sim.test.h.moretaxa.r” generates the RData file “sim.test.h.moretaxa.RData” that mainly contains the data array “h.est.array” used for testing the identifiability of the parameter of the hybrid proportion by increasing the number of taxa.  

R file “Figure3CombineBCO.r” summarizes the results and produces the Fig.3 in the manuscript.

——————————————
Figure 4 
——————————————
R file “bt.v.H.r” generates the RData files “bt.vs.v.H.RData” that mainly contains the data array
“NegLogLike.array” for the purpose of plotting the contour plot of parameter beta and v.H near their MLEs using the rooted phylogenetic network of 19 species of cichlids.  


R file “h.vs.bt.r” generates the RData file “h.vs.bt.RData” that mainly contains the data array
“NegLogLike.array” for the purpose of plotting the contour plot of parameter h and beta near their MLEs using the rooted phylogenetic network of 19 species of cichlids.  


R file “h.vs.v.H.r” generates the RData file “h.vs.v.H.RData” that mainly contains the data array
“NegLogLike.array” for the purpose of plotting the contour plot of parameter h and v.H near their MLEs using the rooted phylogenetic network of 19 species of cichlids.  

R file “Figure4CombineBCO.r” summarizes the results and produces the Fig. 4 in the manuscript.



——————————————
Figure 5
——————————————
R file “sim.test.bt_time_variance.r” generates the RData file  “sim.test.bt_time_variance.RData” which mainly contains the data 
bt.output.array for the hybrid vigor parameter estimates under different time of hybridzation. 

R file “sim.test.v.H_time_variance.r” generates the RData file “sim.test.v.H_time_variance.RData” which mainly contains the data 
“v.H.output.array” for the hybrid variation parameter estimates under different time of hybridization. 


R file “Figure_5.r” summarizes results of “sim.test.bt_time_variance.RData” and “sim.test.v.H_time_variance.RData” and produces the Fig. 5 in the manuscript.


——————————————
Figure S3
——————————————
R file “sim.test.bt.r” generates the RData file “sim.test.bt.RData” which mainly contains the data array “bt.output.array” that represents the bias of the hybrid vigor parameters accessed in different arrangements of hybrid times, number of taxa and proportion of hybrid removed from the network.  
  
R file “sim.test.mu.r” generates the RData file “sim.test.mu.RData” which mainly contains the data array “mu.output.array” that represents the bias of the hybrid vigor parameters accessed in different arrangements of hybrid times, number of taxa and proportion of hybrid removed from the network.  
  
R file “sim.test.sigma_sq.r” generates the RData file “sim.test.sigma_sq.RData” which mainly contains the data array “sigma_sq.output.array” that represents the bias of the hybrid vigor parameters accessed in different arrangements of hybrid times, number of taxa and proportion of hybrid removed from the network. 

R file “sim.test.v.H.r” generates the RData file “sim.test.v.H.RData” which mainly contains the data array “v.H.output.array” that represents the bias of the hybrid vigor parameters accessed in different arrangements of hybrid times, number of taxa and proportion of hybrid removed from the network.  

R file “Figure_S3.r” summarizes above results and produces Fig. S3 in the online supplement file.



——————————————
Table 1
——————————————
R file “Cichlid.Data.se.4models.r” generates the RData file “Cichlid.Data.se.4models.RData” which includes the data array “output.array” of size 4 x 12 and data array “summary.weight” of size 4 x 2.  The row of the table of “output.array” and “summary weight” represents the data analysis under specific model, and the column represents the parameter estimates, negative log likelihood, AICc values, delta AICc values, Akaike weight, cumulative Akaike weight, and weight for parameters estimates while the columns of “summary.weight” represents the weighted averaged and the standard error of parameter estimates using the rooted phylogenetic network of 19 cichlids and the trait data of the body lengths.  

The full result can be seen in Table 1 in the manuscript. 



——————————————
Raw_code_BMhyd
——————————————
R file “bmhyd.r” is the main code for the R-package BMhyd. For description, please see http://cran.at.r-project.org/web/packages/BMhyd/index.html.


contact: Dwueng-Chwuan Jhwueng <djhwueng@umail.iu.edu>




