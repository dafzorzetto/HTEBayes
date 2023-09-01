# CausalHeterogeneity_ConfounderNPBayes

In this repository, we provide the code for the function and simulation study of the paper <a href=https://arxiv.org/abs/2302.11656>_"Confounder-Dependent Bayesian Mixture Model: Characterizing Heterogeneity of Causal Effects in Air Pollution Epidemiology"_ </a> by D. Zorzetto, F.J. Bargagli-Stoffi, A. Canale, and F. Dominici (2023). 

### Examples:

Load functions:
```R
source("src/functions_simulations.R")
source("src/model_CDBMM.R")
```
Information for simulation data:
```R
# sample size
n = 500

#set seed
seed = 1
```
Parameters for the model:
```R
# iterations
R = 3000
R_burnin = 2000

# number of maximum groups for the two treatment levels
L_0 = 20
L_1 = 20
```

**Example 1** (*3 groups*)

Generation dataset and estimation CDBMM:
```R
dataset <- simulation_sample_3groups(seed=seed,
                                     eta_0=c(2,4,6),
                                     eta_1=c(0,3,6),
                                     sigma_0=rep(0.3,3),
                                     sigma_1=rep(0.3,3),
                                     n=n)

cdbmm_results <- CDBMM_Gibbs(c=1,
                             data_sample=list(dataset),
                             n=n)
```

Visualize ITEs and GATEs for discovered groups:
```R
hist(cdbmm_results$tau,
     nclass = 50,
     main = "ITEs")

group_allocation<-paste0(cdbmm_results$partition[,1],"-",cdbmm_results$partition[,1])

sapply(unique(group_allocation), function(g) 
  cdbmm_results$atoms$p_1[as.integer(substr(g, 3, 3))]-
    cdbmm_results$atoms$p_0[as.integer(substr(g, 1, 1))])
```

**Example 2** (*5 groups*)

Generation dataset and estimation CDBMM:
```R
dataset <- simulation_sample_5cov(seed=seed,
                                  eta_0=c(2,2,3,4.5,6.5),
                                  eta_1=c(0,1,2.5,5,7.5),
                                  sigma_0=rep(0.2,5),
                                  sigma_1=rep(0.2,5),
                                  n=n)

cdbmm_results <- CDBMM_Gibbs(c=1,
                             data_sample=list(dataset),
                             n=n)
```
Visualize ITEs and GATEs for discovered groups:
```R
hist(cdbmm_results$tau,
     nclass = 50,
     main = "ITEs")

group_allocation<-paste0(cdbmm_results$partition[,1],"-",cdbmm_results$partition[,1])

sapply(unique(group_allocation), function(g) 
  cdbmm_results$atoms$p_1[as.integer(substr(g, 3, 3))]-
    cdbmm_results$atoms$p_0[as.integer(substr(g, 1, 1))])
```
