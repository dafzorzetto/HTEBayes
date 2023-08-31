# CausalHeterogeneity_ConfounderNPBayes

In this repository, we provide the code for the function and simulation study of the paper <a href=https://arxiv.org/abs/2302.11656>_"Confounder-Dependent Bayesian Mixture Model: Characterizing Heterogeneity of Causal Effects in Air Pollution Epidemiology"_ </a> by D. Zorzetto, F.J. Bargagli-Stoffi, A. Canale, and F. Dominici (2023). 

## Example:
Load functions
```R
source("Rcode/src/functions_simulations.R")

```

Generate 
```R
dataset <- simulation_sample_3groups(seed=1,
                                     eta_0=c(2,4,6),
                                     eta_1=c(0,3,6),
                                     sigma_0=rep(0.3,3),
                                     sigma_1=rep(0.3,3),
                                     n=500)
```
