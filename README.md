# CausalHeterogeneity_ConfounderNPBayes

In this repository, we provide the code for the function and simulation study of the paper <a href=https://arxiv.org/abs/2302.11656>_"Confounder-Dependent Bayesian Mixture Model: Characterizing Heterogeneity of Causal Effects in Air Pollution Epidemiology"_ </a> by D. Zorzetto, F.J. Bargagli-Stoffi, A. Canale, and F. Dominici (2023). 

## Code files:

### general functions
 - "0_functions_simulations.R":
    functions to simulate different settings
 - "0_model_CDBMM.R":
    Gibbs sampler to estimate our proposed model: Confounder-Dependent Bayesian Mixture Model
 - "0_model_BART_BCF_CART.R":
    functions to estimate BART and BCF model
    estimate the groups with CART
 - "0_plots_simulations.R":
    functions to visualize the simulated data

### reproduce the results in Simulation Study Section
 - "1_generating_data_simulations.R":
   generate the seven different settings
 - "2_estimation_models.R":
   estimate CDBMM, BART, and BCF for the seven different settings
 - "3_bias_MSE.R":
   estimate bias and MSE for all the models and settings
 - "4_visualizing_results.R":
   functions to visualize the results: bias and MSE
