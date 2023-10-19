## Code files:

### General functions (in **`src`** folder):
 - **`functions_simulations.R`**:
   
    functions to simulate different settings
 - **`model_CDBMM.R`**:

   Gibbs sampler to estimate our proposed model: Confounder-Dependent Bayesian Mixture Model
 - **`model_BART_BCF_CART.R`**:

   functions to estimate BART and BCF model

   estimate the groups with CART
 - **`plots_simulations.R`**:

   functions to visualize the simulated data

 - **`visualization_GATEs.R`**:

   function to visualize the estimated GATEs (boxplots in Supplementary Material)

### Reproduce the results in the Simulation Study Section (in **`simulation study`** folder):
 - **`1_generating_data_simulations.R`**:

   generate the seven different settings
 - **`2_estimation_models.R`**:

   estimate CDBMM, BART, and BCF for the seven different settings
 - **`3_bias_MSE.R`**:

   estimate bias and MSE for all the models and settings
 - **`4_visualizing_results.R`**:

   functions to visualize the results: bias and MSE

### Supplementary materials (in **`supplementary materials`** folder):
 - **`run_time.R`**:
   
    functions to compare the run times for the models (CDBMM, BART, BCF) for different sample sizes

- **`sensitivity_analysis.R`**:
   
    functions for sensitivity analysis of hyperparameters for regression parameters in probit function in the weights of the mixtures
