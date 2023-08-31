#########################################################################
#       ---     ESTIMATION MODELS     ----
#       --- with CDBMM, BART, and BCF+CART  ----
#########################################################################

# libraries
library(parallel)

# load simulated data
load("data_simulations.RData")

# load functions
source("src/model_CDBMM.R")
source("src/model_BART_BCF_CART.R")

#########################################################################
#    ---  Estimation for the simulated settings with CDBMM    ----
#########################################################################

# iterations
R=3000
R_burnin=2000

# number of maximum groups for the two treatment levels
L_0=12
L_1=12

# we are using the library "parallel"
# mc.cores = number of used cores 
CDBMM_scenario_1=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_1, n=n, mc.cores = 6)
CDBMM_scenario_2=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_2, n=n, mc.cores = 6)
CDBMM_scenario_3=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_3, n=n, mc.cores = 6)
CDBMM_scenario_4=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_4, n=n, mc.cores = 6)
CDBMM_scenario_5=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_5, n=n, mc.cores = 6)
CDBMM_scenario_6=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_6, n=n, mc.cores = 6)
CDBMM_scenario_7=mclapply(1:samples, CDBMM_Gibbs, data_sample=scenario_7, n=n, mc.cores = 6)

#########################################################################
#    ---  Estimation with BART    ----
#########################################################################

# estimation for the simulated samples
BART_scenario_1=mclapply(1:samples, bart_sample, data_sample=scenario_1, mc.cores = 6)
BART_scenario_2=mclapply(1:samples, bart_sample, data_sample=scenario_2, mc.cores = 6)
BART_scenario_3=mclapply(1:samples, bart_sample, data_sample=scenario_3, mc.cores = 6)
BART_scenario_4=mclapply(1:samples, bart_sample, data_sample=scenario_4, mc.cores = 6)
BART_scenario_5=mclapply(1:samples, bart_sample, data_sample=scenario_5, mc.cores = 6)
BART_scenario_6=mclapply(1:samples, bart_sample, data_sample=scenario_6, mc.cores = 6)
BART_scenario_7=mclapply(1:samples, bart_sample, data_sample=scenario_7, mc.cores = 6)


#########################################################################
#    ---  Estimation with BCF    ----
#########################################################################

# estimation for the simulated samples
# we can not parallel BCF 
# due to BCF saves information in file that is called with a name that depend on run time
# running in the same time more samples the file overwrite and we loose information
BCF_scenario_1=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_1))
BCF_scenario_2=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_2))
BCF_scenario_3=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_3))
BCF_scenario_4=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_4))
BCF_scenario_5=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_5))
BCF_scenario_6=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_6))
BCF_scenario_7=lapply(1:samples, function(s) BCF_sample(s, data_sample=scenario_7))


#########################################################################
#    ---  Clustering with CART    ----
#########################################################################

# estimation for the simulated samples
CART_scenario_1=mclapply(1:samples, CART, 
                         data_sample=scenario_1, 
                         estimated_Y=BCF_scenario_1, 
                         mc.cores = 6)
CART_scenario_2=mclapply(1:samples, CART, 
                         data_sample=scenario_2, 
                         estimated_Y=BCF_scenario_2, 
                         mc.cores = 6)
CART_scenario_3=mclapply(1:samples, CART, 
                         data_sample=scenario_3, 
                         estimated_Y=BCF_scenario_3, 
                         mc.cores = 6)
CART_scenario_4=mclapply(1:samples, CART, 
                         data_sample=scenario_4, 
                         estimated_Y=BCF_scenario_4, 
                         mc.cores = 6)
CART_scenario_5=mclapply(1:samples, CART, 
                         data_sample=scenario_5, 
                         estimated_Y=BCF_scenario_5, 
                         mc.cores = 6)
CART_scenario_6=mclapply(1:samples, CART, 
                         data_sample=scenario_6, 
                         estimated_Y=BCF_scenario_6, 
                         mc.cores = 6)
CART_scenario_7=mclapply(1:samples, CART, 
                         data_sample=scenario_7, 
                         estimated_Y=BCF_scenario_7, 
                         mc.cores = 6)

save.image("results_models.RData")
