#########################################################################
#        ---  SIMULATION STUDY    ----
#########################################################################

# load functions
source("src/functions_simulations.R")
source("src/plots_simulations.R")

#########################################################################
# -----    generation data   -----
#########################################################################

n=500          # units for each sample
samples=100    # repetition of each setting (number of samples)

scenario_1=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(2,4,6),
                            eta_1=c(0,3,6),
                            sigma_0=rep(0.3,3),
                            sigma_1=rep(0.3,3),
                            n=n))
scenario_2=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(0,2.2,4.4),
                            eta_1=c(0,0,0),
                            sigma_0=rep(0.2,3),
                            sigma_1=rep(0.2,3),
                            n=n))
scenario_3=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(1,2,3),
                            eta_1=c(0,1.5,3),
                            sigma_0=c(0.2,0.25,0.25),
                            sigma_1=c(0.25,0.3,0.2),
                            n=n))
scenario_4=lapply(1:samples, function(c) 
  simulation_sample_4groups(seed=c,
                            eta_0=c(1,2,3,3),
                            eta_1=c(0,1.5,3,4.5),
                            sigma_0=rep(0.2,4),
                            sigma_1=rep(0.2,4),
                            n=n))

# 5 covariates - 5 groups
scenario_5=lapply(1:samples, function(c) 
  simulation_sample_5cov(seed=c,
                         eta_0=c(2,2,3,4.5,6.5),
                         eta_1=c(0,1,2.5,5,7.5),
                         sigma_0=rep(0.2,5),
                         sigma_1=rep(0.2,5),
                         n=n))

# 3 covariates - 3 groups - closer groups
scenario_6=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(1.5,2,2.5),
                            eta_1=c(1,1.75,2.5),
                            sigma_0=rep(0.3,3),
                            sigma_1=rep(0.3,3),
                            n=n))

# no heterogeneity - no groups
scenario_7=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=rep(2,3),
                            eta_1=rep(3,3),
                            sigma_0=rep(0.5,3),
                            sigma_1=rep(0.5,3),
                            n=n))

save.image("data_simulations.RData")

#########################################################################
# -----    visualization data   -----
#########################################################################

# histogram for the outcomes
hist_sim(data=scenario_1[[1]],
         scenario_n=1)
hist_sim(data=scenario_2[[1]],
         scenario_n=2)
hist_sim(data=scenario_3[[1]],
         scenario_n=3)
hist_sim(data=scenario_4[[1]],
         scenario_n=4)
hist_sim(data=scenario_5[[1]],
         scenario_n=5)
hist_sim(data=scenario_6[[1]],
         scenario_n=6)
hist_sim(data=scenario_7[[1]],
         scenario_n=7)

# histogram for the individual treatment effects
hist_ITE_sim(data=scenario_1[[1]],
             scenario_n=1)
hist_ITE_sim(data=scenario_2[[1]],
         scenario_n=2)
hist_ITE_sim(data=scenario_3[[1]],
         scenario_n=3)
hist_ITE_sim(data=scenario_4[[1]],
         scenario_n=4)
hist_ITE_sim(data=scenario_5[[1]],
         scenario_n=5)
hist_ITE_sim(data=scenario_6[[1]],
         scenario_n=6)
hist_ITE_sim(data=scenario_7[[1]],
         scenario_n=7)
