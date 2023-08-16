#########################################################################
#       ---         COMPARISON        ----
#       --- RUN TIME and SAMPLE SIZE  ----
#########################################################################

# load functions
source("0_functions_simulations.R")
source("0_model_CDBMM.R")
source("0_model_BART_BCF_CART.R")

# Library for performing benchmarking
library(rbenchmark)
library(xtable)

#########################################################################
#  --- generating data changing sample size ---
#########################################################################

# sample sizes
n_ss=c(100,250,500,1000)  #,2000)

scenario_1=lapply(n_ss, function(c) 
  simulation_sample_3groups(seed=1,
                            eta_0=c(2,4,6),
                            eta_1=c(0,3,6),
                            sigma_0=rep(0.3,3),
                            sigma_1=rep(0.3,3),
                            n=c))
scenario_2=lapply(n_ss, function(c) 
  simulation_sample_3groups(seed=1,
                            eta_0=c(0,2.2,4.4),
                            eta_1=c(0,0,0),
                            sigma_0=rep(0.2,3),
                            sigma_1=rep(0.2,3),
                            n=c))
scenario_3=lapply(n_ss, function(c) 
  simulation_sample_3groups(seed=1,
                            eta_0=c(1,2,3),
                            eta_1=c(0,1.5,3),
                            sigma_0=c(0.2,0.25,0.25),
                            sigma_1=c(0.25,0.3,0.2),
                            n=c))
scenario_4=lapply(n_ss, function(c) 
  simulation_sample_4groups(seed=1,
                            eta_0=c(1,2,3,3),
                            eta_1=c(0,1.5,3,4.5),
                            sigma_0=rep(0.2,4),
                            sigma_1=rep(0.2,4),
                            n=c))
scenario_5=lapply(n_ss, function(c) 
  simulation_sample_5cov(seed=1,
                         eta_0=c(2,2,3,4.5,6.5),
                         eta_1=c(0,1,2.5,5,7.5),
                         sigma_0=rep(0.2,5),
                         sigma_1=rep(0.2,5),
                         n=c))
scenario_6=lapply(n_ss, function(c) 
  simulation_sample_3groups(seed=1,
                            eta_0=c(1.5,2,2.5),
                            eta_1=c(1,1.75,2.5),
                            sigma_0=rep(0.3,3),
                            sigma_1=rep(0.3,3),
                            n=c))
scenario_7=lapply(n_ss, function(c) 
  simulation_sample_3groups(seed=1,
                            eta_0=rep(2,3),
                            eta_1=rep(3,3),
                            sigma_0=rep(0.5,3),
                            sigma_1=rep(0.5,3),
                            n=c))


#########################################################################
#  --- comparison ---
#########################################################################

# iterations
R=2000
R_burnin=1000

#useful values
L_0=10
L_1=10

time_comparison_f<-function(sample_size, scenario){

  # time comparison
  time_comp<-benchmark(
    CDBMM = CDBMM_Gibbs(c=sample_size, data_sample=scenario, n=length(scenario[[sample_size]]$data$T)),
    BART = bart_sample(c=sample_size, data_sample=scenario),
    BCF = BCF_sample(c=sample_size, data_sample=scenario),
    columns = c("test", "elapsed", "relative"),
    replications = 1
  )
  
  return(time_comp)
}

time_s1<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_1))
time_s2<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_2))
time_s3<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_3))
time_s4<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_4))
time_s5<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_5))
time_s6<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_6))
time_s7<-lapply(1:(length(n_ss)), function(c) 
  time_comparison_f(sample_size=c,
                    scenario=scenario_7))


xtable(rbind(time_s1,time_s2,time_s3,time_s4,time_s5,time_s6,time_s7))

save.image("run_time.RData")
