#########################################################################
#        ---  SIMULATION STUDY    ----
#########################################################################

#library
library(mvtnorm)
library(ggplot2)
library(ggExtra)

#########################################################################

# --- generation functions ---

n=800          # units for each sample
samples=100    # repetition of each setting (number of samples)

# function with 2 COVARIATES and 3 GROUPS
simulation_sample_3groups<-function(seed,eta_0,eta_1,sigma_0,sigma_1){
  
  # set seed for riproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectivelly
  X=cbind(rbinom(n,1,0.4),
          rbinom(n,1,0.6))
  
  # treatment
  # logit of a covariates function
  reg_T=0+0.4*(X[,1])+0.6*(X[,2])
  logit_T=exp(reg_T)/(1+exp(reg_T))
  T=rbinom(n,1,logit_T)
  
  # cluster allocation
  S_dummy=cbind(apply(X-1,1,prod),
                (X[,1]==1),
                (X[,1]==0 & X[,2]==1))
  S_cluster=S_dummy%*%(1:3)
  
  # outcomes 
  # normal-bivariate 
  # the covariance matrix is diagonal 
  Y=sapply(1:n, function(i) rmvnorm(1,c(eta_0[S_cluster[i]],eta_1[S_cluster[i]]),
                                    c(sigma_0[S_cluster[i]],sigma_1[S_cluster[i]])*diag(2)))
  
  # save all the information as a list
  return(list(data=list(X=X,T=T,Y=Y),                            # simulated data
              parameters=list(eta_0=eta_0, eta_1=eta_1,          # true parameters
                              sigma_0=sigma_0,sigma_1=sigma_1),  
              S_cluster=S_cluster))                              # cluster allocation
}

# function with 2 COVARIATES and 4 GROUPS
simulation_sample_4groups<-function(seed,eta_0,eta_1,sigma_0,sigma_1){
  
  # set seed for riproducibility
  set.seed(seed)
  
  # covariates
  # 2 bernulli with prob=0.4 and prob=0.6 respectivelly
  X=cbind(rbinom(n,1,0.4),
          rbinom(n,1,0.6))
  
  # treatment
  # logit of a covariates function
  reg_T=0+0.4*(X[,1])+0.6*(X[,2])
  logit_T=exp(reg_T)/(1+exp(reg_T))
  T=rbinom(n,1,logit_T)
  
  # cluster allocation
  S_dummy=cbind((X[,1]==0 & X[,2]==1),
                (X[,1]==0 & X[,2]==0),
                (X[,1]==1 & X[,2]==1),
                (X[,1]==1 & X[,2]==0))
  S_cluster=S_dummy%*%(1:4)
  
  # outcomes 
  # normal-bivariate 
  # the covariance matrix is diagonal 
  Y=sapply(1:n, function(i) rmvnorm(1,c(eta_0[S_cluster[i]],eta_1[S_cluster[i]]),
                                    c(sigma_0[S_cluster[i]],sigma_1[S_cluster[i]])*diag(2)))
  
  # save all the information as a list
  return(list(data=list(X=X,T=T,Y=Y),                            # simulated data
              parameters=list(eta_0=eta_0, eta_1=eta_1,          # true parameters
                              sigma_0=sigma_0,sigma_1=sigma_1),  
              S_cluster=S_cluster))                              # cluster allocation
}

# function with 5 COVARIATES and 5 GROUPS
simulation_sample_5cov<-function(seed,eta_0,eta_1,sigma_0,sigma_1){
  
  # set seed for riproducibility
  set.seed(seed)
  
  # covariates
  # 5 bernulli
  X=cbind(rbinom(n,1,0.4),
          rbinom(n,1,0.6),
          rbinom(n,1,0.3),
          rbinom(n,1,0.5),
          rbinom(n,1,0.2))
  
  # treatment
  # logit of a covariates function
  reg_T=0+0.4*(X[,1])+0.6*(X[,2])-0.3*(X[,3])+0.2*(X[,4])*(X[,5])
  logit_T=exp(reg_T)/(1+exp(reg_T))
  T=rbinom(n,1,logit_T)
  
  # cluster allocation
  S_dummy=cbind((X[,1]==1 & X[,2]==1),
                (X[,1]==0 & X[,3]==1),
                (X[,1]==0 & X[,3]==0 & X[,4]==1),
                (X[,1]==0 & X[,3]==0 & X[,4]==0))
  S_dummy=cbind(S_dummy,apply(S_dummy,1,sum)==0)
  S_cluster=S_dummy%*%(1:5)
  
  # outcomes 
  # normal-bivariate 
  # the covariance matrix is diagonal 
  Y=sapply(1:n, function(i) rmvnorm(1,c(eta_0[S_cluster[i]],eta_1[S_cluster[i]]),
                                    c(sigma_0[S_cluster[i]],sigma_1[S_cluster[i]])*diag(2)))
  
  # save all the information as a list
  return(list(data=list(X=X,T=T,Y=Y),                            # simulated data
              parameters=list(eta_0=eta_0, eta_1=eta_1,          # true parameters
                              sigma_0=sigma_0,sigma_1=sigma_1),  
              S_cluster=S_cluster))                              # cluster allocation
}


# --- generation data ---

scenario_1=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(2,4,6),
                            eta_1=c(0,3,6),
                            sigma_0=rep(0.3,3),
                            sigma_1=rep(0.3,3)))
scenario_2=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(0,2.2,4.4),
                            eta_1=c(0,0,0),
                            sigma_0=rep(0.2,3),
                            sigma_1=rep(0.2,3)))
scenario_3=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(1,2,3),
                            eta_1=c(0,1.5,3),
                            sigma_0=c(0.2,0.25,0.25),
                            sigma_1=c(0.25,0.3,0.2)))
scenario_4=lapply(1:samples, function(c) 
  simulation_sample_4groups(seed=c,
                            eta_0=c(1,2,3,3),
                            eta_1=c(0,1.5,3,4.5),
                            sigma_0=rep(0.2,4),
                            sigma_1=rep(0.2,4)))

# new scenarious
# 5 covariates - 5 groups
scenario_5=lapply(1:samples, function(c) 
  simulation_sample_5cov(seed=c,
                         eta_0=c(2,2,3,4.5,6.5),
                         eta_1=c(0,1,2.5,5,7.5),
                         sigma_0=rep(0.2,5),
                         sigma_1=rep(0.2,5)))

# 5 covariates - 5 groups - closer groups
scenario_6=lapply(1:samples, function(c) 
  simulation_sample_5cov(seed=c,
                         eta_0=c(2,2.5,3,3.5,3.5),
                         eta_1=c(1,1.75,2.5,3.5,4),
                         sigma_0=rep(0.2,5),
                         sigma_1=rep(0.2,5)))

# no heterogeneity - no groups
scenario_7=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=rep(2,3),
                            eta_1=rep(3,3),
                            sigma_0=rep(0.5,3),
                            sigma_1=rep(0.5,3)))

# no heterogeneity - 2 groups with same GATE
scenario_8=lapply(1:samples, function(c) 
  simulation_sample_3groups(seed=c,
                            eta_0=c(2,2,3),
                            eta_1=c(3,3,4),
                            sigma_0=rep(0.3,3),
                            sigma_1=rep(0.3,3)))



save.image("data_simulations.RData")

###################################################################################
# ---     plots: histograms of outcome distributions     ---
###################################################################################

setwd("C:/Users/dafne/Desktop/plot temporanei/BNP_for_HTE")

hist_sim<-function(data,scenario_n){
  
  # preparing dataframe for ggplot
  data_hist=as.data.frame(cbind(T=data$data$T,
                                Y=(data$data$T)*data$data$Y[2,]+(1-data$data$T)*data$data$Y[1,]))
  
  # histograms with ggpot2
  # plot is saved as a pdf
  pdf(file=paste0("sim_",scenario_n,".pdf"),width=8, height=6)
  g<-ggplot(data_hist, aes(Y, fill=as.factor(T))) + 
    geom_histogram(alpha = 0.4, position="identity")+
    scale_fill_manual(name = "Treatment \n level",
                      labels = c("0", "1"),
                      values = c("#E69F00", "#56B4E9"))+
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill ='transparent'),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent',
                                           colour='transparent'),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )+ 
    geom_vline(xintercept=data$parameters$eta_0, color = "#F2710C", 
               size=0.4)+
    geom_vline(xintercept=data$parameters$eta_1, color = "#4E60F5", 
               size=0.4,linetype = "dashed")+
    labs(x = "Outcome",y = " ")+
    ggtitle(paste0("    Scenario ",scenario_n))
  print(g)
  dev.off()
}

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
hist_sim(data=scenario_8[[1]],
         scenario_n=8)
