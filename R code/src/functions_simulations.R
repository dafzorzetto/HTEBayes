#########################################################################
#        ---  SIMULATION FUNCTIONS    ----
#########################################################################

#library
library(mvtnorm)

# function with 2 COVARIATES and 3 GROUPS
simulation_sample_3groups<-function(seed,eta_0,eta_1,sigma_0,sigma_1,n){
  
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
simulation_sample_4groups<-function(seed,eta_0,eta_1,sigma_0,sigma_1,n){
  
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
simulation_sample_5cov<-function(seed,eta_0,eta_1,sigma_0,sigma_1,n){
  
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