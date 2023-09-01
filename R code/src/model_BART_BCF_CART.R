#' @title
#' Estimation of BART, BCF, and BCF+CART combo
#'
#' @description
#' 3 function for the estimation of BART, BCF, and BCF+CART combo, respectively
#'
#' @param c : seed 
#' @param data_sample : list of data
#' @param estimated_Y (only for CART function): vector of estimated outcomes
#'
#' @return
#' For BART and BCF functions:
#' - *tau* : vector of individual treatment effect (ITE)
#' For CART function:
#' - *partition* : point estimation of group partition
#'
#' @import bartCause
#' @import bcf
#' @import rpart.plot
#' @import rattle

#########################################################################
# libraries
library(bartCause)
library(bcf)
library(rpart.plot)
library(rattle)

#########################################################################
#    ---     BART  ----
#########################################################################

# function to estimate the BART and to extrapolate the quantities of interest
bart_sample<-function(c,data_sample){
  
  # set seed for riproducibility
  set.seed(c)
  
  # ------   prearing variables   ------
  
  # main variables
  n=length(data_sample[[c]]$data$T)
  T_level=data_sample[[c]]$data$T             # treatment
  X=data_sample[[c]]$data$X                   # covariates-confounders
  Y_obs=unlist(sapply(1:n, function(i)        # observed outcome
    data_sample[[c]]$data$Y[(T_level[i]+1),i]))
  
  # ------   BART estimation   ------
  # estimation
  bart_fit=bartCause::bartc(as.matrix(Y_obs), as.matrix(T_level), as.matrix(X),
                            n.samples = 1000, n.burn = 1000)
  # espcted values for Y -> Y imputed
  Y_imp=matrix(c(apply(bartCause::extract(bart_fit, type = "y.0"),2,mean),
                 apply(bartCause::extract(bart_fit, type = "y.1"),2,mean)),ncol=2)
  
  E_Y_obs=apply(bart_fit$mu.hat.obs[1,,],2,mean)
  E_Y_cf=apply(bart_fit$mu.hat.cf[1,,],2,mean)
  
  tau=rep(NA,n)
  tau[T_level==0]=E_Y_cf[T_level==0]-E_Y_obs[T_level==0]
  tau[T_level==1]=E_Y_obs[T_level==1]-E_Y_cf[T_level==1]
  
  #return(list(Y_mis=sapply(1:n, function(i) Y_imp[i,(2-T_level[i])]),
  #            tau=Y_imp[,2]-Y_imp[,1]))
  return(list(tau=tau))
}

#########################################################################
#    ---    BCF    ----
#########################################################################

# function to estimate the BCF and to extrapolate the quantities of interest
BCF_sample<-function(c,data_sample){
  
  # set seed for riproducibility
  set.seed(c)
  
  # ------   prearing variables   ------
  
  # main variables
  n=length(data_sample[[c]]$data$T)
  T_level=data_sample[[c]]$data$T             # treatment
  X=as.matrix(data_sample[[c]]$data$X)        # covariates-confounders
  Y_obs=unlist(sapply(1:n, function(i)        # observed outcome
    data_sample[[c]]$data$Y[(T_level[i]+1),i]))
  
  # ------   BCF estimation   ------
  
  #propensity score
  p.score <- glm(T_level ~ X,
                 family = binomial,
                 data = as.data.frame(cbind(T_level, X)))
  pihat <- predict(p.score, as.data.frame(X))
  
  # estimation
  bcf_tau <- bcf(Y_obs, T_level, X, X, pihat, 
                 nburn = 1000, nsim = 1000)
  
  return(list(tau=apply(bcf_tau$tau,2,mean)))
}

#########################################################################
#    ---  Clustering with CART    ----
#########################################################################

# function to estimate the BART with BCF results
CART<-function(c, data_sample, estimated_Y){
  
  # set seed for riproducibility
  set.seed(c)
  
  # ------   prearing variables   ------
  
  # main variables
  T_level=data_sample[[c]]$data$T             # treatment
  X=data_sample[[c]]$data$X                   # covariates-confounders
  tau=estimated_Y[[c]]$tau                    # estimation of individual treatment effect
  
  variables= as.data.frame(cbind(tau, X))
  
  # CART
  fit.tree <- rpart(tau ~ ., data= variables, cp=0.01, minsplit=10)
  # extrapolation of group allocation 
  GroupsAllocation_cart <- c(fit.tree$where)
  
  return(list(partition=GroupsAllocation_cart))
}
