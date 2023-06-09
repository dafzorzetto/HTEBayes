#########################################################################
#        ---  ESTIMATING CDBMM    ----
#########################################################################

# libraries
library(mvtnorm)
library(CholWishart)
library(parallel)
library(truncnorm)
library(invgamma)
library(BNPmix)

# load simulated data
load("simulazione_paper_2.RData")

# iterations
R=5000
R_burnin=2000

# number of maximum groups for the two treatment levels
L_0=15
L_1=15


#########################################################################
#    ---  Confonder-Dependent Bayesian Mixture Model    ----
#########################################################################

PSB_var<-function(c,data_sample){
  
  # set seed for riproducibility
  set.seed(c)
  
  # ------   prearing variables   ------
  
  # main variables
  T_level=data_sample[[c]]$data$T             # treatment
  X=data_sample[[c]]$data$TX                  # covariates-confounders
  Y_obs=unlist(sapply(1:n, function(i)        # observed outcome
    data_sample[[c]]$data$Y[(T_level+1),i]))
  
  # number covariates (X)
  n_X=dim(X)[2]+1
  
  # dividing observation in the treatment levels
  # level t=1
  T1=which(T_level==1)  
  n1=length(T1)
  Y_obs_1=Y_obs[T1]
  # level t=0
  T0=which(T_level==0)
  n0=length(T0)
  Y_obs_0=Y_obs[T0]
  
  # ------   hyperparameters   -----
  p_beta=c(0,3)                 # mean and variance of normal distribution for beta param.
  p_sigma=c(5,1)                # shape parameters of inv-gamma for sigma param.
  p_eta=c(0,2)                  # mean and variance of normal distribution for eta param.
  
  # ------   initialitation   -----
  # parameters
  sigma_0=rep(1,L_0)
  sigma_1=rep(1,L_1)
  beta_0=rep(0,n_X*(L_0-1))
  beta_1=rep(0,n_X*(L_1-1))
  eta_0=seq(-3,-3+0.5*(L_0-1),0.5)
  eta_1=seq(-3,-3+0.5*(L_1-1),0.5)
  # cluster allocation variables
  xi_0=sample(1:L_0,n0,replace=TRUE)
  xi_1=sample(1:L_1,n1,replace=TRUE)
  
  # ------   useful quantities   ------   
  # number of units in each cluster
  units_for_clusters_0=table(xi_0)
  units_for_clusters_1=table(xi_1)
  # eta corresponding to each unit
  eta_sogg_0=eta_0[xi_0]
  eta_sogg_1=eta_1[xi_1]
  # latent variable for data augmentation in probit regression 
  Z_0=rep(list(rep(NA,L_0)),n0)
  Z_1=rep(list(rep(NA,L_1)),n1)
  
  # ------   useful functions   -----
  
  # probit stick-breaking weights 
  omega=function(X,beta){
    L=length(beta)/n_X
    alpha=c(sapply(1:L, function(l) pnorm(beta[(l*n_X-(n_X-1)):(l*n_X)]%*%c(1,X))),1)
    return(c(alpha[1],sapply(2:(L+1), function(l) alpha[l]*prod(1-alpha[1:(l-1)]))))
  }
  
  # point estimate parition
  estimation_partition<-function(xi_0,xi_1){
    
    # using Variation of Information as loss function
    clusters_0=partition.BNPdens(list(clust=t(xi_0)),dist = "VI")$partitions[1,]
    clusters_1=partition.BNPdens(list(clust=t(xi_1)),dist = "VI")$partitions[1,]
    
    # alternativly: Binder's loss
    #clusters_0=partition.BNPdens(list(clust=t(xi_0)),dist = "Binder")$partitions
    #clusters_1=partition.BNPdens(list(clust=t(xi_1)),dist = "Binder")$partitions
    
    clusters=cbind(clusters_0,clusters_1)
    dimnames(clusters)=NULL
    return(clusters)
  }
  
  # funtion to impute missing outcomes
  posterior_estimation_ymis<-function(partizione){
    
    # create empty list
    stime=list(stime_0=rep(NA,max(partizione[,1])),
               stime_1=rep(NA,max(partizione[,2])))
    
    # imputing missing outcome
    for (g in 1:max(partizione[,1])){
      unita=intersect(T0,which(partizione[,1]==g))
      stime$stime_0[g]=mean(sapply(1:length(unita), function(x) 
        mean(sapply(1:(R-R_burnin), function(j) 
          post_eta[cluster_allocation_0[unita[x],j],R_burnin+j]))))
    }
    for (g in 1:max(partizione[,2])){
      unita=intersect(T1,which(partizione[,2]==g))
      stime$stime_1[g]=mean(sapply(1:length(unita), function(x) 
        mean(sapply(1:(R-R_burnin), function(j) 
          post_eta[L_0+cluster_allocation_1[unita[x],j],R_burnin+j]))))
    }
    
    y_mis_post_mean=sapply(1:n, function(x) 
      ifelse(t[x]==0, 
             stime$stime_1[partizione[x,2]],
             stime$stime_0[partizione[x,1]]))
    return(y_mis_post_mean)
  }
  
  # ------   saving informations   -----
  
  # empty matrix where save all the informations for each iteration
  post_eta=matrix(NA,ncol=R,nrow=L_0+L_1)
  post_var=matrix(NA,ncol=R,nrow=L_0+L_1)
  cluster_allocation_0=matrix(NA,ncol=R-R_burnin,nrow=n)
  cluster_allocation_1=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_beta=matrix(NA,ncol=R, nrow=length(beta_0)+length(beta_1))
  
  
  # -----   updating parameters and variables at each itearation   -----
  
  for (r in 1:R){
    
    # -----   ETA: mean of normal distributions (atoms)   -----
    # eta_0 = eta for treatment level 0
    eta_0[1]=rnorm(1,mean=1/(units_for_clusters_0[1]/sigma_0[1]+1/p_eta[2])*(sum(Y_obs_0[xi_0==1])/sigma_0[1]+p_eta[1]/p_eta[2]),
                   sd=sqrt(1/(units_for_clusters_0[1]/sigma_0[1]+1/p_eta[2])))
    for(l in 2:L_0){
      eta_0[l]=rtruncnorm(1,a=eta_0[l-1],mean=1/(units_for_clusters_0[l]/sigma_0[l]+1/p_eta[2])*(sum(Y_obs_0[xi_0==l])/sigma_0[l]+p_eta[1]/p_eta[2]),
                          sd=sqrt(1/(units_for_clusters_0[l]/sigma_0[l]+1/p_eta[2])))
    }
    # eta_1 = eta for treatment level 1
    eta_1[1]=rnorm(1,mean=1/(units_for_clusters_1[1]/sigma_1[1]+1/p_eta[2])*(sum(Y_obs_1[xi_1==1])/sigma_1[1]+p_eta[1]/p_eta[2]),
                   sd=sqrt(1/(units_for_clusters_1[1]/sigma_1[1]+1/p_eta[2])))
    for(l in 2:L_1){
      eta_1[l]=rtruncnorm(1,a=eta_1[l-1],mean=1/(units_for_clusters_1[l]/sigma_1[l]+1/p_eta[2])*(sum(Y_obs_1[xi_1==l])/sigma_1[l]+p_eta[1]/p_eta[2]),
                          sd=sqrt(1/(units_for_clusters_1[l]/sigma_1[l]+1/p_eta[2])))
    }
    
    eta_sogg_0=eta_0[xi_0]
    eta_sogg_1=eta_1[xi_1]
    
    # -----   SIGMA: variance of normal distributions (atoms)   -----
    # sigma_0 = sigma for treatment level 0
    sigma_0=sapply(1:L_0, function(l) rinvgamma(1,p_sigma[1]+sum(xi_0==l)/2,
                                                 p_sigma[2]+sum((Y_obs_0[xi_0==l]-eta_0[l])^2)/2))
    # sigma_1 = sigma for treatment level 1
    sigma_1=sapply(1:L_0, function(l) rinvgamma(1,p_sigma[1]+sum(xi_1==l)/2,
                                                 p_sigma[2]+sum((Y_obs_1[xi_1==l]-eta_1[l])^2)/2))
    
    #OMEGA:
    #omega_0
    omega_0=sapply(1:n0, function(i) omega(X=X[T0[i],],beta=beta_0))
    #omega_1
    omega_1=sapply(1:n1, function(i) omega(X=X[T1[i],],beta=beta_1))
    
    #Xi:
    #Xi_0
    dmn_0=sapply(1:n0, function(i) sapply(1:L_0, function(l) dnorm(Y_obs_0[i], eta_0[l], sqrt(sigma_0[l]), log=TRUE)))+log(omega_0)
    dmn_0[which(is.nan(dmn_0))]=-100
    xi=sapply(1:n0, function(i) rmultinom(1,1,exp(dmn_0[,i])))
    xi_0=sapply(1:n0, function(i) xi[,i]%*%(1:L_0))
    units_for_clusters_0=apply(xi, 1, sum)
    
    #Xi_1
    dmn_1=sapply(1:n1, function(i) sapply(1:L_1, function(l) dnorm(Y_obs_1[i], eta_1[l], sqrt(sigma_1[l]), log=TRUE)))+log(omega_1)
    dmn_1[which(is.nan(dmn_1))]=-100
    xi=sapply(1:n1, function(i) rmultinom(1,1,exp(dmn_1[,i])))
    xi_1=sapply(1:n1, function(i) xi[,i]%*%(1:L_1))
    units_for_clusters_1=apply(xi, 1, sum)
    
    # Z
    #pesi_0=exp(dmn_0)
    #somma_0=apply(pesi_0, 2, sum)
    #pesi_0=t(pesi_0)/somma_0
    pesi_0=t(omega_0)
    mu_z=cbind((pesi_0[,1]),(pesi_0[,2]/(1-pesi_0[,1])))
    if (L_0>3){
      mu_z=cbind(mu_z,sapply(3:(L_0-1), function(l) (pesi_0[,l]/(1-apply(pesi_0[,1:(l-1)],1,sum)))))
    }
    mu_z[which(is.nan(mu_z))]=1
    mu_z[which(mu_z>1)]=1
    mu_z=mu_z-9.9e-15*(mu_z>(1-1e-16))
    #Z_0:
    for (i in 1:n0){
      for (l in 1:(min(xi_0[i],L_0-1))) {
        if(l>1){
          if (l<xi_0[i]){
            Z_0[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_0[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }else{
          if (l<xi_0[i]){
            Z_0[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_0[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }
      }
    }
    #Z_1:
    #pesi_1=exp(dmn_1)
    #somma_1=apply(pesi_1, 2, sum)
    #pesi_1=t(pesi_1)/somma_1
    pesi_1=t(omega_1)
    mu_z=cbind((pesi_1[,1]),(pesi_1[,2]/(1-pesi_1[,1])))
    if (L_1>3){
      mu_z=cbind(mu_z,sapply(3:(L_1-1), function(l) (pesi_1[,l]/(1-apply(pesi_1[,1:(l-1)],1,sum)))))
    }
    mu_z[which(is.nan(mu_z))]=1
    mu_z[which(mu_z>1)]=1
    mu_z=mu_z-9.9e-15*(mu_z>(1-1e-16))
    for (i in 1:n1){
      for (l in 1:(min(xi_1[i],L_1-1))) {
        if(l>1){
          if (l<xi_1[i]){
            Z_1[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_1[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }else{
          if (l<xi_1[i]){
            Z_1[[i]][l]=rtruncnorm(1,b=0,mean=qnorm(mu_z[i,l]))
          }else{
            Z_1[[i]][l]=rtruncnorm(1,a=0,mean=qnorm(mu_z[i,l]))
          }
        }
      }
    }
    
    #BETA:
    #beta_0
    gruppi=which(units_for_clusters_0!=0)
    if (max(gruppi)==L_0)
      gruppi=gruppi[-length(gruppi)]
    for (l in gruppi){
      val=which(xi_0>=l)
      z_tilde=unlist(sapply(val, function(i) Z_0[[i]][l]))
      x_tilde=cbind(rep(1,length(val)),matrix(X[T0[val],],ncol=n_X-1))
      V=solve(diag(n_X)/p_beta[2]+t(x_tilde)%*%x_tilde)
      beta_0[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,V%*%(1/p_beta[2]*(diag(n_X)%*%rep(p_beta[1],n_X))+t(x_tilde)%*%z_tilde),V)[,1:n_X]
    }
    vuoti=which(units_for_clusters_0==0)
    if (length(vuoti)>0){
      if (max(vuoti)==L_1)
        vuoti=vuoti[-length(vuoti)]
    }
    for (l in vuoti){
      beta_0[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,rep(p_beta[1],n_X),diag(n_X)/p_beta[2])[,1:n_X]
    }
    #beta_1
    gruppi=which(units_for_clusters_1!=0)
    if (max(gruppi)==L_1)
      gruppi=gruppi[-length(gruppi)]
    for (l in gruppi){
      val=which(xi_1>=l)
      z_tilde=unlist(sapply(val, function(i) Z_1[[i]][l]))
      x_tilde=cbind(rep(1,length(val)),matrix(X[T1[val],],ncol=n_X-1))
      V=solve(diag(n_X)/p_beta[2]+t(x_tilde)%*%x_tilde)
      beta_1[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,V%*%(1/p_beta[2]*(diag(n_X)%*%rep(p_beta[1],n_X))+t(x_tilde)%*%z_tilde),V)[,1:n_X]
    }
    vuoti=which(units_for_clusters_1==0)
    if (length(vuoti)>0){
      if (max(vuoti)==L_1)
        vuoti=vuoti[-length(vuoti)]
    }
    for (l in vuoti){
      beta_1[(l*n_X-(n_X-1)):(l*n_X)]=rmvnorm(1,rep(p_beta[1],n_X),diag(n_X)/p_beta[2])[,1:n_X]
    }
    
    # -----   missing outcome   -----
    
    # estimate cluster allocation for eqch iteration
    if(r>R_burnin){
      
      # level t=0 observed --> level T=1 missing
      om_0=sapply(T0, function(i) omega(X=X[i,],beta=beta_1))
      om_0[which(is.nan(om_0))]=0
      clusters_Ymiss_0=sapply(1:n0, function(i) (1:L_1)%*%(rmultinom(1,1,om_0[,i])))
      
      # level t=1 observed --> level T=0 missing
      om_1=sapply(T1, function(i) omega(X=X[i,],beta=beta_0))
      om_1[which(is.nan(om_1))]=0
      clusters_Ymiss_1=sapply(1:n1, function(i) (1:L_0)%*%(rmultinom(1,1,om_1[,i])))
    }
    
    # -----   saving information   -----
    # parameters
    post_eta[,r]=c(eta_0,eta_1)
    post_var[,r]=c(sigma_0,sigma_1)
    post_beta[,r]=c(beta_0,beta_1)
    # cluster allocation for each iteration
    if(r>R_burnin){
      cluster_allocation_0[T0,r-R_burnin]=xi_0
      cluster_allocation_0[T1,r-R_burnin]=clusters_Ymiss_1
      cluster_allocation_1[T1,r-R_burnin]=xi_1
      cluster_allocation_1[T0,r-R_burnin]=clusters_Ymiss_0
    }
    
    #
    if (r%%500==0) print(r)
  }
  
  # -----   point estimation partition   -----
  # partition: cluster allocation
  partition=estimation_partition(cluster_allocation_0,cluster_allocation_1)
  # posterior mean for missing outcome
  Y_mis=posterior_estimation_ymis(partition)

  
  print(paste0("campione n ",c))
  
  return(list(post_eta=post_eta,
              post_var=post_var,
              post_beta=post_beta,
              partition=partition,
              Y_mis=Y_mis))
}


