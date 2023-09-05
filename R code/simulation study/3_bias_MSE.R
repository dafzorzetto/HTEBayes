#########################################################################
#        ---  MODEL COMPARISON    ----
#
#   ---- ATE among CDBMM, BART, and BCF   ----
#     ---- GATE between CDBMM and CART ----
#########################################################################

# load data
load("results_models.RData")

#liraries
library(mcclust)

#########################################################################

# simulated individual treatment effect (ITE)
simulated_tau_1=sapply(1:samples, function(s) scenario_1[[s]]$data$Y[2,]-scenario_1[[s]]$data$Y[1,])
simulated_tau_2=sapply(1:samples, function(s) scenario_2[[s]]$data$Y[2,]-scenario_2[[s]]$data$Y[1,])
simulated_tau_3=sapply(1:samples, function(s) scenario_3[[s]]$data$Y[2,]-scenario_3[[s]]$data$Y[1,])
simulated_tau_4=sapply(1:samples, function(s) scenario_4[[s]]$data$Y[2,]-scenario_4[[s]]$data$Y[1,])
simulated_tau_5=sapply(1:samples, function(s) scenario_5[[s]]$data$Y[2,]-scenario_5[[s]]$data$Y[1,])
simulated_tau_6=sapply(1:samples, function(s) scenario_6[[s]]$data$Y[2,]-scenario_6[[s]]$data$Y[1,])
simulated_tau_7=sapply(1:samples, function(s) scenario_7[[s]]$data$Y[2,]-scenario_7[[s]]$data$Y[1,])

#########################################################################
# --- bias and MSE for Average Treatment Effect (ATE) ---
#########################################################################

# ATE computed by CDBMM
bias_ATE_CDBMM_1=apply(sapply(1:samples, function(s) CDBMM_scenario_1[[s]]$tau-simulated_tau_1[,s]),2,mean)
bias_ATE_CDBMM_2=apply(sapply(1:samples, function(s) CDBMM_scenario_2[[s]]$tau-simulated_tau_2[,s]),2,mean)
bias_ATE_CDBMM_3=apply(sapply(1:samples, function(s) CDBMM_scenario_3[[s]]$tau-simulated_tau_3[,s]),2,mean)
bias_ATE_CDBMM_4=apply(sapply(1:samples, function(s) CDBMM_scenario_4[[s]]$tau-simulated_tau_4[,s]),2,mean)
bias_ATE_CDBMM_5=apply(sapply(1:samples, function(s) CDBMM_scenario_5[[s]]$tau-simulated_tau_5[,s]),2,mean)
bias_ATE_CDBMM_6=apply(sapply(1:samples, function(s) CDBMM_scenario_6[[s]]$tau-simulated_tau_6[,s]),2,mean)
bias_ATE_CDBMM_7=apply(sapply(1:samples, function(s) CDBMM_scenario_7[[s]]$tau-simulated_tau_7[,s]),2,mean)

mse_ATE_CDBMM_1=apply(sapply(1:samples, function(s) (CDBMM_scenario_1[[s]]$tau-simulated_tau_1[,s])^2),2,mean)
mse_ATE_CDBMM_2=apply(sapply(1:samples, function(s) (CDBMM_scenario_2[[s]]$tau-simulated_tau_2[,s])^2),2,mean)
mse_ATE_CDBMM_3=apply(sapply(1:samples, function(s) (CDBMM_scenario_3[[s]]$tau-simulated_tau_3[,s])^2),2,mean)
mse_ATE_CDBMM_4=apply(sapply(1:samples, function(s) (CDBMM_scenario_4[[s]]$tau-simulated_tau_4[,s])^2),2,mean)
mse_ATE_CDBMM_5=apply(sapply(1:samples, function(s) (CDBMM_scenario_5[[s]]$tau-simulated_tau_5[,s])^2),2,mean)
mse_ATE_CDBMM_6=apply(sapply(1:samples, function(s) (CDBMM_scenario_6[[s]]$tau-simulated_tau_6[,s])^2),2,mean)
mse_ATE_CDBMM_7=apply(sapply(1:samples, function(s) (CDBMM_scenario_7[[s]]$tau-simulated_tau_7[,s])^2),2,mean)

# ATE computed by BART
bias_ATE_BART_1=apply(sapply(1:samples, function(s) BART_scenario_1[[s]]$tau-simulated_tau_1[,s]),2,mean)
bias_ATE_BART_2=apply(sapply(1:samples, function(s) BART_scenario_2[[s]]$tau-simulated_tau_2[,s]),2,mean)
bias_ATE_BART_3=apply(sapply(1:samples, function(s) BART_scenario_3[[s]]$tau-simulated_tau_3[,s]),2,mean)
bias_ATE_BART_4=apply(sapply(1:samples, function(s) BART_scenario_4[[s]]$tau-simulated_tau_4[,s]),2,mean)
bias_ATE_BART_5=apply(sapply(1:samples, function(s) BART_scenario_5[[s]]$tau-simulated_tau_5[,s]),2,mean)
bias_ATE_BART_6=apply(sapply(1:samples, function(s) BART_scenario_6[[s]]$tau-simulated_tau_6[,s]),2,mean)
bias_ATE_BART_7=apply(sapply(1:samples, function(s) BART_scenario_7[[s]]$tau-simulated_tau_7[,s]),2,mean)

mse_ATE_BART_1=apply(sapply(1:samples, function(s) (BART_scenario_1[[s]]$tau-simulated_tau_1[,s])^2),2,mean)
mse_ATE_BART_2=apply(sapply(1:samples, function(s) (BART_scenario_2[[s]]$tau-simulated_tau_2[,s])^2),2,mean)
mse_ATE_BART_3=apply(sapply(1:samples, function(s) (BART_scenario_3[[s]]$tau-simulated_tau_3[,s])^2),2,mean)
mse_ATE_BART_4=apply(sapply(1:samples, function(s) (BART_scenario_4[[s]]$tau-simulated_tau_4[,s])^2),2,mean)
mse_ATE_BART_5=apply(sapply(1:samples, function(s) (BART_scenario_5[[s]]$tau-simulated_tau_5[,s])^2),2,mean)
mse_ATE_BART_6=apply(sapply(1:samples, function(s) (BART_scenario_6[[s]]$tau-simulated_tau_6[,s])^2),2,mean)
mse_ATE_BART_7=apply(sapply(1:samples, function(s) (BART_scenario_7[[s]]$tau-simulated_tau_7[,s])^2),2,mean)

# ATE computed by BCF
bias_ATE_BCF_1=apply(sapply(1:samples, function(s) BCF_scenario_1[[s]]$tau-simulated_tau_1[,s]),2,mean)
bias_ATE_BCF_2=apply(sapply(1:samples, function(s) BCF_scenario_2[[s]]$tau-simulated_tau_2[,s]),2,mean)
bias_ATE_BCF_3=apply(sapply(1:samples, function(s) BCF_scenario_3[[s]]$tau-simulated_tau_3[,s]),2,mean)
bias_ATE_BCF_4=apply(sapply(1:samples, function(s) BCF_scenario_4[[s]]$tau-simulated_tau_4[,s]),2,mean)
bias_ATE_BCF_5=apply(sapply(1:samples, function(s) BCF_scenario_5[[s]]$tau-simulated_tau_5[,s]),2,mean)
bias_ATE_BCF_6=apply(sapply(1:samples, function(s) BCF_scenario_6[[s]]$tau-simulated_tau_6[,s]),2,mean)
bias_ATE_BCF_7=apply(sapply(1:samples, function(s) BCF_scenario_7[[s]]$tau-simulated_tau_7[,s]),2,mean)

mse_ATE_BCF_1=apply(sapply(1:samples, function(s) (BCF_scenario_1[[s]]$tau-simulated_tau_1[,s])^2),2,mean)
mse_ATE_BCF_2=apply(sapply(1:samples, function(s) (BCF_scenario_2[[s]]$tau-simulated_tau_2[,s])^2),2,mean)
mse_ATE_BCF_3=apply(sapply(1:samples, function(s) (BCF_scenario_3[[s]]$tau-simulated_tau_3[,s])^2),2,mean)
mse_ATE_BCF_4=apply(sapply(1:samples, function(s) (BCF_scenario_4[[s]]$tau-simulated_tau_4[,s])^2),2,mean)
mse_ATE_BCF_5=apply(sapply(1:samples, function(s) (BCF_scenario_5[[s]]$tau-simulated_tau_5[,s])^2),2,mean)
mse_ATE_BCF_6=apply(sapply(1:samples, function(s) (BCF_scenario_6[[s]]$tau-simulated_tau_6[,s])^2),2,mean)
mse_ATE_BCF_7=apply(sapply(1:samples, function(s) (BCF_scenario_7[[s]]$tau-simulated_tau_7[,s])^2),2,mean)

#########################################################################
# --- Group Average Treatment Effect (GATE) ---
#           --- CDBMM vs CART ---
#########################################################################

# RAND index for clusters obtained with CDBMM
rand_CDBMM_1=sapply(1:samples, function(s) arandi(scenario_1[[s]]$S_cluster,CDBMM_scenario_1[[s]]$partitio[,1]*(CDBMM_scenario_1[[s]]$partitio[,2]+L_0)))
rand_CDBMM_2=sapply(1:samples, function(s) arandi(scenario_2[[s]]$S_cluster,CDBMM_scenario_2[[s]]$partitio[,1]*(CDBMM_scenario_2[[s]]$partitio[,2]+L_0)))
rand_CDBMM_3=sapply(1:samples, function(s) arandi(scenario_3[[s]]$S_cluster,CDBMM_scenario_3[[s]]$partitio[,1]*(CDBMM_scenario_3[[s]]$partitio[,2]+L_0)))
rand_CDBMM_4=sapply(1:samples, function(s) arandi(scenario_4[[s]]$S_cluster,CDBMM_scenario_4[[s]]$partitio[,1]*(CDBMM_scenario_4[[s]]$partitio[,2]+L_0)))
rand_CDBMM_5=sapply(1:samples, function(s) arandi(scenario_5[[s]]$S_cluster,CDBMM_scenario_5[[s]]$partitio[,1]*(CDBMM_scenario_5[[s]]$partitio[,2]+L_0)))
rand_CDBMM_6=sapply(1:samples, function(s) arandi(scenario_6[[s]]$S_cluster,CDBMM_scenario_6[[s]]$partitio[,1]*(CDBMM_scenario_6[[s]]$partitio[,2]+L_0)))
rand_CDBMM_7=sapply(1:samples, function(s) arandi(scenario_7[[s]]$S_cluster,rep(1,n)))

# RAND index for clusters obtained with CART
rand_CART_1=sapply(1:samples, function(s) arandi(scenario_1[[s]]$S_cluster,CART_scenario_1[[s]]$partition))
rand_CART_2=sapply(1:samples, function(s) arandi(scenario_2[[s]]$S_cluster,CART_scenario_2[[s]]$partition))
rand_CART_3=sapply(1:samples, function(s) arandi(scenario_3[[s]]$S_cluster,CART_scenario_3[[s]]$partition))
rand_CART_4=sapply(1:samples, function(s) arandi(scenario_4[[s]]$S_cluster,CART_scenario_4[[s]]$partition))
rand_CART_5=sapply(1:samples, function(s) arandi(scenario_5[[s]]$S_cluster,CART_scenario_5[[s]]$partition))
rand_CART_6=sapply(1:samples, function(s) arandi(scenario_6[[s]]$S_cluster,CART_scenario_6[[s]]$partition))
rand_CART_7=sapply(1:samples, function(s) arandi(scenario_7[[s]]$S_cluster,rep(1,n)))


#########################################################################
save.image("estimands.RData")
