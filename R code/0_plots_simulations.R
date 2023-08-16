###################################################################################
# ---     plots: histograms of outcome distributions     ---
###################################################################################

load("data_simulations.RData")

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
