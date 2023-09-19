#' @title
#' Visualization estimated GATEs 
#'
#' @description
#' Boxplots for the estimated GATEs
#'
#' @param GATEs : matrix of estimated GATEs
#' @param scenario_n : name of the considered scenario 
#' @param n_groups : number of groups
#'
#' @return
#' pdf file with the boxplot
#'
#' @import ggplot2


###################################################################################
# ---     plots: histograms of outcome distributions     ---
###################################################################################

#create the dataset for the boxplots
order_cluster<-function(partition,atoms, n_groups){
  
  clusters_0=length(atoms$p_0)
  clusters_1=length(atoms$p_1)
  groups=matrix(NA,clusters_0,clusters_1)
  units=matrix(NA,clusters_0,clusters_1)
  for (x in 1:clusters_0){
    for (y in 1:clusters_1) {
      groups[x,y]=atoms$p_1[y]-atoms$p_0[x]
      units[x,y]=sum(partition[,1]==x & partition[,2]==y)
    }
  }
  order_g=order(units, decreasing = TRUE)
  first_groups=(groups)[order_g][1:n_groups]
  
  return(sort(first_groups))
}

#libraries
library(ggplot2)

cbPalette <- c("#F4C000", "#D45C0B", "#EB002F", "#7C0BD4", "#0C51F7")

boxplots_gates<-function(GATEs, scenario_n, n_groups){
  
  boxplot_df=as.data.frame(cbind(Xi=c(GATEs),
                                 Q=(rep(1:n_groups,samples))))
  boxplot_df$Q=as.character(boxplot_df$Q)
  
  g<-ggplot(boxplot_df, aes(x=Q, y=Xi, fill=Q)) + 
    scale_fill_manual(values=cbPalette, name="groups")+
    geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
    theme(panel.background = element_rect(fill='white'),
          plot.background = element_rect(fill ="white"),
          axis.title = element_text(size=14),
          legend.text=element_text(size=10),
          plot.title = element_text(hjust = 0.2),
          title =element_text(size=12),
          legend.background = element_rect(fill='transparent'),
          panel.grid.major = element_line(color = "grey",size = 0.1))+
    ylab("GATEs") +
    xlab("") 
  
  pdf(file=paste0("GATEs_",scenario_n,".pdf"),width=7, height=5)
  print(g)
  dev.off()
  
}


