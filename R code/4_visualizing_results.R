#########################################################################
#     ---  MODEL COMPARISON   ----
#      --- tables and plots ----
#########################################################################

# load data
load("estimands.RData")

#libraries
library(xtable)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)

#########################################################################

cbPalette <- c("#007399", "#00ace6", "#80dfff")

bias_boxplot=as.data.frame(cbind(Xi=c(bias_ATE_CDBMM_1,bias_ATE_BART_1,bias_ATE_BCF_1,
                                      bias_ATE_CDBMM_2,bias_ATE_BART_2,bias_ATE_BCF_2,
                                      bias_ATE_CDBMM_3,bias_ATE_BART_3,bias_ATE_BCF_3,
                                      bias_ATE_CDBMM_4,bias_ATE_BART_4,bias_ATE_BCF_4,
                                      bias_ATE_CDBMM_5,bias_ATE_BART_5,bias_ATE_BCF_5,
                                      bias_ATE_CDBMM_6,bias_ATE_BART_6,bias_ATE_BCF_6,
                                      bias_ATE_CDBMM_7,bias_ATE_BART_7,bias_ATE_BCF_7),
                                 Q=(rep(c(rep("CDBMM",samples),rep("BART",samples),
                                          rep("BCF",samples)),7)),
                                 cov=paste0("scenario ",rep(1:7,each=samples*3))))
bias_boxplot$cov=as.character(bias_boxplot$cov)
bias_boxplot$Q=as.character(bias_boxplot$Q)
bias_boxplot$Xi=as.numeric(bias_boxplot$Xi)

pdf(file="bias_sim.pdf",width=10, height=5)
ggplot(bias_boxplot, aes(x=cov, y=Xi, fill=Q)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  geom_hline(yintercept = 0, col="#D90224", size=0.4) +
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        #panel.grid.minor = element_line(color = "grey"),
        axis.title = element_text(size=14),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("Bias") +
  xlab("")
dev.off()

mse_boxplot=as.data.frame(cbind(Xi=c(mse_ATE_CDBMM_1,mse_ATE_BART_1,mse_ATE_BCF_1,
                                     mse_ATE_CDBMM_2,mse_ATE_BART_2,mse_ATE_BCF_2,
                                     mse_ATE_CDBMM_3,mse_ATE_BART_3,mse_ATE_BCF_3,
                                     mse_ATE_CDBMM_4,mse_ATE_BART_4,mse_ATE_BCF_4,
                                     mse_ATE_CDBMM_5,mse_ATE_BART_5,mse_ATE_BCF_5,
                                     mse_ATE_CDBMM_6,mse_ATE_BART_6,mse_ATE_BCF_6,
                                     mse_ATE_CDBMM_7,mse_ATE_BART_7,mse_ATE_BCF_7),
                                Q=(rep(c(rep("CDBMM",samples),rep("BART",samples),
                                         rep("BCF",samples)),7)),
                                cov=paste0("scenario ",rep(1:7,each=samples*3))))
mse_boxplot$cov=as.character(mse_boxplot$cov)
mse_boxplot$Q=as.character(mse_boxplot$Q)
mse_boxplot$Xi=as.numeric(mse_boxplot$Xi)

pdf(file="mse_sim.pdf",width=10, height=5)
ggplot(mse_boxplot, aes(x=cov, y=Xi, fill=Q)) + 
  scale_fill_manual(values=cbPalette, name="")+
  geom_boxplot(lwd=0.3,fatten = 1.5, outlier.size = 0.3)+
  theme(panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill ="white"),
        axis.title = element_text(size=14),
        legend.text=element_text(size=10),
        plot.title = element_text(hjust = 0.2),
        title =element_text(size=18),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_line(color = "grey",size = 0.1))+
  ylab("MSE") +
  xlab("") +
  coord_cartesian(ylim = c(0.25, 1.2))
dev.off()


#########################################################################

#table rand index
rand_index<-cbind(matrix(apply(cbind(rand_CDBMM_1,rand_CDBMM_2,rand_CDBMM_3,rand_CDBMM_4,
                                     rand_CDBMM_5,rand_CDBMM_6,rand_CDBMM_7),
                               2,mean)),
                  matrix(apply(cbind(rand_CDBMM_1,rand_CDBMM_2,rand_CDBMM_3,rand_CDBMM_4,
                                     rand_CDBMM_5,rand_CDBMM_6,rand_CDBMM_7),
                               2,sd)),
                  matrix(apply(cbind(rand_CDBMM_1,rand_CART_2,rand_CART_3,rand_CART_4,
                                     rand_CART_5,rand_CART_6,rand_CART_7),
                               2,mean)),
                  matrix(apply(cbind(rand_CART_1,rand_CART_2,rand_CART_3,rand_CART_4,
                                     rand_CART_5,rand_CART_6,rand_CART_7),
                               2,sd)))

row.names(rand_index)<-paste0("scenario ",1:7)
colnames(rand_index)<-rep(c("mean","sd"),2)

xtable(rand_index, digits=4)