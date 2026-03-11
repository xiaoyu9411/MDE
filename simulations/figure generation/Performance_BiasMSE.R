library(gridExtra)
library(reshape2)

csn<-read.table("resultbetaS3.txt",sep=" ")

row<-nrow(csn)/25
csn<-csn[c(1:nrow(csn)),]
csn = structure(.Data=csn,.Dim=c(25, row))
csn<-t(csn)
csn<-as.data.frame(csn)
colnames(csn)<-c("N","rbmluo","rbsdshi",
                 "rbmqe","rbsdqe","rbmwqe","rbsdwqe",
                "rbmmdex","rbsdmdex",
                "msenluo","msesdshi","msemqe","msesdqe",
                "msemwqe","msesdwqe","msemmdex","msesdmdex",
                "msemsample","msesdsample","pqe","pwqe","pmdex",
                "disqe","diswqe","dismdex")


#relative bias_mean
rbmean<-csn[,c(1,2,4,6,8)]
colnames(rbmean)<-c("n","Luo","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Relative Bias")
p1 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = `Relative Bias`, group = Method, 
                colour = Method,linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c("dotted", "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "blue2", "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  coord_cartesian(ylim = c(-0.05, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Relative Bias in Mean") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14,  margin = margin(b = 15)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom",  # Move legend to bottom
    #legend.box.spacing = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),  # Subtle grid lines
    plot.margin = margin(15, 15, 15, 15)  # Add margin around plot
  ) +
  labs(
    x = "n",
    y = "Relative Bias"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p1


#relative bias_sd 
rbmean<-csn[,c(1,3,5,7,9)]
colnames(rbmean)<-c("n","Shi","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Relative Bias")
p2 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = `Relative Bias`, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c("dotted", "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "blue2", "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  coord_cartesian(ylim = c(-0.1, 0.1)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Relative Bias in SD") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14,  margin = margin(b = 15)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom",  # Move legend to bottom
    #legend.box.spacing = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),  # Subtle grid lines
    plot.margin = margin(15, 15, 15, 15)  # Add margin around plot
  ) +
  labs(
    x = "n",
    y = "Relative Bias"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p2



#RMSE_mean

rbmean<-csn[,c(1,10,12,14,16,18)]
colnames(rbmean)<-c("n","Luo","QE","weighted QE","MDEx","sample")
rbmean<-as.data.frame(cbind(rbmean$n,rbmean$Luo/rbmean$sample,rbmean$QE/rbmean$sample,
                            rbmean$`weighted QE`/rbmean$sample,rbmean$MDEx/rbmean$sample))
colnames(rbmean)<-c("n","Luo","QE","wQE","MDE")
#rbmean$MDE<-lnm
#rbmean$MDE[2:98]<-rbmean$MDE[2:98]+0.3
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","RMSE")
p3 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = RMSE, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c("dotted", "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "blue2", "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.025, 0.025)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("RMSE in Mean") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14,  margin = margin(b = 15)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom",  # Move legend to bottom
    #legend.box.spacing = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),  # Subtle grid lines
    plot.margin = margin(15, 15, 15, 15)  # Add margin around plot
  ) +
  labs(
    x = "n",
    y = "RMST"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p3


#RMSE_sd 
rbmean<-csn[,c(1,11,13,15,17,19)]
colnames(rbmean)<-c("n","Shi","QE","weighted QE","MDEx","sample")
rbmean<-as.data.frame(cbind(rbmean$n,rbmean$Shi/rbmean$sample,rbmean$QE/rbmean$sample,
                            rbmean$`weighted QE`/rbmean$sample,rbmean$MDEx/rbmean$sample))
colnames(rbmean)<-c("n","Shi","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","RMSE")

p4 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = RMSE, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c("dotted", "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "blue2", "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.025, 0.025)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("RMSE in SD") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14,  margin = margin(b = 15)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom",  # Move legend to bottom
    #legend.box.spacing = unit(0.3, "cm"),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_blank(),  # Subtle grid lines
    plot.margin = margin(15, 15, 15, 15)  # Add margin around plot
  ) +
  labs(
    x = "n",
    y = "RMSE"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p4


grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)



