csn<-read.table("resultbetaS1.txt",sep=" ")

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

rbmean<-csn[,c(1,23,24,25)]
colnames(rbmean)<-c("n","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Proportion")

p1 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = Proportion, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c( "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.1, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Beta distribution") +
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
    y = "Proportion"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p1



csn<-read.table("resultnormalS1.txt",sep=" ")

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

rbmean<-csn[,c(1,23,24,25)]
colnames(rbmean)<-c("n","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Proportion")

p2 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = Proportion, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c( "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.1, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Normal distribution") +
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
    y = "Proportion"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p2



csn<-read.table("resultlognormalS1.txt",sep=" ")

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

rbmean<-csn[,c(1,23,24,25)]
colnames(rbmean)<-c("n","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Proportion")

p3 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = Proportion, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c( "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.1, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Lognormal distribution") +
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
    y = "Proportion"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p3


csn<-read.table("resultgammaS1.txt",sep=" ")

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

rbmean<-csn[,c(1,23,24,25)]
colnames(rbmean)<-c("n","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Proportion")

p4 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = Proportion, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c( "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.1, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Gamma distribution") +
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
    y = "Proportion"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p4

csn<-read.table("resultweibullS1.txt",sep=" ")

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

rbmean<-csn[,c(1,23,24,25)]
colnames(rbmean)<-c("n","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Proportion")

p5 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = Proportion, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c( "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.1, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Weibull distribution") +
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
    y = "Proportion"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p5




csn<-read.table("resultweibullS1o.txt",sep=" ")

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

rbmean<-csn[,c(1,23,24,25)]
colnames(rbmean)<-c("n","QE","wQE","MDE")
rbmean<-melt(rbmean,id="n")
colnames(rbmean)<-c("n","Method","Proportion")

p6 <- ggplot() +
  geom_line(data = rbmean, 
            aes(x = n, y = Proportion, group = Method, 
                colour = Method, linetype = Method),
            size = 0.9,  # Thicker lines for better visibility
            alpha = 0.8) +  # Slight transparency
  scale_linetype_manual(values = c( "longdash", "solid", "dotdash")) +
  scale_colour_manual(values = c( "#4DAF4A","#E41A1C", "black")) +  # Colorblind-friendly colors
  #coord_cartesian(ylim = c(-0.1, 0.05)) +  # Better than ylim() as it doesn't remove data
  theme_classic() +
  ggtitle("Exponential distribution") +
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
    y = "Proportion"
  ) +
  guides(
    colour = guide_legend(nrow = 1, title.position = "top"),  # Single row legend
    linetype = guide_legend(nrow = 1, title.position = "top")
  )

p6



grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
