library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')
library(scales)

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/AbsMag_Bin12_Bin34_LRG",head=T)
data$MagR[which(data$Bin=="bin34")] = data$MagR[which(data$Bin=="bin34")] - 0.8
frame=data.frame(data)

p<-ggplot(data,aes(x=Z,y=MagR)) + geom_density2d(aes(colour=Bin)) + geom_point(alpha=1/100)
p<-p + scale_color_discrete("Galaxy Sample",breaks=c("LRG","bin12","bin34"),labels=c("Luminous Red Galaxies","CMASS 0.43 < z < 0.55","CMASS 0.55 < z < 0.7"))
p<-p + theme_bw() + ggtitle("Redshift vs Absolute Magnitude") + xlab("Redshift") + ylab("Absolute Magnitude")
p<-p + scale_y_reverse(limits=c(-20,-23)) + theme(legend.position=c(0.25,0.9),legend.key.size = unit(0.7, "cm"),legend.text=element_text(size=15),legend.title=element_text(size=15))
