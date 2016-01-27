library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library(scales)

data=read.table("/Users/jpiscionere/Boss/AbsMag_Bin12_Bin34_LRG",head=T)
data$MagR[which(data$Bin=="bin34")] = data$MagR[which(data$Bin=="bin34")] - 0.8
frame=data.frame(data)

p<-ggplot(data,aes(x=Z,y=MagR)) + geom_density2d(aes(colour=Bin),alpha=1/1000) + geom_point(alpha=1/150)
p<-p + scale_color_discrete("Galaxy Sample",breaks=c("aLRG","bin12","bin34"),labels=c("Luminous Red Galaxies","CMASS 0.43 < z < 0.55","CMASS 0.55 < z < 0.7"))
p<-p + theme_bw() + ggtitle("Redshift vs Absolute Magnitude") + xlab("Redshift") + ylab("Absolute Magnitude")
p<-p + scale_y_reverse(limits=c(-20,-23)) + theme(legend.position=c(0.25,0.9),legend.key.size = unit(0.7, "cm"),legend.text=element_text(size=15),legend.title=element_text(size=15))

p<-p + annotate("rect",xmin=0.16,xmax=0.36,ymin=-23,ymax=-20,alpha=0.2,fill="red")
p<-p + annotate("rect",xmin=0.43,xmax=0.55,ymin=-23,ymax=-20,alpha=0.2,fill="green")
p<-p + annotate("rect",xmin=0.55,xmax=0.7,ymin=-23,ymax=-20,alpha=0.2,fill="blue")
p<-p + guides(colour=guide_legend(override.aes=list(alpha=1)))
p<-p + theme(axis.title.y=element_text(size=16,angle=90,vjust=-0.2),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=15),panel.border=element_rect(size=1))

