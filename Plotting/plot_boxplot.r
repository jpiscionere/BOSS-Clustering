library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/Redshift_Magnitude.out")







redshift=data$V3[which(data$V3 >= 0.43 & data$V3 <= 0.7)]
magnitude=data$V4[which(data$V3 >= 0.43 & data$V3 <= 0.7)]
bin=c(1:length(magnitude))
bin=bin*0
bin[which(redshift >= 0.43 & redshift < 0.5)]="Bin1"
bin[which(redshift >= 0.5 & redshift < 0.55)]="Bin2"
bin[which(redshift >= 0.55 & redshift < 0.6)]="Bin3"
bin[which(redshift >= 0.6 & redshift < 0.7)]="Bin4"

frame=data.frame(redshift=redshift,magnitude=magnitude,bin=bin)





f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.1573, 0.5, 0.84, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


bins<-ggplot(frame,aes(factor(bin),magnitude)) +theme_bw()
bins<- bins + theme(axis.title.y=element_text(size=20,angle=90,hjust=0.5,vjust=1),axis.title.x=element_text(size=25,vjust=-2),axis.text=element_text(size=12),panel.border=element_rect(size=1.5),legend.position="none")
bins<- bins +stat_summary(fun.data = f, geom="errorbar",aes(fill=bin,width=0.5)) +stat_summary(fun.data = f, geom="boxplot",aes(fill=bin,width=0.5),position="dodge") + xlab("Redshift Bin") + ylab("Magnitude")
bins<- bins + scale_x_discrete(breaks=c("Bin1","Bin2","Bin3","Bin4"),labels=c("0.43 < z < 0.5","0.5 < z < 0.55","0.55 < z < 0.6","0.6 < z < 0.7") ) + theme(plot.margin=unit(c(0.1,0.1,1,0.5),"cm"))
bins<- bins + scale_y_reverse()
bins

