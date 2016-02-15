library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')
library(scales)

#data_bin12=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin12_wp_dr11v2_10bins.out")
data_bin1=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin1_wp_dr11v2_new")
data_bin1$V5="bin1"
data_bin2=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin2_wp_dr11v2_new")
data_bin2$V5="bin2"
data_bin2$V3=data_bin2$V3 + 9
data_bin3=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin3_wp_dr11v2_new")
data_bin3$V5="bin3"
data_bin3$V3=data_bin3$V3 + 7
data_bin4=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin4_wp_dr11v2_new")
data_bin4$V5="bin4"
data_bin4$V3=data_bin4$V3




df=data.frame(rbind(data_bin1,data_bin2,data_bin3,data_bin4))

df$ucl=df$V3+df$V4
df$lcl=df$V3-df$V4



ribbon<-ggplot(df,aes(x=V1,y=V3,color=factor(V5))) + scale_x_log10(limits=c(0.007,10),breaks=c(0.01,0.1,1.0,10) ) + scale_y_log10( ) +  annotation_logticks() + 
	xlab(expression(r[p]* (h^-1*Mpc))) +
	ylab(expression(w[p]*(r[p])*(h^-1*Mpc))) 
ribbon<-ribbon + geom_smooth(aes(ymin = lcl, ymax = ucl,fill=factor(V5)),lwd=1.2, data=df, stat="identity") + theme_bw() + 
	theme(legend.position=c(.75, .85),legend.key.size = unit(0.7, "cm"),legend.text=element_text(size=15),legend.title=element_text(size=15))
ribbon<-ribbon +
	scale_fill_discrete(guide = guide_legend(title = "Galaxy Sample")
	,labels=c("0.43 < z < 0,5","0.5 < z < 0.55","0.55 < z < 0.6","0.6 < z < 0.7")
	,breaks=c("bin1","bin2","bin3","bin4")) + 
	scale_color_discrete(guide = guide_legend(title = "Galaxy Sample")
        ,labels=c("0.43 < z < 0,5","0.5 < z < 0.55","0.55 < z < 0.6","0.6 < z < 0.7")
        ,breaks=c("bin1","bin2","bin3","bin4")) 

 

ribbon<- ribbon + geom_point(size=2) + theme(axis.title.y=element_text(size=16,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=15),panel.border=element_rect(size=1))

ribbon
