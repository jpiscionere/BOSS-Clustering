library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')
library(scales)


data_bin12=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin12_Nden_Matched_With_Errors.out")
data_bin12$V5="bin12"

data_dr7=read.table("/hd0/Research/Clustering/Boss/DR7/dr7_wp_full.out")
data_dr7$V5="LRG"

mean_wp=0
data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test",pattern="wp.*.fibcol_nearest_neighbor",full.names=T)
for( i in 1:length(data_files)){a=read.table(data_files[i]);mean_wp=mean_wp + a$V3}
mean_wp=mean_wp/length(data_files)
std_err=c(1:length(mean_wp))*0
for( i in 1:length(data_files)){a=read.table(data_files[i]);std_err = std_err + sqrt((a$V3-mean_wp)^2)}
std_err=std_err/sqrt(length(data_files))

data_model=data.frame(V1=a$V1,V2=mean_wp,V3=mean_wp,V4=std_err,V5="Model")


df=data.frame(rbind(data_bin12,data_dr7,data_model))

df$ucl=df$V3+df$V4
df$lcl=df$V3-df$V4

cols=c("#619CFF","#F8766D","black")


ribbon<-ggplot(df,aes(x=V1,y=V3,color=factor(V5))) + scale_x_log10(limits=c(0.007,10) ) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) + annotation_logticks() +
        xlab(expression(r[p]* (h^-1*Mpc))) +
        ylab(expression(w[p]*(r[p])*(h^-1*Mpc)))
ribbon<-ribbon + geom_smooth(aes(ymin = lcl, ymax = ucl,fill=factor(V5)),lwd=1.2, data=df, stat="identity") + theme_bw() +
        theme(legend.position=c(.75, .85),legend.key.size = unit(0.7, "cm"),legend.text=element_text(size=15),legend.title=element_text(size=15))
ribbon<-ribbon +
        scale_fill_manual(guide = guide_legend(title = "Galaxy Sample")
        ,labels=c("LRG: 0.16 < z < 0.36","CMASS: 0.43 < z < 0.55","Model Assuming NFW")
        ,breaks=c("LRG","bin12","Model")
        ,values=cols) +
        scale_color_manual(guide = guide_legend(title = "Galaxy Sample")
        ,labels=c("LRG: 0.16 < z < 0.36","CMASS: 0.43 < z < 0.55","Model Assuming NFW")
        ,breaks=c("LRG","bin12","Model")
        ,values=cols)

ribbon<- ribbon + theme(axis.title.y=element_text(size=16,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=15),panel.border=element_rect(size=1))

