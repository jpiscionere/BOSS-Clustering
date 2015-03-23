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



m<-ggplot(data,aes(x=data$V3,y=data$V4)) +ylab(expression(paste(italic(M[g]))))+ xlab("Redshift") + ylim(-20,-24) + xlim(0.35,0.75) + geom_point(color="grey",alpha=0.4,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none")

m<- m + stat_density2d()


gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols=gg_color_hue(4)


bin1<-data.frame(y=c(-20.5,-20.5,-23,-23), x=c(0.43,0.5,0.5,0.43))
bin2<-data.frame(y=c(-20.5,-20.5,-23,-23), x=c(0.5,0.55,0.55,0.5))
bin3<-data.frame(y=c(-20.5,-20.5,-23,-23), x=c(0.55,0.6,0.6,0.55))
bin4<-data.frame(y=c(-20.5,-20.5,-23,-23), x=c(0.6,0.7,0.7,0.6))

m<- m + geom_polygon(data=bin1,mapping=aes(x=x,y=y),alpha=0.25,fill=cols[1])
m<- m + geom_polygon(data=bin2,mapping=aes(x=x,y=y),alpha=0.25,fill=cols[2])
m<- m + geom_polygon(data=bin3,mapping=aes(x=x,y=y),alpha=0.25,fill=cols[3])
m<- m + geom_polygon(data=bin4,mapping=aes(x=x,y=y),alpha=0.25,fill=cols[4])



dd<-with(density(data$V3),data.frame(x,y))
ggplot(data = dd, mapping = aes(x = x, y = y)) +
geom_line(color="black") + xlim(0.3,0.8) + 
layer(data = dd, mapping = aes(x=ifelse(x < 0.5 & x > 0.43,x,0), y=y), geom = "area", geom_params=list(fill=cols[1],alpha=.5)) + 
layer(data = dd, mapping = aes(x=ifelse(x < 0.55 & x > 0.5,x,0), y=y), geom = "area", geom_params=list(fill=cols[2],alpha=.5)) + 
layer(data = dd, mapping = aes(x=ifelse(x < 0.6 & x > 0.55,x,0), y=y), geom = "area", geom_params=list(fill=cols[3],alpha=.5)) + 
layer(data = dd, mapping = aes(x=ifelse(x < 0.7 & x > 0.6,x,0), y=y), geom = "area", geom_params=list(fill=cols[4],alpha=.5)) + 
scale_y_continuous(limits = c(0,max(dd$y)), name=" ",labels=NULL) + theme_bw() + xlab("Redshift")+ theme(axis.title.x=element_text(size=25,vjust=-0.2),axis.text.x=element_text(size=20),panel.border=element_rect(size=1.1),axis.ticks.y=element_blank(),axis.text.y=element_blank()) + geom_line(color="black")


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






