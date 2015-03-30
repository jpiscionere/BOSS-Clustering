library('mapproj')
library("ggplot2")
library('gridExtra')
library('permute')

short_data=read.table('/hd0/Research/Clustering/Boss/dr11/dr11v2/jackknife_regions_2',nrows=5)
classes=sapply(short_data,class)
data=read.table('/hd0/Research/Clustering/Boss/dr11/dr11v2/jackknife_regions_2',colClasses=classes)
rands=data.frame(RA=data$V1,Dec=data$V2,jackid=data$V4)
cols=rainbow(length(unique(rands$jackid)),s=0.6,v=0.9)
cols=sample(cols,replace=FALSE)

breaks=seq(0,359,by=20)
labels=seq(0,359,by=20)
labels[c(1:length(labels))]=""
labels[0]="0"
labels[which(breaks==90)]=90
labels[which(breaks==180)]=180
labels[which(breaks==270)]=270


skyplot2<-ggplot(rands,aes(x=RA,y=Dec,group=jackid)) + geom_point(aes(colour=factor(jackid))) + scale_colour_manual(values=cols) 
skyplot2<-skyplot2 + ylab("Declination") + xlab("Right Ascension")  + theme_bw() + ylim(-20,89.999) + xlim(0,360) + coord_map(projection="aitoff",orientation=c(89.999,180,0)) 
skyplot2<-skyplot2 + scale_x_continuous(breaks=(0:8)*45,limits=c(0,360), labels=c("0","","","","180","","","","360")) + scale_y_continuous(breaks=c(0,30,60),labels=c("0","30","60")) 
skyplot2<-skyplot2 + theme(legend.position="none",panel.grid.major = element_line(colour = "grey90", size = 0.5))



