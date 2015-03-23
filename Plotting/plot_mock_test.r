library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')

ylims=c(50,50000)
data=read.table("../DR7/masjedi_lrg.out")
mag_correction=data$V6
mag_correction[16]=1
mag_correction[17]=1

plot.new()
plot.window(xlab=NA,ylab=NA,ylim=ylims, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5,xlim=c(0.01,10),log="xy")

#data_files<-list.files(".",pattern="wp_dr3_mag_cut_",full.names=T)

#for(i in 1:length(data_files)){a=read.table(data_files[i]);lines(bins2,a$V3/mag_correction,type="l",col=rgb(150,150,150,75,maxColorValue=255))}
#data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin3_wp_dr11v2_full_map_jack100.out")
data=read.table("wp_dr3_mag_cut")
wp=data.frame(data)
bins2=(data$V1[c(1:17)]+data$V1[1+c(1:17)])/2

plot(bins2,data$V3/mag_correction,log="xy",axes=FALSE,type="l",col="dodgerblue",xlab=NA,ylab=NA,ylim=ylims, cex.lab=1,lwd=3, cex.axis=3, cex.main=1.5, cex.sub=1.5,main="DR3 LRG Test")
#errbar(wp$V1,wp$V2,wp$V2+wp$V4,wp$V2-wp$V4,log="xy",add=TRUE,col="dodgerblue")

bins2=(data$V1[c(1:17)]+data$V1[1+c(1:17)])/2
#lines(bins2,data$V3/mag_correction,lty=4)

data=read.table("../DR7/wp_dr7_mag_cut")
lines(bins2,data$V3,col="green3",lwd=3)


data=read.table("../DR7/masjedi_lrg.out")
lines(bins2[c(2:16)],data$V2,col="red",lwd=2)
errbar(bins2[c(2:16)],data$V2,data$V2+data$V3,data$V2-data$V3,log="xy",add=TRUE,col="red")



data=read.table("zehavi_wp.out")
lines(data$V1,data$V3,col="magenta",lwd=2,lty=2)


#data=read.table("wp_dr3_mag_cut")
#lines(data$V1,data$V3,col="forestgreen")








magaxis(log='xy')
box()
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p]) *r[p]*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=2000)


#data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test/Redshift_Test/LittleH",pattern="sphere.0.test.out",full.names=T)
#for( i in 1:length(data_files)){a=read.table(data_files[i]);lines(a$V1,a$V4,type="l",col=rgb(50,50,50,50,maxColorValue=255))}
#data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test/",pattern="norm_2.out",full.names=T)
#for( i in 1:length(data_files)){a=read.table(data_files[i]);lines(a$V1,a$V3,type="l",col=rgb(250,50,50,50,maxColorValue=255))}

#data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test/Redshift_Test/LittleH",pattern="wp_nojack_mstar_norm2_full",full.names=T)
#data_files=data_files[c(1:length(data_files)-1)]
#mean=c(1:length(data_files))*0
#for( i in 1:length(data_files)){a=read.table(data_files[i]);mean=mean+a$V3;lines(a$V1,a$V3,type="l",col=rgb(250,50,50,75,maxColorValue=255))} 
#lines(a$V1,mean/length(data_files),col="red3",lwd=3)



#data=read.table("hong_dr10.out")
#lines(data$V1,data$V2,col="green2",lwd=2,lty=2)
#lines(data$V1,data$V4,col="green3",lwd=2,lty=2)
#lines(data$V1,data$V6,col="green4",lwd=2,lty=2)

names=c("Me-DR3 Mag Cut","DR7","linear bin centers","Masjedi- Sample14 LRG")
colors=c("dodgerblue","green3","black","red")
legend("topright",names,col=colors,pch=c(1,1,-1,1),lty=c(1,1,0,1))
