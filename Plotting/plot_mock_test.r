library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')

ylims=c(20,10000)

plot.new()
plot.window(xlab=NA,ylab=NA,ylim=ylims, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5,xlim=c(0.01,10),log="xy")

#data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin3_wp_dr11v2_full_map_jack100.out")
data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin_all_wp_dr11v2_new")
wp=data.frame(data)
plot(V2~V1,data=wp,log="xy",axes=FALSE,type="l",col="dodgerblue",xlab=NA,ylab=NA,ylim=ylims, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5,main="HighZ Sphere")
errbar(wp$V1,wp$V2,wp$V2+wp$V4,wp$V2-wp$V4,log="xy",add=TRUE,col="dodgerblue")


magaxis(log='xy')
box()
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p])*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=500)


data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test",pattern="wp.*.fibcol.out",full.names=T)
for( i in 1:length(data_files)){a=read.table(data_files[i]);lines(a$V1,a$V4,type="l",col=rgb(50,50,50,50,maxColorValue=255))}
#data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test/",pattern="norm_2.out",full.names=T)
#for( i in 1:length(data_files)){a=read.table(data_files[i]);lines(a$V1,a$V3,type="l",col=rgb(250,50,50,50,maxColorValue=255))}

data_files<-list.files("/home/piscioja/Clustering/Boss/Mock_Test",pattern="wp.*.nofibcol.out",full.names=T)
for( i in 1:length(data_files)){a=read.table(data_files[i]);lines(a$V1,a$V3,type="l",col=rgb(250,50,50,50,maxColorValue=255))} 

names=c("BOSS Data","Fiber Collisions","No Fiber Collisions")
colors=c("dodgerblue","darkgrey","red")
legend("topright",names,col=colors,pch=1,lty=1)
