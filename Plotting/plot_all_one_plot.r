library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')

ylims=c(50,50000)

plot.new()
plot.window(xlab=NA,ylab=NA,ylim=ylims, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5,xlim=c(0.01,10),log="xy")





data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin1_nojack_selection.out")
wp=data.frame(data)
#bins2=(data$V1[c(1:17)]+data$V1[1+c(1:17)])/2

plot(data$V1,data$V3,log="xy",axes=FALSE,type="l",col="deeppink",xlab=NA,ylab=NA,ylim=ylims, cex.lab=1,lwd=2, cex.axis=3, cex.main=1.5, cex.sub=1.5,main="BOSS CMASS")
#errbar(data$V1,wp$V3,wp$V3+wp$V4,wp$V3-wp$V4,log="xy",add=TRUE,col="dodgerblue")

#data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin1_nojack.out")
#lines(data$V1,data$V3*data$V1,col="deeppink",lwd=2,lty=2)

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin2_nojack_selection.out")
lines(data$V1,data$V3,col="dodgerblue",lwd=2)

#data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin2_nojack.out")
#lines(data$V1,data$V3*data$V1,col="dodgerblue",lwd=2,lty=2)

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin3_nojack_selection.out")
lines(data$V1,data$V3,col="purple",lwd=2)

#data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin3_nojack.out")
#lines(data$V1,data$V3*data$V1,col="purple",lwd=2,lty=2)


data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin4_nojack_selection.out")
lines(data$V1,data$V3,col="seagreen2",lwd=2)

#data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin4_nojack.out")
#lines(data$V1,data$V3*data$V1,col="seagreen2",lwd=2,lty=2)


data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin_all_nojack.out")
lines(data$V1,data$V3,lwd=2)


data=read.table("/hd0/Research/Clustering/Boss/DR3/wp_dr3_mag_cut")
lines(data$V1,data$V3,lwd=2,col="darkgreen",lty=2)

data=read.table("/hd0/Research/Clustering/Boss/DR7/wp_dr7_mag_cut")
lines(data$V1,data$V3,lwd=2,col="red",lty=2)




magaxis(log='xy')
box()
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p]) *r[p]*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=2000)



names=c("0.43 < z < 0.5","0.5 < z < 0.55","0.55 < z < 0.6","0.6 < z < 0.7","All","DR3","DR7")
colors=c("deeppink","dodgerblue","purple","seagreen2","black","darkgreen","red")
legend("topright",names,col=colors,pch=1,lty=c(1,1,1,1,1,2,2))
