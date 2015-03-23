library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')

par(mfrow=c(2,2), mgp = c(1.5, 1, 0) ,mar=c(2.55, 3.25, 2.75, 1.5))
ylims=c(30,5000)

#hong=read.table("cmassdr9_hong.dat")
#hong_rp=hong$V1




########################Bin1################################


data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin1_wp_dr11v2_full_map_jack150_littleh.out")
wp=data.frame(data)
plot(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="deeppink",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p])*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=500)


magaxis(side=1:4,log='xy')
box()


errbar(wp$V1,wp$V2,wp$V2+wp$V4,wp$V2-wp$V4,log="xy",add=TRUE,col="deeppink")

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin1_wp_dr11v2_full_map_jack150_new_bins.out")
wp=data.frame(data)
points(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="grey",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/wp_bin1_imaging_test.out")
lines(data$V1,data$V3,log="xy",axes=FALSE,type="p",col="grey",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)


names=c("0.43 < z < 0.5")
colors=c("deeppink")
legend("topright",names,col=colors,lty=1,pch=1,bty="n",cex=1.25)


fibcol=6.082*62/1000
abline(v=fibcol,lty=5,col="grey")

deblend=0.0086
abline(v=deblend,lty=5,col="gold")


#points(hong$V1,hong$V2)
#errbar(hong$V1,hong$V2,hong$V2+hong$V3,hong$V2-hong$V3,log="xy",add=TRUE,col="black")



###########################MAIN19##############################



data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin2_wp_dr11v2_full_map_jack150_littleh.out")
wp=data.frame(data)
plot(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="dodgerblue",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p])*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=500)

magaxis(side=1:4,log='xy')
box()


errbar(wp$V1,wp$V2,wp$V2+wp$V4,wp$V2-wp$V4,log="xy",add=TRUE,col="dodgerblue")


data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin2_wp_dr11v2_full_map_jack150_new_bins.out")
wp=data.frame(data)
points(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="grey",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)



names=c("0.5 < z < 0.55")
colors=c("dodgerblue")
legend("topright",names,col=colors,lty=1,pch=1,bty="n",cex=1.25)

fibcol=6.344*62/1000
abline(v=fibcol,lty=5,col="grey")


deblend=0.009
abline(v=deblend,lty=5,col="gold")


#points(hong$V1,hong$V4)
#errbar(hong$V1,hong$V4,hong$V4+hong$V5,hong$V4-hong$V5,log="xy",add=TRUE,col="black")

##########################MAIN20###############################




data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin3_wp_dr11v2_full_map_jack150_littleh.out")
wp=data.frame(data)
plot(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="purple",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p])*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=500)

magaxis(side=1:4,log='xy')
box()


errbar(wp$V1,wp$V2,wp$V2+wp$V4,wp$V2-wp$V4,log="xy",add=TRUE,col="purple")

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin3_wp_dr11v2_full_map_jack150_new_bins.out")
wp=data.frame(data)
points(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="grey",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)




names=c("0.55 < z < 0.6")
colors=c("purple")
legend("topright",names,col=colors,lty=1,pch=1,bty="n",cex=1.25)

fibcol=6.673*62/1000
abline(v=fibcol,lty=5,col="grey")

deblend=0.0094
abline(v=deblend,lty=5,col="gold")


#points(hong$V1,hong$V6)
#errbar(hong$V1,hong$V6,hong$V6+hong$V7,hong$V6-hong$V7,log="xy",add=TRUE,col="black")


##########################MAIN21##################################


data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin4_wp_dr11v2_full_map_jack150_littleh.out")
wp=data.frame(data)
plot(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="seagreen2",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)
mtext(expression(r[p]* (h^-1*Mpc)),side=1,line=1.75,cex=0.81,at=0.45)
mtext(expression(w[p]*(r[p])*(h^-1*Mpc)),side=2,line=1.75,las=0,cex=0.81,at=500)

magaxis(side=1:4,log='xy')
box()


errbar(wp$V1,wp$V2,wp$V2+wp$V4,wp$V2-wp$V4,log="xy",add=TRUE,col="seagreen2")

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/bin4_wp_dr11v2_full_map_jack150_new_bins.out")
wp=data.frame(data)
points(V2~V1,data=wp,log="xy",axes=FALSE,type="p",col="grey",ylim=ylims,xlab=NA,ylab=NA, cex.lab=1, cex.axis=3, cex.main=1.5, cex.sub=1.5)



names=c("0.6 < z < 0.7")
colors=c("seagreen2")
legend("topright",names,col=colors,lty=1,pch=1,bty="n",cex=1.25)



fibcol=7.14*62/1000
abline(v=fibcol,lty=5,col="grey")

deblend=0.010
abline(v=deblend,lty=5,col="gold")



#points(hong$V1,hong$V8)
#errbar(hong$V1,hong$V8,hong$V8+hong$V9,hong$V8-hong$V9,log="xy",add=TRUE,col="black")


