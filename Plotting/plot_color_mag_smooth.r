data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_all.out",head=T)


g_r=data$modelflux_g-data$modelflux_r
r_i=data$modelflux_r-data$modelflux_i

dperp=r_i - (g_r)/8.0

redshift=data$V3

bin=c(1:length(data$redshift))*0

bin[which(data$redshift > 0.43 & data$redshift <= 0.5 )]="Bin1"
bin[which(data$redshift > 0.5 & data$redshift <= 0.55 )]="Bin2"
bin[which(data$redshift > 0.55 & data$redshift <= 0.6 )]="Bin3"
bin[which(data$redshift > 0.6 & data$redshift <= 0.7)]="Bin4"


frame=data.frame(g_r,r_i,i_fib,dperp,bin,redshift=data$redshift,Mag_i=data$Mag_i,i_cmod=data$i_cmod)
frame$bin[which(frame$bin==0)]=NA

i_dperp<-ggplot(subset(frame,!is.na(bin)),aes(x=i_cmod,y=dperp,group=bin))  + ylab(expression(d[perp])) + xlab(expression(i[cmod])) + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=15,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1)) + ylim(0.2,1.5)


i_dperp<- i_dperp + geom_hline(aes(yintercept=0.55))
i_dperp<- i_dperp + geom_abline(intercept=-11.575,slope=0.625)
i_dperp<- i_dperp + stat_smooth(aes(fill=bin,col=bin)) +scale_colour_hue(guide="none")  + scale_fill_discrete(name="Redshift Sample\n",breaks=c("Bin1","Bin2","Bin3","Bin4"),labels=c("0.5 > z > 0.43","0.55 > z > 0.5","0.6 > z > 0.55","0.7 > z > 0.6" ))
i_dperp<- i_dperp + theme(legend.position=c(.25, .85),legend.text=element_text(size=15),legend.title=element_text(size=15))


gr_ri<-ggplot(subset(frame,!is.na(bin)),aes(x=r_i,y=g_r)) + ylab(expression(r[mod] - i[mod])) + xlab(expression(g[mod] - r[mod])) + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none") + ylim(0.2,1.5)

gr_ri<- gr_ri + geom_abline(intercept=0.6,slope=1/8)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols=gg_color_hue(4)




gr_ri_bin1<- gr_ri + stat_density2d(data=subset(frame,bin=="Bin1"),colour=cols[1],show_guide=FALSE)
gr_ri_bin2<- gr_ri + stat_density2d(data=subset(frame,bin=="Bin2"),colour=cols[2],show_guide=FALSE)
gr_ri_bin3<- gr_ri + stat_density2d(data=subset(frame,bin=="Bin3"),colour=cols[3],show_guide=FALSE)
gr_ri_bin4<- gr_ri + stat_density2d(data=subset(frame,bin=="Bin4"),colour=cols[4],show_guide=FALSE)
 

grid.arrange(gr_ri_bin1,gr_ri_bin2,gr_ri_bin3,gr_ri_bin4)

grid.arrange(gr_ri + theme(plot.margin=unit(c(0,1,0.5,2),"cm")),i_dperp + theme(plot.margin=unit(c(0,-1,0.5,0.5),"cm")),nrow=2,ncol=2)

