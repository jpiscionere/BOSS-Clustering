data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_all.out",head=T)


g_r=data$modelflux_g-data$modelflux_r
r_i=data$modelflux_r-data$modelflux_i

dperp=r_i - (g_r)/8.0

redshift=data$V3

bin=c(1:length(data$redshift))*0

bin[which(data$redshift > 0.43 & data$redshift < 0.5 )]="Bin1"
bin[which(data$redshift > 0.5 & data$redshift < 0.55 )]="Bin2"
bin[which(data$redshift > 0.55 & data$redshift < 0.6 )]="Bin2"
bin[which(data$redshift > 0.6 & data$redshift < 0.7)]="Bin4"


frame=data.frame(g_r,r_i,i_fib,dperp,bin,redshift=data$redshift,Mag_i=data$Mag_i,i_cmod=data$i_cmod)

gr_ri<-ggplot(frame,aes(x=r_i,y=g_r)) + ylab(expression(r[mod] - i[mod])) + xlab(expression(g[mod] - r[mod])) + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none") + ylim(0.2,1.5)

gr_ri<- gr_ri + stat_density2d()
gr_ri<- gr_ri + geom_abline(intercept=0.6,slope=1/8)


z_dperp<-ggplot(frame,aes(x=redshift,y=dperp)) + ylab(expression(d[perp])) + xlab("Redshift") + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none") + ylim(0.2,1.5) + xlim(0.35,0.75)

z_dperp<- z_dperp + stat_density2d()
z_dperp<- z_dperp + geom_hline(aes(yintercept=0.55))


i_dperp<-ggplot(frame,aes(x=i_cmod,y=dperp)) + ylab(expression(d[perp])) + xlab(expression(i[cmod])) + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none") + ylim(0.2,1.5)

i_dperp<- i_dperp + stat_density2d()
i_dperp<- i_dperp + geom_hline(aes(yintercept=0.55))
i_dperp<- i_dperp + geom_abline(intercept=-11.575,slope=0.625)




z_mag<-ggplot(frame,aes(x=redshift,y=Mag_i)) +ylab(expression(paste(italic(M[i]))))+ xlab("Redshift") + ylim(-20,-24) + xlim(0.35,0.75) + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none")

z_mag<- z_mag + stat_density2d()

grid.arrange(gr_ri,z_dperp,i_dperp,z_mag)


