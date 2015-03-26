data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_all.out")

g_mod=data$V5
r_mod=data$V6
i_mod=data$V7
i_fib=data$V8

g_r=g_mod-r_mod
r_i=r_mod-i_mod

dperp=r_i - (g_r)/8.0

redshift=data$V3

bin=c(1:length(redshift))*0

bin[which(redshift > 0.43 & redshift < 0.5 )]="Bin1"
bin[which(redshift > 0.5 & redshift < 0.55 )]="Bin2"
bin[which(redshift > 0.55 & redshift < 0.6 )]="Bin2"
bin[which(redshift > 0.6 & redshift < 0.7)]="Bin4"


data=data.frame(g_r,r_i,i_fib,dperp,bin,redshift)

gr_ri<-ggplot(data,aes(x=r_i,y=g_r)) + ylab("r_model - i_model") + xlab("g_model - r_model") + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none") + ylim(0.2,1.5)

gr_ri<- gr_ri + stat_density2d()
gr_ri<- gr_ri + geom_abline(intercept=0.6,slope=1/8)


z_dperp<-ggplot(data,aes(x=redshift,y=dperp)) + ylab("d_perp") + xlab("Redshift") + geom_point(color="grey",alpha=0.1,shape=".") + theme_bw() + theme(axis.title.y=element_text(size=20,angle=90,face="italic"),axis.title.x=element_text(size=20,vjust=-0.2),axis.text=element_text(size=20),panel.border=element_rect(size=1.1),legend.position="none") + ylim(0.2,1.5)

z_dperp<- z_dperp + stat_density2d()
z_dperp<- z_dperp + geom_hline(aes(yintercept=0.55))
