data=read.table("/hd0/Research/Clustering/Boss/DR7/cross_match.out",head=T)
frame=data.frame(petro=data$petroMag_G,diff=(data$petroMag_G-data$Mag_G)/data$Mag_G)
c<-ggplot(frame,aes(petro,diff))
c + stat_smooth() + geom_point(alpha=0.2,aes(colour=diff)) + ylim(-0.05,0.025) + xlim(-23.2,-20.8)+ scale_colour_gradient(low = "blue",high="red",name="Difference") + xlab("Calculated Mag_G") + ylab("Fractional Difference") + labs(title="Comparing Given Mag G to Calculated M_G DR7") + theme(legend.position=c(0,-21)) + geom_vline(xintercept=-21.2) + theme_bw()
