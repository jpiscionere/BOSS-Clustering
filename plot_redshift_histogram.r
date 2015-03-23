library("MASS")
library("ggplot2")
library("gridExtra")
library("Hmisc")
library('magicaxis')
library("foreign")
library("gplots")
library('gridExtra')

data=read.table("/hd0/Research/Clustering/Boss/dr11/dr11v2/Redshift_Magnitude.out")




gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols=gg_color_hue(4)

dd<-with(density(data$V3),data.frame(x,y))
ggplot(data = dd, mapping = aes(x = x, y = y)) +
geom_line(color="black") + xlim(0.3,0.8) +
layer(data = dd, mapping = aes(x=ifelse(x < 0.5 & x > 0.43,x,0), y=y), geom = "area", geom_params=list(fill=cols[1],alpha=.5)) +
layer(data = dd, mapping = aes(x=ifelse(x < 0.55 & x > 0.5,x,0), y=y), geom = "area", geom_params=list(fill=cols[2],alpha=.5)) +
layer(data = dd, mapping = aes(x=ifelse(x < 0.6 & x > 0.55,x,0), y=y), geom = "area", geom_params=list(fill=cols[3],alpha=.5)) +
layer(data = dd, mapping = aes(x=ifelse(x < 0.7 & x > 0.6,x,0), y=y), geom = "area", geom_params=list(fill=cols[4],alpha=.5)) +
scale_y_continuous(limits = c(0,max(dd$y)), name=" ",labels=NULL) + theme_bw() + xlab("Redshift")+ theme(axis.title.x=element_text(size=25,vjust=-0.2),axis.text.x=element_text(size=20),panel.border=element_rect(size=1.1),axis.ticks.y=element_blank(),axis.text.y=element_blank()) + geom_line(color="black")
















