#R-script to analyze average dna propterties
#show the fluctuations of DNA
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(gridExtra)
###################
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))
##############
###############
##Loading data frames
df<-read.csv('../analysis_data/dna2_param_df_md.csv',na.strings=c("NA",'---'))
# nf=length(table(dna$Time))

# ???df=dna[dna$Time %in% seq(250,nf,10),c('Time','BPnum','x','y','z')]
# write.csv(df, file="../analysis_data/dna2_param_df_md_xyz.csv")

# df<-read.csv('../analysis_data/dna2_param_df_md_xyz.csv',na.strings=c("NA",'---'))


# print(nf)
head(df)


dna_cryst<-read.csv('../analysis_data/dna2_param_df_cryst.csv')
dna_cryst=dna_cryst[,c('BPnum','x','y','z')]
dna_cryst$Time=1
df=df[,c('Time','BPnum','x','y','z')]
df$data='MD'
dna_cryst$data='X-ray'
d=rbind(df,dna_cryst)

d$BPnum=d$BPnum-94

d1=subset(d,BPnum>0)
d2=subset(d,BPnum<0)
d2$BPnum=d2$BPnum*(-1)

d=merge(d1,d2,by=c('Time','data','BPnum'))

d$dist=sqrt((d$x.x-d$x.y)**2+(d$y.x-d$y.y)**2+(d$z.x-d$z.y)**2)
d$distJ=sqrt((d$x.x+d$x.x)**2+(d$z.x+d$z.x)**2)
d$distI=sqrt((d$x.y+d$x.y)**2+(d$z.y+d$z.y)**2)


dav=ddply(d,c('data','BPnum'),summarize, dist_av=mean(dist))

theme_set(theme_bw(base_size = 18))

#We want pure geometry of DNA, z and x-y

q<-ggplot(data=subset(d,(data=='MD') & (BPnum==78)),aes(x=Time,y=dist))+
#geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + 
# ggtitle("XY geometry") + 
xlab("Time, ns")+ylab('Distance, Å')+
# geom_point(aes(color=data),size=3)+
# scale_color_gradient(low = "green", high = "blue", guide = "colourbar")+
# scale_fill_gradient2(limits=c(0,100),low="blue",mid="green", high="red",midpoint=50,guide_legend(title="Distance, A"))+
# scale_y_continuous(breaks = round(seq(6,10, by = 0.5),1),limits=c(7,10))+
# scale_x_continuous(breaks = round(seq(0,1000, by = 100),0),limits=c(0,1000))+
geom_line()+
geom_line(aes(x=Time,y=distI),color='red')+
geom_line(aes(x=Time,y=distJ),color='green')

# scale_color_manual(values=c('MD'='blue','Initial'='red','green'))+
# theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none")+
# coord_fixed()+geom_text(data=annot[annot$BPnum<=0,],aes(x=x*1.1,y=y*1.1,label=BPnum),color='black')
# geom_point(data=annot,aes(x=x,y=y),size=4,color='red')


ggsave(filename="../analysis_data/pub_dna2_fret_5bp.png",plot=q,width=12,height=8)

summary(subset(d,(data=='MD') & (BPnum==78) & (Time>250)))
summary(subset(d,(data=='MD') & (BPnum==83) & (Time>250)))
sd(subset(d,(data=='MD') & (BPnum==78) & (Time>250))$dist)

q()




q<-ggplot(data=subset(dav,BPnum==78),aes(x=BPnum/10,y=dist_av,color=data))+
#geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + 
# ggtitle("XY geometry") + 
xlab("SHL")+ylab('Distance')+geom_line()

ggsave(filename="../analysis_data/pub_dna2_fret_5bp.png",plot=q,width=12,height=8)


dcor=ddply(subset(d,data=='MD'),c('BPnum'),summarize, cor=(cov(x.x,x.y)+cov(y.x,y.y)+cov(z.x,z.y))/(sqrt((var(x.x)+var(y.x)+var(z.x))*((var(x.y)+var(y.y)+var(z.y))))))


q<-ggplot(data=subset(dcor,BPnum>70),aes(x=BPnum/10,y=cor))+
#geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + 
# ggtitle("XY geometry") + 
xlab("SHL")+ylab('Correlation')+geom_line()

ggsave(filename="../analysis_data/pub_dna2_geom_linker_cov.png",plot=q,width=12,height=8)



q()

z<-ggplot(data=d[d$variable=='z',],aes(x=BPnum-74,y=value)) + ggtitle("Z axis distance") + ylab('Z, A') +
xlab("Base pair")+#geom_line(aes(color=data),size=1)+
geom_point(aes(color=data),size=3)+geom_point(data=d[d$variable=='z' & d$data=='Initial',],aes(x=BPnum-74,y=value),color="red",size=3)+
scale_y_continuous(limits=c(-33,31),breaks = round(seq(-30,30, by = 10),1))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('MD'='blue','Initial'='red'),name='',breaks=c('MD','Initial'),labels=c('MD snapshots','Initial'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=-33,label=sequence),size=3)

# z=z+scale_color_manual(data=idealz, values=c('blue','red','green'),name='',breaks=c('MD average','X-ray'),labels=c('MD average','X-ray'))
# z=z+
# xyp=xyp+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2)
# xym=xym+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2)


q<-arrangeGrob(xyp,xym,ncol=2)
z=z+theme(legend.position = "bottom")
q2<-arrangeGrob(q,z,ncol=1)

ggsave(filename="../analysis_data/pub_dna2_geom_fluct_linkers.png",plot=q2,width=15,height=13)

quit()

#Let's try to tack shifts
###########

xyp<-ggplot(data=dxy[dxy$BPnum<=74,],aes(x=x,y=y)) + ggtitle("XY geometry, Z>0") + 
xlab("X")+ylab('Y')+geom_path(aes(alpha=data),size=1)+geom_point(aes(shape=data,color=factor((BPnum-74)%%4)),size=3)+
scale_y_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+scale_color_manual(values=c('blue','red','green','yellow'),name='track')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "bottom")+
coord_fixed()+geom_text(data=annot[annot$BPnum<=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))+scale_alpha_manual(values=c(1,0.5),name='data')

xym<-ggplot(data=dxy[dxy$BPnum>=74,],aes(x=x,y=y)) + ggtitle("XY geometry, Z<0") + 
xlab("X")+ylab('Y')+geom_path(aes(alpha=data),size=1)+geom_point(aes(shape=data,color=factor((BPnum-74)%%4)),size=3)+
scale_y_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+scale_color_manual(values=c('blue','red','green','yellow'),name='track')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "bottom")+
coord_fixed()+geom_text(data=annot[annot$BPnum>=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))+scale_alpha_manual(values=c(1,0.5),name='data')
q<-arrangeGrob(xyp,xym,ncol=2)
# ggsave(filename="../analysis_data/dna2_geom_avr_shifts.png",plot=q,width=17,height=9)


########################





