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
# dna<-read.csv('../analysis_data/dna2_param_df_md.csv',na.strings=c("NA",'---'))
# nf=length(table(dna$Time))

# ???df=dna[dna$Time %in% seq(250,nf,10),c('Time','BPnum','x','y','z')]
# write.csv(df, file="../analysis_data/dna2_param_df_md_xyz.csv")

dna_avr<-read.csv('../analysis_data/dna2_param_df_md_avr.csv')
# dna_avr=subset(dna_avr,BPnum<168&BPnum>20)


dna_avr=rename(dna_avr, c("x_av"='x',"y_av"='y',"z_av"='z', "Roll_av"="Roll",'Twist_av'='Twist','Slide_av'='Slide','P_1_av'='P_1','P_2_av'='P_2','chi_1_av'='chi_1','chi_2_av'='chi_2'))
dna_avr$Time=1


df<-read.csv('../analysis_data/dna2_param_df_md_xyz.csv',na.strings=c("NA",'---'))


# print(nf)
head(df)


dna_cryst<-read.csv('../analysis_data/dna2_param_df_cryst.csv')
dna_cryst$Time=1

d_avr=melt(dna_avr,id.vars=c('Time','BPnum'),measure.vars=c('x','y','z'))

d_c=melt(dna_cryst,id.vars=c('Time','BPnum'),measure.vars=c('x','y','z'))
d_md=melt(df,id.vars=c('Time','BPnum'),measure.vars=c('x','y','z'))
d_c$data='Initial'
d_md$data='MD'
d_avr$data='Average'

d=rbind(d_c,d_md,d_avr)

d$BPnum=d$BPnum-20
# d=d[d$BPnum>0&d$BPnum<148,]

head(d)


theme_set(theme_bw(base_size = 18))

#We want pure geometry of DNA, z and x-y
dxy=dcast(d,BPnum+data+Time~variable)

head(dxy)

dxy=arrange(dxy,Time,BPnum)

head(dxy)


#make annotations
annot=dxy[dxy$data=='Initial',]
annot$BPnum=annot$BPnum-74
annot=annot[annot$BPnum %in% seq(-90,90,5),]

print(max(dxy$Time))


xyp<-ggplot(data=dxy[dxy$BPnum<=165 & dxy$BPnum>-18 & dxy$data=="MD",],aes(x=x,y=y,group=Time,color=Time))+
#geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + 
# ggtitle("XY geometry") + 
xlab("X, Å")+ylab('Y, Å')+geom_path(size=1)+
scale_alpha_continuous(range=c(1,1))+
# geom_point(aes(color=data),size=3)+
scale_color_gradient(low = "green", high = "blue", guide = "colourbar")+
geom_path(data=dxy[dxy$BPnum<=165 & dxy$BPnum>-18 & dxy$data=="Initial",],aes(x=x,y=y),color="red",size=2)+
# geom_path(data=dxy[dxy$BPnum<=165 & dxy$BPnum>-18 & dxy$data=="Average",],aes(x=x,y=y),color="blue",size=2)+

scale_y_continuous(breaks = round(seq(-50,110, by = 10),1),limits=c(-50,115))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-60,60))+
# scale_color_manual(values=c('MD'='blue','Initial'='red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none")+
# coord_fixed()+geom_text(data=annot[annot$BPnum<=0,],aes(x=x*1.1,y=y*1.1,label=BPnum),color='black')
geom_point(data=annot,aes(x=x,y=y),size=4,color='red')

xym<-ggplot(data=dxy[dxy$BPnum<=165 & dxy$BPnum>-18 & dxy$data=="MD",],aes(x=z,y=y,group=Time))+
#geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + 
# ggtitle("ZY geometry") + 
xlab("Z, Å")+ylab('Y, Å')+geom_path(aes(color=Time),size=1)+
scale_color_gradient(low = "green", high = "blue", guide = "colourbar",limits=c(0,1000),name='Time, ns')+
# scale_alpha_continuous(range=c(0.1,1))+
# geom_point(aes(color=data),size=3)+
geom_path(data=dxy[dxy$BPnum<=165 & dxy$BPnum>-18 & dxy$data=="Initial",],aes(x=z,y=y),color="red",size=2)+
scale_y_continuous(breaks = round(seq(-50,110, by = 10),1),limits=c(-50,115))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-50,50))+
# scale_color_manual(values=c('MD'='blue','Initial'='red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "right")+
# coord_fixed()+geom_text(data=annot[annot$BPnum<=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))
geom_point(data=annot,aes(x=z,y=y),size=4,color='red')

q<-arrangeGrob(xyp,xym,ncol=2)


ggsave(filename="../analysis_data/pub_dna2_geom_fluct_linkers.png",plot=q,width=12,height=8)

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





