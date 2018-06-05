#R-script to analyze average dna propterties
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
#DNA-protein interactions
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
# dna<-read.csv('../analysis_data/dna2_param_df_md.csv',na.strings=c("NA",'---'))
# nf=length(table(dna$Time))
# df=dna[dna$Time %in% seq(1,nf,100),c('Time','BPnum','x','y','z')]
# write.csv(df, file="../analysis_data/dna2_param_df_md_xyz.csv")

df<-read.csv('../analysis_data/dna2_param_df_md_xyz.csv',na.strings=c("NA",'---'))


# print(nf)
head(df)


dna_cryst<-read.csv('../analysis_data/dna2_param_df_cryst.csv')
dna_cryst$Time=1

d_c=melt(dna_cryst,id.vars=c('Time','BPnum'),measure.vars=c('x','y','z'))
d_md=melt(df,id.vars=c('Time','BPnum'),measure.vars=c('x','y','z'))
d_c$data='X-ray'
d_md$data='MD'
d=rbind(d_c,d_md)

head(d)


theme_set(theme_bw(base_size = 18))

#We want pure geometry of DNA, z and x-y
dxy=dcast(d,BPnum+data+Time~variable)



#make annotations
annot=dxy[dxy$data=='X-ray',]
annot$BPnum=annot$BPnum-74
annot=annot[annot$BPnum %in% seq(-80,80,5),]


#Now let's do comparisson with ideal superhelix
#Assuming the porjection on X-Y plane should be a circle
#since we have already done the minimization when orienting the initial structure
#we can just the average distance from center.
rad=sqrt(mean(dxy[dxy$data=='X-ray','x']**2+dxy[dxy$data=='X-ray','y']**2))
slope=lm(z~BPnum,dxy)$coefficients[[2]]
idealxy=data.frame(x=rad*cos(seq(0,2*pi,0.01)),y=rad*sin(seq(0,2*pi,0.01)))
idealz=data.frame(BPnum=seq(-73,73,1),z=slope*seq(-73,73,1))


xyp<-ggplot(data=dxy[dxy$BPnum<=74,],aes(x=x,y=y))+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + ggtitle("XY geometry, Z>0") + 
xlab("X, A")+ylab('Y, A')+#geom_path(aes(color=data),size=1)+
geom_point(aes(color=data),size=3)+geom_point(data=dxy[dxy$BPnum<=74 & dxy$data=="X-ray",],aes(x=x,y=y),color="red",size=3)+
scale_y_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none")+
coord_fixed()+geom_text(data=annot[annot$BPnum<=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))

xym<-ggplot(data=dxy[dxy$BPnum>=74,],aes(x=x,y=y))+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2) + ggtitle("XY geometry, Z<0") + 
xlab("X, A")+ylab('Y, A')+#geom_path(aes(color=data),size=1)+
geom_point(aes(color=data),size=3)+geom_point(data=dxy[dxy$BPnum>=74 & dxy$data=="X-ray",],aes(x=x,y=y),color="red",size=3)+
scale_y_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none")+
coord_fixed()+geom_text(data=annot[annot$BPnum>=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))

a=data.frame(BPnum=c(74),variable=c('z'),value=c(0),data=c('Ideal fit'),Time=c(1))
d=rbind(a,d)


z<-ggplot(data=d[d$variable=='z',],aes(x=BPnum-74,y=value)) +geom_line(data=idealz,aes(x=BPnum,y=z),color='green',size=2)+ ggtitle("Z axis distance") + ylab('Z, A') +
xlab("Base pair")+#geom_line(aes(color=data),size=1)+
geom_point(aes(color=data),size=3)+geom_point(data=d[d$variable=='z' & d$data=='X-ray',],aes(x=BPnum-74,y=value),color="red",size=3)+
scale_y_continuous(limits=c(-33,31),breaks = round(seq(-30,30, by = 10),1))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('MD'='blue','X-ray'='red','Ideal fit'='green'),name='',breaks=c('MD','X-ray','Ideal fit'),labels=c('MD snapshots','X-ray','Ideal superhelix fit'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=-33,label=sequence),size=3)

# z=z+scale_color_manual(data=idealz, values=c('blue','red','green'),name='',breaks=c('MD average','X-ray'),labels=c('MD average','X-ray'))
# z=z+
# xyp=xyp+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2)
# xym=xym+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2)


q<-arrangeGrob(xyp,xym,ncol=2)
z=z+theme(legend.position = "bottom")
q2<-arrangeGrob(q,z,ncol=1)

ggsave(filename="../analysis_data/pub_dna2_geom_fluct_ideal.png",plot=q2,width=15,height=13)

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





