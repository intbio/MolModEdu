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
dna_avr<-read.csv('../analysis_data/dna2_param_df_md_avr.csv')
dna_avr=rename(dna_avr, c("x_av"='x',"y_av"='y',"z_av"='z', "Roll_av"="Roll",'Twist_av'='Twist','Slide_av'='Slide','P_1_av'='P_1','P_2_av'='P_2','chi_1_av'='chi_1','chi_2_av'='chi_2'))

dna_cryst<-read.csv('../analysis_data/dna2_param_df_cryst.csv')

d_c=melt(dna_cryst,id.vars=c('BPnum'),measure.vars=c('x','y','z','Roll','Twist','Slide','chi_1','chi_2'))
d_md=melt(dna_avr,id.vars=c('BPnum'),measure.vars=c('x','y','z','Roll','Twist','Slide','chi_1','chi_2'))
d_c$data='X-ray'
d_md$data='MD average'
d=rbind(d_c,d_md)


theme_set(theme_bw(base_size = 18))
###Roll and z-profiles
r<-ggplot(data=d[d$variable=='Roll'&d$BPnum<147,],aes(x=BPnum-74.5,y=value)) + ggtitle("Roll") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Roll'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)

z<-ggplot(data=d[d$variable=='z',],aes(x=BPnum-74,y=value)) + ggtitle("Z axis distance") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='z',]$value,na.rm=TRUE)*1.1,label=sequence),size=3)

q<-arrangeGrob(r,z,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_z.png",plot=q,width=15,height=10)
#adding RMSD
rd=r+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>2,],aes(x=BPnum-74.5,y=Roll,ymin=Roll-Roll_sd,ymax=Roll+Roll_sd))+
ggtitle("Roll and RMSD")
re=r+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>2,],aes(x=BPnum-74.5,y=Roll,ymin=Roll-Roll_er,ymax=Roll+Roll_er))+
ggtitle("Roll and error estimate for mean")

q<-arrangeGrob(rd,re,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_z_fluct.png",plot=q,width=15,height=10)

#We want pure geometry of DNA, z and x-y
dxy=dcast(d,BPnum+data~variable)
#make annotations
annot=dxy[dxy$data=='X-ray',]
annot$BPnum=annot$BPnum-74
annot=annot[annot$BPnum %in% seq(-80,80,5),]

xyp<-ggplot(data=dxy[dxy$BPnum<=74,],aes(x=x,y=y)) + ggtitle("XY geometry, Z>0") + 
xlab("X")+ylab('Y')+geom_path(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_y_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none")+
coord_fixed()+geom_text(data=annot[annot$BPnum<=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))

xym<-ggplot(data=dxy[dxy$BPnum>=74,],aes(x=x,y=y)) + ggtitle("XY geometry, Z<0") + 
xlab("X")+ylab('Y')+geom_path(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_y_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+
scale_x_continuous(breaks = round(seq(-50,50, by = 10),1),limits=c(-45,45))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none")+
coord_fixed()+geom_text(data=annot[annot$BPnum>=0,],aes(x=x*1.1,y=y*1.1,label=BPnum))

q<-arrangeGrob(xyp,xym,ncol=2)
z=z+theme(legend.position = "bottom")
q2<-arrangeGrob(q,z,ncol=1)

ggsave(filename="../analysis_data/dna2_geom_avr.png",plot=q2,width=15,height=13)



#Now let's do comparisson with ideal superhelix
#Assuming the porjection on X-Y plane should be a circle
#since we have already done the minimization when orienting the initial structure
#we can just the average distance from center.
rad=sqrt(mean(dxy[dxy$data=='X-ray','x']**2+dxy[dxy$data=='X-ray','y']**2))
slope=lm(z~BPnum,dxy)$coefficients[[2]]
idealxy=data.frame(x=rad*cos(seq(0,2*pi,0.01)),y=rad*sin(seq(0,2*pi,0.01)))
idealz=data.frame(BPnum=seq(-73,73,1),z=slope*seq(-73,73,1))

z=z+geom_line(data=idealz,aes(x=BPnum,y=z),color='green',size=2)
xyp=xyp+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2)
xym=xym+geom_path(data=idealxy,aes(x=x,y=y),color='green',size=2)


q<-arrangeGrob(xyp,xym,ncol=2)
# z=z+theme(legend.position = "bottom")
q2<-arrangeGrob(q,z,ncol=1)

ggsave(filename="../analysis_data/dna2_geom_avr_ideal.png",plot=q2,width=15,height=13)


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
ggsave(filename="../analysis_data/dna2_geom_avr_shifts.png",plot=q,width=17,height=9)


########################


quit()



