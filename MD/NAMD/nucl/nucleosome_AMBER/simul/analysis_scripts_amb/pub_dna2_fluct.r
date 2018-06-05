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

dna_avr<-read.csv('../analysis_data/dna2_param_df_md_avr.csv')
# dna_avr=subset(dna_avr,BPnum<168&BPnum>20)


dna_avr=rename(dna_avr, c("x_av"='x',"y_av"='y',"z_av"='z', "Roll_av"="Roll",'Twist_av'='Twist','Slide_av'='Slide','P_1_av'='P_1','P_2_av'='P_2','chi_1_av'='chi_1','chi_2_av'='chi_2'))

dna_avr$r_sd=sqrt(dna_avr$x_sd**2+dna_avr$y_sd**2+dna_avr$z_sd**2)

# print(nf)
head(dna_avr)

dna_avr$BPnum=dna_avr$BPnum-74
# d=d[d$BPnum>0&d$BPnum<148,]
df=dna_avr

dfp=subset(df,BPnum>=0)
dfm=subset(df,BPnum<=0)
dfm$BPnum=(-1)*dfm$BPnum
dfp$turn='positive'
dfm$turn='negative'
dfd=rbind(dfp,dfm)


df_dna<-read.csv("../analysis_data/dna_rot_df_avr.csv",header=TRUE,check.name=FALSE)
# df_cryst<-read.csv("../analysis_data/dna_rot_df_cryst.csv",header=TRUE,check.name=FALSE)

# df_avr$DATA='Average'
# head(df_avr)
# df_cryst$DATA='X-ray'
# head(df_cryst)
# df=rbind(df_avr,df_cryst)
dfdd=df_dna

dfp=subset(dfdd,Basepair>=0)
dfm=subset(dfdd,Basepair<=0)
dfm$Basepair=(-1)*dfm$Basepair
dfp$turn='positive'
dfm$turn='negative'
dfdd=rbind(dfp,dfm)

dfdd[dfdd$Angle>0,'Angle']=dfdd[dfdd$Angle>0,'Angle']-180


theme_set(theme_bw(base_size = 18))


q1<-ggplot(data=dfd,aes(x=BPnum/10,y=r_sd,color=turn)) + ggtitle("MD simulations, NCPamb-model, RMSF of base pair center positions") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_line(size=2)+
# facet_grid(turn~.,scales='fixed')+
ylab('RMSF, Å')+
ylim(0,4)+
scale_x_continuous(limits=c(-0.05,9.2),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','','±1.0','','±2.0','','±3.0','','±4.0','','±5.0','','±6.0','','±7.0','','±8.0','','±9.0'),expand=c(0,0),name='Superhelix location (SHL)')+
# scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
# scale_color_manual(breaks=c(''),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
scale_color_manual(breaks=c('negative','positive'),labels=c('SHL < 0','SHL > 0'),values=c('negative'='indianred1','positive'='green4'),name='SHL')+
geom_line(data=dfdd,aes(x=Basepair/10,y=-Angle/180,alpha=turn),fill='black',color='black',size=1)+
geom_point(data=dfdd,aes(x=Basepair/10,y=-Angle/180,alpha=turn),fill='black',color='black',size=2)+
scale_alpha_manual(breaks=c('negative','positive'),labels=c('SHL<0','SHL>0'),values=c('negative'=0.3,'positive'=0.6),name='SHL')


q2<-ggplot(data=dfd,aes(x=BPnum/10,y=r_sd,color=turn)) + ggtitle("MD simulations, NCPch-model, RMSF of base pair center positions") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_line(size=2)+
# facet_grid(turn~.,scales='fixed')+
ylab('RMSF, Å')+
ylim(0,15)+
# scale_y_log10()+
scale_x_continuous(limits=c(0.05,9.2),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','','±1.0','','±2.0','','±3.0','','±4.0','','±5.0','','±6.0','','±7.0','','±8.0','','±9.0'),expand=c(0,0),name='Superhelix location (SHL)')+
# scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
# scale_color_manual(breaks=c(''),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
scale_color_manual(breaks=c('negative','positive'),labels=c('SHL < 0','SHL > 0'),values=c('negative'='indianred1','positive'='green4'),name='SHL')+
# geom_line(data=dfdd,aes(x=Basepair/10,y=-Angle/180,alpha=turn),fill='black',color='black',size=1)+
# geom_point(data=dfdd,aes(x=Basepair/10,y=-Angle/180,alpha=turn),fill='black',color='black',size=2)+
scale_alpha_manual(breaks=c('negative','positive'),labels=c('SHL<0','SHL>0'),values=c('negative'=0.3,'positive'=0.6),name='SHL')


# geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
# scale_y_continuous(limits=c(0,25),breaks = c(0,5,10,15),labels=c('0','5','10','15'),name='Stable contacts, total',expand=c(0,0))+

# geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=14,xmax=x+0.12,ymax=19,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)+
# geom_line(data=dfd,aes(x=Basepair/10,y=-Angle/18,alpha=turn),fill='black',color='black',size=1)+
# geom_point(data=dfd,aes(x=Basepair/10,y=-Angle/18,alpha=turn),fill='black',color='black',size=2)+
# scale_alpha_discrete(range = c(0.3, 0.6))
# scale_y_continuous(limits=c(0,25),breaks = c(0,5,10,15),labels=c('0','90','180','270'),name='Angle, deg',expand=c(0,0))

q<-arrangeGrob(q2,q1,ncol=1)

ggsave(filename="../analysis_data/pub_dna2_fluct_zoomed.png",plot=q1,width=11,height=5)
ggsave(filename="../analysis_data/pub_dna2_fluct_full.png",plot=q2,width=11,height=5)

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





