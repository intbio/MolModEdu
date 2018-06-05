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
dna_avr=rename(dna_avr, c('Shear_av'='Shear','Stretch_av'='Stretch','Stagger_av'='Stagger','Buckle_av'='Buckle','Prop.Tw_av'='Prop.Tw','Opening_av'='Opening','Shift_av'='Shift','Tilt_av'='Tilt','Rise_av'='Rise','Pairing_av'='Pairing',"x_av"='x',"y_av"='y',"z_av"='z', "Roll_av"="Roll",'Twist_av'='Twist','Slide_av'='Slide','P_1_av'='P_1','P_2_av'='P_2','chi_1_av'='chi_1','chi_2_av'='chi_2'))

dna_cryst<-read.csv('../analysis_data/dna2_param_df_cryst.csv')

d_c=melt(dna_cryst,id.vars=c('BPnum'),measure.vars=c('Shear','Stretch','Stagger','Buckle','Prop.Tw','Opening','Shift','Tilt','Rise','x','y','z','Roll','Twist','Slide','chi_1','chi_2'))
d_md=melt(dna_avr,id.vars=c('BPnum'),measure.vars=c('Shear','Stretch','Stagger','Buckle','Prop.Tw','Opening','Shift','Tilt','Rise','x','y','z','Roll','Twist','Slide','chi_1','chi_2'))
d_c$data='X-ray'
d_md$data='MD average'
d=rbind(d_c,d_md)


theme_set(theme_bw(base_size = 18))

#Pairing before all
p<-ggplot(data=dna_avr,aes(x=BPnum-74,y=Pairing*100)) + ggtitle("Pairing") + 
xlab("Base pair")+ylab("Percentage")+geom_bar(stat='identity',size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=-1.1,label=sequence),size=3)

ggsave(filename="../analysis_data/dna2_param_pairing.png",plot=p,width=15,height=6)



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



t<-ggplot(data=d[d$variable=='Twist'&d$BPnum<147&d$BPnum>2,],aes(x=BPnum-74.5,y=value)) + ggtitle("Twist") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Twist'&d$BPnum<147&d$BPnum>3,]$value,na.rm=TRUE)*0.95,label=sequence),size=3)

s<-ggplot(data=d[d$variable=='Slide'&d$BPnum<147,],aes(x=BPnum-74.5,y=value)) + ggtitle("Slide") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Slide'&d$BPnum<147,]$value,na.rm=TRUE)*1.2,label=sequence),size=3)


q<-arrangeGrob(t,s,ncol=1)
ggsave(filename="../analysis_data/dna2_param_twist_slide.png",plot=q,width=15,height=10)

q<-arrangeGrob(r,t,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_twist.png",plot=q,width=15,height=10)


q<-arrangeGrob(r,s,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_slide.png",plot=q,width=15,height=10)


#######All 12 parameters

#we already have roll, twist, slide and their versions with RMSD
r=r
rd=rd+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+ggtitle('Roll')
t=t+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))
td=t+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74.5,y=Twist,ymin=Twist-Twist_sd,ymax=Twist+Twist_sd))
s=s+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))
sd=s+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74.5,y=Slide,ymin=Slide-Slide_sd,ymax=Slide+Slide_sd))


sh<-ggplot(data=d[d$variable=='Shift'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74.5,y=value)) + ggtitle("Shift") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Shift'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
shd=sh+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74.5,y=Shift,ymin=Shift-Shift_sd,ymax=Shift+Shift_sd))


tl<-ggplot(data=d[d$variable=='Tilt'&d$BPnum<147&d$BPnum>3,],aes(x=BPnum-74.5,y=value)) + ggtitle("Tilt") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Tilt'&d$BPnum<147&d$BPnum>3,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
tld=tl+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74.5,y=Tilt,ymin=Tilt-Tilt_sd,ymax=Tilt+Tilt_sd))


ri<-ggplot(data=d[d$variable=='Rise'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74.5,y=value)) + ggtitle("Rise") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Rise'&d$BPnum<147&d$BPnum>3,]$value,na.rm=TRUE)*0.9,label=sequence),size=3)
rid=ri+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74.5,y=Rise,ymin=Rise-Rise_sd,ymax=Rise+Rise_sd))



she<-ggplot(data=d[d$variable=='Shear'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74,y=value)) + ggtitle("Shear") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Shear'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
shed=she+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74,y=Shear,ymin=Shear-Shear_sd,ymax=Shear+Shear_sd))


st<-ggplot(data=d[d$variable=='Stretch'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74,y=value)) + ggtitle("Stretch") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Stretch'&d$BPnum<147&d$BPnum>3,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
std=st+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74,y=Stretch,ymin=Stretch-Stretch_sd,ymax=Stretch+Stretch_sd))


sta<-ggplot(data=d[d$variable=='Stagger'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74,y=value)) + ggtitle("Stagger") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Stagger'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
stad=sta+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74,y=Stagger,ymin=Stagger-Stagger_sd,ymax=Stagger+Stagger_sd))


b<-ggplot(data=d[d$variable=='Buckle'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74,y=value)) + ggtitle("Buckle") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Buckle'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
bd=b+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74,y=Buckle,ymin=Buckle-Buckle_sd,ymax=Buckle+Buckle_sd))


pt<-ggplot(data=d[d$variable=='Prop.Tw'&d$BPnum<147&d$BPnum>0,],aes(x=BPnum-74,y=value)) + ggtitle("Prop.Tw") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Prop.Tw'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
ptd=pt+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74,y=Prop.Tw,ymin=Prop.Tw-Prop.Tw_sd,ymax=Prop.Tw+Prop.Tw_sd))


o<-ggplot(data=d[d$variable=='Opening'&d$BPnum<147&d$BPnum>3,],aes(x=BPnum-74,y=value)) + ggtitle("Opening") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Opening'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
od=o+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>4,],aes(x=BPnum-74,y=Opening,ymin=Opening-Opening_sd,ymax=Opening+Opening_sd))


###And  all x y z 
z=z+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))
zd=z+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>2,],aes(x=BPnum-74,y=z,ymin=z-z_sd,ymax=z+z_sd))

x<-ggplot(data=d[d$variable=='x',],aes(x=BPnum-74,y=value)) + ggtitle("X axis distance") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='x',]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
xd=x+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>2,],aes(x=BPnum-74,y=x,ymin=x-x_sd,ymax=x+x_sd))


y<-ggplot(data=d[d$variable=='y',],aes(x=BPnum-74,y=value)) + ggtitle("Y axis distance") + 
xlab("Base pair")+geom_line(aes(color=data),size=1)+geom_point(aes(color=data),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),plot.margin = unit(c(-0.15,0.2,-0.15,0), "cm"))+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='y',]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
yd=y+geom_errorbar(data=dna_avr[dna_avr$BPnum<147&dna_avr$BPnum>2,],aes(x=BPnum-74,y=y,ymin=y-y_sd,ymax=y+y_sd))


##Plot all profiles
#Manual tweeks
shed=shed+ylim(-2,1.5)
od=od+ylim(-20,20)
td=td+ylim(10,45)

#all on one
q<-arrangeGrob(shed,shd,std,sd,stad,rid,bd,tld,ptd,rd,od,td,ncol=2)
ggsave(filename="../analysis_data/dna2_param_12.png",plot=q,width=25,height=15)


q<-arrangeGrob(shed,std,stad,ncol=1)
ggsave(filename="../analysis_data/dna2_param_bp_tr.png",plot=q,width=15,height=10)


q<-arrangeGrob(bd,ptd,od,ncol=1)
ggsave(filename="../analysis_data/dna2_param_bp_ang.png",plot=q,width=15,height=10)


q<-arrangeGrob(shd,sd,rid,ncol=1)
ggsave(filename="../analysis_data/dna2_param_ds_tr.png",plot=q,width=15,height=10)


q<-arrangeGrob(tld,rd,td,ncol=1)
ggsave(filename="../analysis_data/dna2_param_ds_ang.png",plot=q,width=15,height=10)


q<-arrangeGrob(xd,yd,zd,ncol=1)
ggsave(filename="../analysis_data/dna2_param_xyz.png",plot=q,width=15,height=10)


quit()


#################
#################
#Distribution of roll/twist/slide
rd<-ggplot(data=d[d$variable=='Roll'&d$BPnum<147,],aes(x=value,color=data)) + ggtitle("Roll distribution along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_freqpoly(binwidth=2,size=2)+scale_x_continuous(breaks = round(seq(-30,30, by = 5),1),limits=c(-30,30))

ggsave(filename="../analysis_data/dna2_param_roll_D.png",plot=rd,width=10,height=7)

td<-ggplot(data=d[d$variable=='Twist'&d$BPnum<147,],aes(x=value,color=data)) + ggtitle("Twist distribution along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_freqpoly(binwidth=2,size=2)+scale_x_continuous(breaks = round(seq(0,50, by = 5),1),limits=c(0,50))

ggsave(filename="../analysis_data/dna2_param_twist_D.png",plot=td,width=10,height=7)


sd<-ggplot(data=d[d$variable=='Slide'&d$BPnum<147,],aes(x=value,color=data)) + ggtitle("Slide distribution along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_freqpoly(binwidth=0.1,size=2)+scale_x_continuous(breaks = round(seq(-5,5, by = 1),1),limits=c(-3,3))

ggsave(filename="../analysis_data/dna2_param_slide_D.png",plot=sd,width=10,height=7)



################
################
#Distributions of glycosidic angle

gd<-ggplot(data=d[(d$variable %in% c('chi_1','chi_2')) & (d$BPnum<147),],aes(x=value,color=data)) + ggtitle("Glycosidic angle distribution along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_freqpoly(binwidth=2,size=2)#+scale_x_continuous(breaks = round(seq(-30,30, by = 5),1),limits=c(-30,30))

ggsave(filename="../analysis_data/dna2_param_chi_D.png",plot=gd,width=10,height=7)

#P-sugar and chi correlation

d1_md=dna_avr[,c('chi_1','P_1')]
d2_md=dna_avr[,c('chi_2','P_2')]
d1_md=rename(d1_md, c("chi_1"='chi','P_1'='P'))
d2_md=rename(d2_md, c("chi_2"='chi','P_2'='P'))
d_md=rbind(d1_md,d2_md)
d_md$data='MD average'

d1_c=dna_cryst[,c('chi_1','P_1')]
d2_c=dna_cryst[,c('chi_2','P_2')]
d1_c=rename(d1_c, c("chi_1"='chi','P_1'='P'))
d2_c=rename(d2_c, c("chi_2"='chi','P_2'='P'))
d_c=rbind(d1_c,d2_c)
d_c$data='X-ray'

d=rbind(d_md,d_c)

gpc<-ggplot(data=d,aes(x=chi,y=P,color=data)) + ggtitle("CHI-P correlation along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_point(size=2)#+scale_x_continuous(breaks = round(seq(-30,30, by = 5),1),limits=c(-30,30))

ggsave(filename="../analysis_data/dna2_param_pucker_chi_C.png",plot=gpc,width=10,height=7)

################
################
#Distributions of P-sugar
#add seq columns to X-ray
dna_cryst$SEQ_1=dna_avr$SEQ_1
dna_cryst$SEQ_2=dna_avr$SEQ_2
#melt
d1_c=melt(dna_cryst,id.vars=c('SEQ_1'),measure.vars=c('P_1'))
d2_c=melt(dna_cryst,id.vars=c('SEQ_2'),measure.vars=c('P_2'))
d1_c=rename(d1_c, c("SEQ_1"='SEQ'))
d2_c=rename(d2_c, c("SEQ_2"='SEQ'))
d1_c$data='X-ray'
d2_c$data='X-ray'

d1_md=melt(dna_avr,id.vars=c('SEQ_1'),measure.vars=c('P_1'))
d2_md=melt(dna_avr,id.vars=c('SEQ_2'),measure.vars=c('P_2'))
d1_md=rename(d1_md, c("SEQ_1"='SEQ'))
d2_md=rename(d2_md, c("SEQ_2"='SEQ'))
d1_md$data='MD average'
d2_md$data='MD average'

d=rbind(d1_c,d2_c,d1_md,d2_md)

b_to_t=data.frame(SEQ=c('A','T','G','C'),TYPE=c('R','Y','R','Y'))
d=merge(d,b_to_t)

pdr<-ggplot(data=d[d$TYPE=='R',],aes(x=value,color=data)) + ggtitle("Purine (R): Pucker distribution along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_freqpoly(binwidth=2,size=2)#+scale_x_continuous(breaks = round(seq(0,50, by = 5),1),limits=c(0,50))

pdy<-ggplot(data=d[d$TYPE=='Y',],aes(x=value,color=data)) + ggtitle("Pyrimidine (Y): Pucker distribution along sequence(!)") + 
scale_color_manual(values=c('blue','red','green'))+
geom_freqpoly(binwidth=2,size=2)#+scale_x_continuous(breaks = round(seq(0,50, by = 5),1),limits=c(0,50))

q<-arrangeGrob(pdr,pdy,ncol=2)

ggsave(filename="../analysis_data/dna2_param_pucker_D.png",plot=q,width=16,height=7)




quit()

# pl[[i]]<-ggplot(data=data,aes_string(x='X',y=i))+
# geom_point(size=3)+
# #geom_line(size=3)+
# geom_errorbar(ymin=data[[i]]-data$sd,ymax=data[[i]]+data$sd, color="blue",size=3)+
# geom_errorbar(ymin=data[[i]]-data$se,ymax=data[[i]]+data$se, color="black",size=3)+
# xlab("Base pair / base pair step")+geom_line(data=dfcryst_t,aes_string(x='X',y=i),color="red",size=3)+
# ggtitle(paste(i,"profile"))

# q<-ggplot(data=int[(int$DNA_chain %in% c('CHI','CHJ')) &(int$DNA_part %in% c('sugar','phosphate')),],aes(x=DNA_resid,y=min_dist,color=DNA_chain)) + ggtitle("X-ray: DNA-protein interactions, minimum distance to CZ of canonical ARG") + 
# xlab("Base pair")+geom_line(size=1)+geom_point(subset = .(min_dist > 0),size=3)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red','magenta'))+
# facet_grid(DNA_part~.,scales='free')+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())





