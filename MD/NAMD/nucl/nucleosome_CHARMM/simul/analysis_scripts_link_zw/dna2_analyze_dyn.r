#R-script to analyze 
#dynamical properties of DNA
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
dna_act<-read.csv('../analysis_data/dna2_param_df_md_act.csv')
dna_act=subset(dna_act,BPnum<168&BPnum>20)
d_md=melt(dna_act,id.vars=c('BPnum'),measure.vars=c('x','y','z','Roll','Twist','Slide','chi_1','chi_2'))
d=d_md

theme_set(theme_bw(base_size = 18))
###Roll and z-profiles ACT
r<-ggplot(data=d[d$variable=='Roll'&d$BPnum<187,],aes(x=BPnum-94.5,y=value)) + ggtitle("Roll") + 
xlab("Base pair")+ylab('Auto-correlation time, ns')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Roll'&d$BPnum<187,]$value,na.rm=TRUE)*0.1,label=sequence),size=3)

z<-ggplot(data=d[d$variable=='z',],aes(x=BPnum-94,y=value)) + ggtitle("Z axis distance") + 
xlab("Base pair")+ylab('Auto-correlation time, ns')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='z',]$value,na.rm=TRUE)*0.1,label=sequence),size=3)

q<-arrangeGrob(r,z,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_z_act.png",plot=q,width=15,height=10)


t<-ggplot(data=d[d$variable=='Twist'&d$BPnum<187,],aes(x=BPnum-94.5,y=value)) + ggtitle("Twist") + 
xlab("Base pair")+ylab('Auto-correlation time, ns')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=-1.3,label=sequence),size=3)


s<-ggplot(data=d[d$variable=='Slide'&d$BPnum<187,],aes(x=BPnum-94.5,y=value)) + ggtitle("Slide") + 
xlab("Base pair")+ylab('Auto-correlation time, ns')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=-1.3,label=sequence),size=3)


q<-arrangeGrob(t,s,ncol=1)
ggsave(filename="../analysis_data/dna2_param_twist_slide_act.png",plot=q,width=15,height=10)

q<-arrangeGrob(r,t,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_twist_act.png",plot=q,width=15,height=10)


q<-arrangeGrob(r,s,ncol=1)
ggsave(filename="../analysis_data/dna2_param_roll_slide_act.png",plot=q,width=15,height=10)


#Let's go for cross-correlation analysis
dna_cor<-read.csv('../analysis_data/dna2_param_df_md_cor.csv')
dna_cor=subset(dna_cor,BPnum1<168&BPnum1>20&BPnum2<168&BPnum2>20)
dna_cor=subset(dna_cor,BPnum1!=24+20+74&BPnum1!=25+20+74&BPnum2!=24+20+74&BPnum2!=25+20+74)


d=dna_cor


#Now we are ready to calculate correlations	
# <dTW(i), dTW(j)> 
# <dRL(i), dRL(j)> 
# <dSL(i), dSL(j)> 
# For all pairs (i,j), i.e., matrices 146 x 146. 
# [here d is for delta; I hope we understand each other]
# If possible, I would like to have access to numeric data.
# I expect that TW(i) and TW(i+1) should be anti-correlated, but I don’t know what to expect for TW(i) and TW(i+5)… the same for ROLL & SLIDE.

# Also, please calculate
# <dTW(i), dRL(i)> 
# <dTW(i), dSL(i)> 
# <dSL(i), dRL(i)> 

#Let's start with simple correlations

tr<-ggplot(data=d[d$variable1=='Twist'&d$variable2=='Roll'&d$BPnum1==d$BPnum2,],aes(x=BPnum1-94.5,y=corr)) + ggtitle("Twist-roll correlation") + 
xlab("Base pair")+ylab('Correlation')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable1=='Twist'&d$variable2=='Roll'&d$BPnum1==d$BPnum2,]$corr,na.rm=TRUE)*1.1,label=sequence),size=3)

ggsave(filename="../analysis_data/dna2_param_twist_roll_cor.png",plot=tr,width=15,height=7)


ts<-ggplot(data=d[d$variable1=='Twist'&d$variable2=='Slide'&d$BPnum1==d$BPnum2,],aes(x=BPnum1-94.5,y=corr)) + ggtitle("Twist-slide correlation") + 
xlab("Base pair")+ylab('Correlation')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable1=='Twist'&d$variable2=='Slide'&d$BPnum1==d$BPnum2,]$corr,na.rm=TRUE)*1.1,label=sequence),size=3)


sr<-ggplot(data=d[d$variable1=='Slide'&d$variable2=='Roll'&d$BPnum1==d$BPnum2,],aes(x=BPnum1-94.5,y=corr)) + ggtitle("Slide-roll correlation") + 
xlab("Base pair")+ylab('Correlation')+geom_line(size=1)+geom_point(size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
geom_text(data=seqdf,aes(x=X,y=min(d[d$variable1=='Slide'&d$variable2=='Roll'&d$BPnum1==d$BPnum2,]$corr,na.rm=TRUE)*1.1,label=sequence),size=3)

q<-arrangeGrob(tr,ts,sr,ncol=1)

ggsave(filename="../analysis_data/dna2_param_cor.png",plot=q,width=15,height=15)


#Let's go for correlation tile plots

ttm<-ggplot(data=d[d$variable1=='Twist'&d$variable2=='Twist',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("Twist-twist cross correlation") + 
xlab("Base pair (twist)")+ylab('Base pair (twist)')+geom_tile(aes(fill=corr))+scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())

ggsave(filename="../analysis_data/dna2_param_twist_twist_cor_mat.png",plot=ttm,width=12,height=10)


zzm<-ggplot(data=d[d$variable1=='z'&d$variable2=='z',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("Z-z cross correlation") + 
xlab("Base pair (z)")+ylab('Base pair (z)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_z_z_cor_mat.png",plot=zzm,width=12,height=10)


rrm<-ggplot(data=d[d$variable1=='Roll'&d$variable2=='Roll',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("Roll-roll cross correlation") + 
xlab("Base pair (Roll)")+ylab('Base pair (Roll)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_roll_roll_cor_mat.png",plot=rrm,width=12,height=10)


ssm<-ggplot(data=d[d$variable1=='Slide'&d$variable2=='Slide',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("Slide-slide cross correlation") + 
xlab("Base pair (Slide)")+ylab('Base pair (Slide)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_slide_slide_cor_mat.png",plot=ssm,width=12,height=10)


xxm<-ggplot(data=d[d$variable1=='x'&d$variable2=='x',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("X-x cross correlation") + 
xlab("Base pair (x)")+ylab('Base pair (x)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_x_x_cor_mat.png",plot=xxm,width=12,height=10)


yym<-ggplot(data=d[d$variable1=='y'&d$variable2=='y',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("y-y cross correlation") + 
xlab("Base pair (y)")+ylab('Base pair (y)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_y_y_cor_mat.png",plot=yym,width=12,height=10)


xym<-ggplot(data=d[d$variable1=='x'&d$variable2=='y',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("x-y cross correlation") + 
xlab("Base pair (x)")+ylab('Base pair (y)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_x_y_cor_mat.png",plot=xym,width=12,height=10)


xzm<-ggplot(data=d[d$variable1=='x'&d$variable2=='z',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("x-z cross correlation") + 
xlab("Base pair (x)")+ylab('Base pair (z)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_x_z_cor_mat.png",plot=xzm,width=12,height=10)


yzm<-ggplot(data=d[d$variable1=='y'&d$variable2=='z',],aes(x=BPnum1-94.5,y=BPnum2-94.5)) + ggtitle("y-z cross correlation") + 
xlab("Base pair (y)")+ylab('Base pair (z)')+geom_tile(aes(fill=corr))+scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+
scale_y_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_fill_gradient2(low='red',high='green',mid='white')+
theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+coord_fixed()

ggsave(filename="../analysis_data/dna2_param_y_z_cor_mat.png",plot=yzm,width=12,height=10)




#we want to output frames for Vitya
#Profiles of cross corr along sequence

prof_tr=d[d$variable1=='Twist'&d$variable2=='Roll'&d$BPnum1==d$BPnum2,c('BPnum1','corr')]
prof_tr=rename(prof_tr, c("BPnum1"='BPnum','corr'='<dTwist(i)*dRoll(i)>'))

prof_ts=d[d$variable1=='Twist'&d$variable2=='Slide'&d$BPnum1==d$BPnum2,c('BPnum1','corr')]
prof_ts=rename(prof_ts, c("BPnum1"='BPnum','corr'='<dTwist(i)*dSlide(i)>'))

prof_sr=d[d$variable1=='Slide'&d$variable2=='Roll'&d$BPnum1==d$BPnum2,c('BPnum1','corr')]
prof_sr=rename(prof_sr, c("BPnum1"='BPnum','corr'='<dSlide(i)*dRoll(i)>'))
profs_cc=merge(prof_tr,prof_ts)
profs_cc=merge(profs_cc,prof_sr)

write.csv(profs_cc, file="../analysis_data/dna2_param_tsr_cros_corr_prof.csv")

#matrix
mat_tt=d[d$variable1=='Twist'&d$variable2=='Twist',c('BPnum1','BPnum2','corr')]
mat_tt=rename(mat_tt, c('corr'='<dTwist(i)*dTwist(j)>'))

mat_rr=d[d$variable1=='Roll'&d$variable2=='Roll',c('BPnum1','BPnum2','corr')]
mat_rr=rename(mat_rr, c('corr'='<dRoll(i)*dRoll(j)>'))

mat_ss=d[d$variable1=='Slide'&d$variable2=='Slide',c('BPnum1','BPnum2','corr')]
mat_ss=rename(mat_ss, c('corr'='<dSlide(i)*dSlide(j)>'))

mat=merge(mat_tt,mat_rr)
mat=merge(mat,mat_ss)

write.csv(mat, file="../analysis_data/dna2_param_tsr_cros_corr_mat.csv")

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






j<-ggplot(data=subset(tot_md,INT_chain=='CHJ')) + ggtitle("Number of iteractions between DNA and ions, total (C+IP+HB), chain J: MD") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))
j<-j+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=3)



####


dna_prot<-subset(dna_prot,type!='Z')
dna_prot_cryst<-subset(dna_prot_cryst,type!='Z')
#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('IP','HB','VdW','IM','WM')

#Let's assign level order
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)




####General statistics section:
theme_set(theme_gray(base_size = 15))

##########Histograms with DNA_part classification

#---Cryst
a2<-ggplot(data=dna_prot_cryst[dna_prot_cryst$PROT_resname!='LYS' & dna_prot_cryst$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=DNA_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple", "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("blue", "purple", "dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot_cryst,aes(x=PROT_resname,alpha=DNA_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("red", "blue", "purple", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in crystal')
a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart_cryst.png",plot=a1,width=10,height=5)

##------MD

a2<-ggplot(data=dna_prot[dna_prot$PROT_resname!='LYS' & dna_prot$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=DNA_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c( "blue", "purple", "grey", "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c( "blue", "purple", "grey", "dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot,aes(x=PROT_resname,alpha=DNA_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("red", "blue", "purple", "grey", "dark green"))+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in MD simulations (average count)')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=100,ymax=620)
ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart.png",plot=a1,width=10,height=5)


###Histograms with PROT_part classification
#---Cryst
a2<-ggplot(data=dna_prot_cryst[dna_prot_cryst$PROT_resname!='LYS' & dna_prot_cryst$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple", "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("blue", "purple", "dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot_cryst,aes(x=PROT_resname,alpha=PROT_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("red", "blue", "purple", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in crystal')
a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
ggsave(filename="../analysis_data/int_dna_prot_hist_protpart_cryst.png",plot=a1,width=10,height=5)

##------MD

a2<-ggplot(data=dna_prot[dna_prot$PROT_resname!='LYS' & dna_prot$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c( "blue", "purple", "grey", "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c( "blue", "purple", "grey", "dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot,aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("red", "blue", "purple", "grey", "dark green"))+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in MD simulations (average count)')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=100,ymax=620)
ggsave(filename="../analysis_data/int_dna_prot_hist_protpart.png",plot=a1,width=10,height=5)


####Histogram of contacts with bases
theme_set(theme_gray(base_size = 12))

#---Cryst
c<-ggplot(data=dna_prot_cryst[dna_prot_cryst$DNA_part=='base',],aes(x=PROT_resname,alpha=PROT_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple",  "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("blue", "purple", "dark green"))
c<-c+ggtitle('Interactions between DNA bases and protein in crystal')

##------MD

md<-ggplot(data=dna_prot[dna_prot$DNA_part=='base',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "purple",'grey', "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("blue", "purple",'grey', "dark green"))
md<-md+ggtitle('Interactions between DNA bases and protein in MD simulations (average count)')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
q<-arrangeGrob(c,md)
ggsave(filename="../analysis_data/int_dna_base_prot_hist.png",plot=q,width=10,height=5)




###Histograms of contact types
#---Cryst

a1<-ggplot(data=dna_prot_cryst[dna_prot_cryst$type=='VdW',],aes(x=PROT_resname,fill=C_type))+
geom_bar(aes(color=C_type,y=..count..),position='stack')+
scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("red", "blue", "purple", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in crystal: contact types')
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
ggsave(filename="../analysis_data/int_dna_prot_hist_ctype_cryst.png",plot=a1,width=10,height=5)

##------MD

a1<-ggplot(data=dna_prot[dna_prot$type=='VdW',],aes(x=PROT_resname,fill=C_type,weight=av_num))+
geom_bar(aes(color=C_type,y=..count..),position='stack')+
scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c("red", "blue", "purple", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in MD: contact types')
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
ggsave(filename="../analysis_data/int_dna_prot_hist_ctype.png",plot=a1,width=10,height=5)

##----


##############################
#Now let's make profiles along DNA sequence

# theme_set(theme_grey(base_size=30))

# dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
# seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
# seqdf_t<-seqdf
# seqdf_t$X<-seqdf_t$X-73.0
# seqdf_t$Y<-rep(0,length(seqdf_t$X))
# # #Let;s average interactions for every resid.

################Total number of interactions
q=seq(73,-73,-1)

###Crystal
d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","type"),function(df) c(num=nrow(df)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
#Let's add zeros also

d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
t=data.frame(type=c('VdW','WM','IP','IM','HB'))
d1z=merge(d1zt,t)

d1<-rbind(d1i,d1j,d1z)

tot_cryst=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))

c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('VdW','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
o<-ggplot(data=tot_cryst[tot_cryst$type %in% c('IP','HB'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple", "grey","dark green"))

# q<-arrangeGrob(c,o)
# ggsave(filename="../analysis_data/int_dna_prot_prof.png",plot=q,width=15,height=5)


##MD----------
d1=ddply(dna_prot,c("DNA_chain","DNA_resid","type"),function(df) c(num=sum(df$av_num)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
d1<-rbind(d1i,d1j,d1z)
tot_md=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))


cmd<-ggplot(data=tot_md[tot_md$type %in% c('VdW','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
omd<-ggplot(data=tot_md[tot_md$type %in% c('IP','HB','IM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple","grey", "dark green"))


q<-arrangeGrob(c,cmd,o,omd,ncol=1)
ggsave(filename="../analysis_data/int_dna_prot_gen_prof.png",plot=q,width=15,height=10)


################Profile of interactions with bases, sugars, phosphates from MD
theme_set(theme_gray(base_size = 18))


d1=ddply(dna_prot,c("DNA_chain","DNA_resid","type","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid","type",'DNA_part','DNA_chain'),summarize,number=sum(num))

#CHAIN I
#By int type
chi_dna_part_cwm<-ggplot(data=int_md[(int_md$type %in% c('VdW','WM')) & (int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, VdW & WM") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
facet_grid(DNA_part~.,scales='free')

chi_dna_part_iphbim<-ggplot(data=int_md[(int_md$type %in% c('IP','HB','IM')) & (int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, IP & HB & IM") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple",'grey'))+
facet_grid(DNA_part~.,scales='free')

ggsave(filename="../analysis_data/int_dna_prot_dp_chi_cwm.png",plot=chi_dna_part_cwm,width=15,height=10)
ggsave(filename="../analysis_data/int_dna_prot_dp_chi_iphbim.png",plot=chi_dna_part_iphbim,width=15,height=10)

#By contact type

d1=ddply(subset(dna_prot,type=='VdW'),c("DNA_chain","DNA_resid","C_type","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1zp=merge(d1zt,t)
t2=data.frame(C_type=c('PN','PP','NN'))
d1z=merge(d1zp,t2)

d1<-rbind(d1,d1z)
int_md=ddply(d1,c("DNA_resid","DNA_chain",'DNA_part','C_type'),summarize,number=sum(num))

chi_dna_part_c_type<-ggplot(data=int_md[(int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, contacts by type") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=C_type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_chi_c_type.png",plot=chi_dna_part_c_type,width=15,height=10)


####Let's specifically study ARGinines

#MD
d1=ddply(subset(dna_prot,PROT_resname=='ARG' & type %in% c('IP','HB','VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (VdW+IP+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_arg.png",plot=chi_dna_part_arg_int,width=15,height=10)

#Crystal
d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG' & type %in% c('IP','HB','VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (VdW+IP+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_arg_cryst.png",plot=chi_dna_part_arg_int_cryst,width=15,height=10)


####Let's specifically study LYSinines


d1=ddply(subset(dna_prot,PROT_resname=='LYS' & type %in% c('IP','HB','VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_lys_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with LYS (VdW+IP+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_lys.png",plot=chi_dna_part_lys_int,width=15,height=10)



####Misc


####Let's specifically study ARGinines

#MD
d1=ddply(subset(dna_prot,PROT_resname=='ARG' & type %in% c('VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (contacts only)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')

#Crystal
d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG' & type %in% c('VdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (contacts only)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
q<-arrangeGrob(chi_dna_part_arg_int,chi_dna_part_arg_int_cryst,ncol=1)

ggsave(filename="../analysis_data/int_dna_prot_dp_arg_md_vs_cryst.png",plot=q,width=15,height=20)


