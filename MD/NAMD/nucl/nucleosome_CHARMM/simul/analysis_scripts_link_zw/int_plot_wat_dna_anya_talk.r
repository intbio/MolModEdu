#R-script analysis of water dna interactins 
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(gridExtra)
##############
###############
##Loading data frames
#DNA-protein interactions
# wat_dna<-read.csv('../analysis_data/wat_dna_raw_df.csv')
wat_dna_cryst<-read.csv('../analysis_data/wat_dna_avr_df_cryst.csv')
wat_dna<-read.csv('../analysis_data/wat_dna_avr_df.csv')

#Let's define levels order for better graphics
# prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('HB','C')

#Let's assign level order
wat_dna_cryst$type<-factor(wat_dna_cryst$type,levels=type_lev)
wat_dna$type<-factor(wat_dna$type,levels=type_lev)

wat_dna$DNA_part<-factor(wat_dna$DNA_part,levels=c('phosphate','sugar','base'))
wat_dna_cryst$DNA_part<-factor(wat_dna_cryst$DNA_part,levels=c('phosphate','sugar','base'))


#How to calculate number of water molecules
#remember that WAT_ind - identifies unique molecules


####General statistics section:
theme_set(theme_gray(base_size = 15))

##########Histograms of water molecules numbers

#Let's get for crystal number of waters in contact with DNA, with sugars, with phosphates,
#with bases

oc=subset(wat_dna_cryst,type=='C')
wnall=length(unique(oc$WAT_ind))
wnph=length(unique(subset(oc,DNA_part=='phosphate')$WAT_ind))
wns=length(unique(subset(oc,DNA_part=='sugar')$WAT_ind))
wnb=length(unique(subset(oc,DNA_part=='base')$WAT_ind))

wnumdf_c=data.frame(data='crystal',cutoff='3.9A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna,type=='C')
wnall=mean(ddply(oc,'Time',summarize,num=length(unique(WAT_ind)))$num)
wnph=mean(ddply(subset(oc,DNA_part=='phosphate'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wns=mean(ddply(subset(oc,DNA_part=='sugar'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnb=mean(ddply(subset(oc,DNA_part=='base'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnumdf_m=data.frame(data='MD',cutoff='3.9A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna_cryst,type=='C' & param1 < 3.5)
wnall=length(unique(oc$WAT_ind))
wnph=length(unique(subset(oc,DNA_part=='phosphate')$WAT_ind))
wns=length(unique(subset(oc,DNA_part=='sugar')$WAT_ind))
wnb=length(unique(subset(oc,DNA_part=='base')$WAT_ind))

wnumdf_c2=data.frame(data='crystal',cutoff='3.5A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna,type=='C' & param1 < 3.5)
wnall=mean(ddply(oc,'Time',summarize,num=length(unique(WAT_ind)))$num)
wnph=mean(ddply(subset(oc,DNA_part=='phosphate'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wns=mean(ddply(subset(oc,DNA_part=='sugar'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnb=mean(ddply(subset(oc,DNA_part=='base'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnumdf_m2=data.frame(data='MD',cutoff='3.5A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))


wnumdf=rbind(wnumdf_c,wnumdf_m,wnumdf_c2,wnumdf_m2)


#---Plot


#Let's try to discriminate between major and minor groove

# A: major: N6, C6, C5, N7, C8, minor: N1, C2,N3,C4
# G: major: O6, C6, C5, N7, C8, minor: N2, C2,N3,C4
# T: major: C6, C5, C4, C5M, O4 minor: N3, C2,O2
# C: major: C6, C5, C4, N4 minor: N3, C2,O2

#Major
wnall=0
wnph=0
wns=0
wnb=0

oc=subset(wat_dna_cryst,type=='C' & groove=='major')
wnall=length(unique(oc$WAT_ind))
wnb=length(unique(subset(oc,DNA_part=='base')$WAT_ind))

wnumdf_c=data.frame(data='crystal',cutoff='3.9A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna,type=='C'& groove=='major')
wnall=mean(ddply(oc,'Time',summarize,num=length(unique(WAT_ind)))$num)
wnb=mean(ddply(subset(oc,DNA_part=='base'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnumdf_m=data.frame(data='MD',cutoff='3.9A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna_cryst,type=='C' & param1 < 3.5& groove=='major')
wnall=length(unique(oc$WAT_ind))
wnb=length(unique(subset(oc,DNA_part=='base')$WAT_ind))

wnumdf_c2=data.frame(data='crystal',cutoff='3.5A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna,type=='C' & param1 < 3.5& groove=='major')
wnall=mean(ddply(oc,'Time',summarize,num=length(unique(WAT_ind)))$num)
wnb=mean(ddply(subset(oc,DNA_part=='base'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnumdf_m2=data.frame(data='MD',cutoff='3.5A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))


wnumdf=rbind(wnumdf_c,wnumdf_m,wnumdf_c2,wnumdf_m2)


#---Plot

# a1<-ggplot(data=wnumdf,aes(x=DNA_part,y=num_waters,fill=data,alpha=cutoff,color=data))+
# geom_bar(stat='identity',position='dodge')
# a1<-a1+ggtitle('Number of water molecules interacting with different parts of DNA in major groove:\n
# the criterion for interaction is distance between water oxygen and\n
# any of the heavy atoms of DNA (C,N,O,P)\n
# less than 3.5 or 3.9 A')

# ggsave(filename="../analysis_data/int_wat_dna_hist_dnapart_major.png",plot=a1,width=10,height=5)

major=wnumdf


#Minor


oc=subset(wat_dna_cryst,type=='C' & groove=='minor')
wnall=length(unique(oc$WAT_ind))
wnb=length(unique(subset(oc,DNA_part=='base')$WAT_ind))

wnumdf_c=data.frame(data='crystal',cutoff='3.9A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna,type=='C'& groove=='minor')
wnall=mean(ddply(oc,'Time',summarize,num=length(unique(WAT_ind)))$num)
wnb=mean(ddply(subset(oc,DNA_part=='base'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnumdf_m=data.frame(data='MD',cutoff='3.9A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna_cryst,type=='C' & param1 < 3.5& groove=='minor')
wnall=length(unique(oc$WAT_ind))
wnb=length(unique(subset(oc,DNA_part=='base')$WAT_ind))

wnumdf_c2=data.frame(data='crystal',cutoff='3.5A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))

oc=subset(wat_dna,type=='C' & param1 < 3.5& groove=='minor')
wnall=mean(ddply(oc,'Time',summarize,num=length(unique(WAT_ind)))$num)
wnb=mean(ddply(subset(oc,DNA_part=='base'),'Time',summarize,num=length(unique(WAT_ind)))$num)
wnumdf_m2=data.frame(data='MD',cutoff='3.5A',DNA_part=c('all','phosphate','sugar','base'),num_waters=c(wnall,wnph,wns,wnb))


wnumdf=rbind(wnumdf_c,wnumdf_m,wnumdf_c2,wnumdf_m2)


minor=wnumdf

minor$groove='minor'
major$groove='major'

wnumdf=rbind(minor,major)

#---Plot

# a1<-ggplot(data=wnumdf,aes(x=DNA_part,y=num_waters,fill=data,alpha=cutoff,color=data))+
# geom_bar(stat='identity',position='dodge')
# a1<-a1+ggtitle('Number of water molecules interacting with different parts of DNA in minor groove:\n
# the criterion for interaction is distance between water oxygen and\n
# any of the heavy atoms of DNA (C,N,O,P)\n
# less than 3.5 or 3.9 A')

# ggsave(filename="../analysis_data/int_wat_dna_hist_dnapart_minor.png",plot=a1,width=10,height=5)

a1<-ggplot(data=subset(wnumdf,data=='MD' & cutoff=='3.9A'),aes(x=groove,y=num_waters))+
geom_bar(stat='identity',position='dodge')
a1<-a1+ggtitle('Number of water molecules interacting with DNA grooves')+ylab('Water molecules number')

ggsave(filename="../analysis_data/int_wat_dna_hist_dnapart_anya_talk.png",plot=a1,width=10,height=5)

quit()





##------MD

# a2<-ggplot(data=wat_dna[wat_dna$PROT_resname!='LYS' & wat_dna$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=DNA_part,fill=type,weight=av_num))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c( "blue", "purple", "grey", "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c( "blue", "purple", "grey", "dark green"))+theme(legend.position="none")

# a1<-ggplot(data=wat_dna,aes(x=PROT_resname,alpha=DNA_part,fill=type,weight=av_num))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c("red", "blue", "purple", "grey", "dark green"))+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
# a1<-a1+ggtitle('Interactions between DNA and protein in MD simulations (average count)')
# # png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# # ggplot(data=wat_dna_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=100,ymax=620)
# ggsave(filename="../analysis_data/int_wat_dna_hist_dnapart.png",plot=a1,width=10,height=5)


# ###Histograms with PROT_part classification
# #---Cryst
# a2<-ggplot(data=wat_dna_cryst[wat_dna_cryst$PROT_resname!='LYS' & wat_dna_cryst$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c("blue", "purple", "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+
# scale_color_manual(values=c("blue", "purple", "dark green"))+theme(legend.position="none")

# a1<-ggplot(data=wat_dna_cryst,aes(x=PROT_resname,alpha=PROT_part,fill=type))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+
# scale_color_manual(values=c("red", "blue", "purple", "dark green"))
# a1<-a1+ggtitle('Interactions between DNA and protein in crystal')
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
# ggsave(filename="../analysis_data/int_wat_dna_hist_protpart_cryst.png",plot=a1,width=10,height=5)

# ##------MD

# a2<-ggplot(data=wat_dna[wat_dna$PROT_resname!='LYS' & wat_dna$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c( "blue", "purple", "grey", "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c( "blue", "purple", "grey", "dark green"))+theme(legend.position="none")

# a1<-ggplot(data=wat_dna,aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c("red", "blue", "purple", "grey", "dark green"))+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
# a1<-a1+ggtitle('Interactions between DNA and protein in MD simulations (average count)')
# # png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# # ggplot(data=wat_dna_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=100,ymax=620)
# ggsave(filename="../analysis_data/int_wat_dna_hist_protpart.png",plot=a1,width=10,height=5)


# ####Histogram of contacts with bases
# theme_set(theme_gray(base_size = 12))

# #---Cryst
# c<-ggplot(data=wat_dna_cryst[wat_dna_cryst$DNA_part=='base',],aes(x=PROT_resname,alpha=PROT_part,fill=type))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c("blue", "purple",  "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+
# scale_color_manual(values=c("blue", "purple", "dark green"))
# c<-c+ggtitle('Interactions between DNA bases and protein in crystal')

# ##------MD

# md<-ggplot(data=wat_dna[wat_dna$DNA_part=='base',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c("blue", "purple",'grey', "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+
# scale_color_manual(values=c("blue", "purple",'grey', "dark green"))
# md<-md+ggtitle('Interactions between DNA bases and protein in MD simulations (average count)')
# # png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# # ggplot(data=wat_dna_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
# q<-arrangeGrob(c,md)
# ggsave(filename="../analysis_data/int_dna_base_prot_hist.png",plot=q,width=10,height=5)




# ###Histograms of contact types
# #---Cryst

# a1<-ggplot(data=wat_dna_cryst[wat_dna_cryst$type=='C',],aes(x=PROT_resname,fill=C_type))+
# geom_bar(aes(color=C_type,y=..count..),position='stack')+
# scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+
# scale_color_manual(values=c("red", "blue", "purple", "dark green"))
# a1<-a1+ggtitle('Interactions between DNA and protein in crystal: contact types')
# # a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
# ggsave(filename="../analysis_data/int_wat_dna_hist_ctype_cryst.png",plot=a1,width=10,height=5)

# ##------MD

# a1<-ggplot(data=wat_dna[wat_dna$type=='C',],aes(x=PROT_resname,fill=C_type,weight=av_num))+
# geom_bar(aes(color=C_type,y=..count..),position='stack')+
# scale_fill_manual(values=c("red", "blue", "purple",  "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+
# scale_color_manual(values=c("red", "blue", "purple", "dark green"))
# a1<-a1+ggtitle('Interactions between DNA and protein in MD: contact types')
# # a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
# ggsave(filename="../analysis_data/int_wat_dna_hist_ctype.png",plot=a1,width=10,height=5)

# ##----


# ##############################
# #Now let's make profiles along DNA sequence

# # theme_set(theme_grey(base_size=30))

# # dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
# # seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
# # seqdf_t<-seqdf
# # seqdf_t$X<-seqdf_t$X-73.0
# # seqdf_t$Y<-rep(0,length(seqdf_t$X))
# # # #Let;s average interactions for every resid.

# ################Total number of interactions
# q=seq(73,-73,-1)

# ###Crystal
# d1=ddply(wat_dna_cryst,c("DNA_chain","DNA_resid","type"),function(df) c(num=nrow(df)))
# d1i=d1[d1$DNA_chain=='CHI',]
# d1j=d1[d1$DNA_chain=='CHJ',]
# d1j$DNA_resid<-q[d1j$DNA_resid+74]
# #Let's add zeros also

# d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
# t=data.frame(type=c('C','WM','IP','IM','HB'))
# d1z=merge(d1zt,t)

# d1<-rbind(d1i,d1j,d1z)

# tot_cryst=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))

# c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('C','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
# o<-ggplot(data=tot_cryst[tot_cryst$type %in% c('IP','HB'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple", "grey","dark green"))

# # q<-arrangeGrob(c,o)
# # ggsave(filename="../analysis_data/int_wat_dna_prof.png",plot=q,width=15,height=5)


# ##MD----------
# d1=ddply(wat_dna,c("DNA_chain","DNA_resid","type"),function(df) c(num=sum(df$av_num)))
# d1i=d1[d1$DNA_chain=='CHI',]
# d1j=d1[d1$DNA_chain=='CHJ',]
# d1j$DNA_resid<-q[d1j$DNA_resid+74]
# d1<-rbind(d1i,d1j,d1z)
# tot_md=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))


# cmd<-ggplot(data=tot_md[tot_md$type %in% c('C','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
# omd<-ggplot(data=tot_md[tot_md$type %in% c('IP','HB','IM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple","grey", "dark green"))


# q<-arrangeGrob(c,cmd,o,omd,ncol=1)
# ggsave(filename="../analysis_data/int_wat_dna_gen_prof.png",plot=q,width=15,height=10)


# ################Profile of interactions with bases, sugars, phosphates from MD
# theme_set(theme_gray(base_size = 18))


# d1=ddply(wat_dna,c("DNA_chain","DNA_resid","type","DNA_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# d1zp=merge(d1z,t)
# d2zp=d1zp
# d2zp[,'DNA_chain']=factor(c('CHJ'))
# d1<-rbind(d1,d1zp,d2zp)
# int_md=ddply(d1,c("DNA_resid","type",'DNA_part','DNA_chain'),summarize,number=sum(num))

# #CHAIN I
# #By int type
# chi_dna_part_cwm<-ggplot(data=int_md[(int_md$type %in% c('C','WM')) & (int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, C & WM") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
# facet_grid(DNA_part~.,scales='free')

# chi_dna_part_iphbim<-ggplot(data=int_md[(int_md$type %in% c('IP','HB','IM')) & (int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, IP & HB & IM") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple",'grey'))+
# facet_grid(DNA_part~.,scales='free')

# ggsave(filename="../analysis_data/int_wat_dna_dp_chi_cwm.png",plot=chi_dna_part_cwm,width=15,height=10)
# ggsave(filename="../analysis_data/int_wat_dna_dp_chi_iphbim.png",plot=chi_dna_part_iphbim,width=15,height=10)

# #By contact type

# d1=ddply(subset(wat_dna,type=='C'),c("DNA_chain","DNA_resid","C_type","DNA_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# d1zp=merge(d1zt,t)
# t2=data.frame(C_type=c('PN','PP','NN'))
# d1z=merge(d1zp,t2)

# d1<-rbind(d1,d1z)
# int_md=ddply(d1,c("DNA_resid","DNA_chain",'DNA_part','C_type'),summarize,number=sum(num))

# chi_dna_part_c_type<-ggplot(data=int_md[(int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, contacts by type") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=C_type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(DNA_part~.,scales='free')
# ggsave(filename="../analysis_data/int_wat_dna_dp_chi_c_type.png",plot=chi_dna_part_c_type,width=15,height=10)


# ####Let's specifically study ARGinines

# #MD
# d1=ddply(subset(wat_dna,PROT_resname=='ARG' & type %in% c('IP','HB','C') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
# d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

# chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (C+IP+HB)") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(DNA_part~.,scales='free')
# ggsave(filename="../analysis_data/int_wat_dna_dp_arg.png",plot=chi_dna_part_arg_int,width=15,height=10)

# #Crystal
# d1=ddply(subset(wat_dna_cryst,PROT_resname=='ARG' & type %in% c('IP','HB','C') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
# #Add zero frames
# d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
# d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

# chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (C+IP+HB)") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(DNA_part~.,scales='free')
# ggsave(filename="../analysis_data/int_wat_dna_dp_arg_cryst.png",plot=chi_dna_part_arg_int_cryst,width=15,height=10)


# ####Let's specifically study LYSinines


# d1=ddply(subset(wat_dna,PROT_resname=='LYS' & type %in% c('IP','HB','C') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
# d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

# chi_dna_part_lys_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with LYS (C+IP+HB)") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(DNA_part~.,scales='free')
# ggsave(filename="../analysis_data/int_wat_dna_dp_lys.png",plot=chi_dna_part_lys_int,width=15,height=10)



# ####Misc


# ####Let's specifically study ARGinines

# #MD
# d1=ddply(subset(wat_dna,PROT_resname=='ARG' & type %in% c('C') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
# d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

# chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (contacts only)") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(DNA_part~.,scales='free')

# #Crystal
# d1=ddply(subset(wat_dna_cryst,PROT_resname=='ARG' & type %in% c('C') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
# #Add zero frames
# d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
# d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
# t=data.frame(DNA_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

# chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (contacts only)") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(DNA_part~.,scales='free')
# q<-arrangeGrob(chi_dna_part_arg_int,chi_dna_part_arg_int_cryst,ncol=1)

# ggsave(filename="../analysis_data/int_wat_dna_dp_arg_md_vs_cryst.png",plot=q,width=15,height=20)


