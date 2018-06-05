#R-script analysis of prot protein interactions
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
#PROT1-protein interactions
# prot_prot<-read.csv('../analysis_data/prot_prot_raw_df.csv')
prot_prot_cryst<-read.csv('../analysis_data/prot_prot_avr_df_cryst.csv')
prot_prot<-read.csv('../analysis_data/prot_prot_avr_df.csv')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')

# prot_prot_cryst=subset(prot_prot_cryst,type!='SC')
prot_prot=subset(prot_prot,type!='SC')

#exclude intrachain interactions

# prot_prot_cryst=subset(prot_prot_cryst,PROT1_chain!=PROT2_chain)
# prot_prot=subset(prot_prot,PROT1_chain!=PROT2_chain)


type_lev=c('SB','HB','vdW','IM','WM')

#residues important for H2A.Z
	# vdw.changeSelection("segname CHC and noh and resid 30 33 35 37 38 39 40 41 48 62 71 73 75 76 79 89 94 95 97 98 100 104 110 117")
#Let's concentrate on them
h2azimp=c(30,33,35,37,38,39,40,41,48,62,71,73,75,76,79,89,94,95,97,98,100,104,110,117)
prot_prot_cryst=subset(prot_prot_cryst,((PROT1_chain=='CHC')&(PROT1_resid %in% h2azimp))|((PROT2_chain=='CHC')&(PROT2_resid %in% h2azimp)))

prot_prot=subset(prot_prot,((PROT1_chain=='CHC')&(PROT1_resid %in% h2azimp))&((PROT2_chain=='CHG')&(PROT2_resid %in% h2azimp)))


#Let's assign level order
prot_prot_cryst$type<-factor(prot_prot_cryst$type,levels=type_lev)
prot_prot$type<-factor(prot_prot$type,levels=type_lev)

prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
prot_prot_cryst$PROT1_part<-factor(prot_prot_cryst$PROT1_part,levels=c('backbone','side chain'))


prot_prot_cryst$PROT1_resname<-factor(prot_prot_cryst$PROT1_resname,levels=prot_rn_lev)
prot_prot$PROT1_resname<-factor(prot_prot$PROT1_resname,levels=prot_rn_lev)

prot_prot_cryst$PROT2_resname<-factor(prot_prot_cryst$PROT2_resname,levels=prot_rn_lev)
prot_prot$PROT2_resname<-factor(prot_prot$PROT2_resname,levels=prot_rn_lev)



####General statistics section:
theme_set(theme_gray(base_size = 15))

##############################
#Now let's make a contact map
# Idea 1. is to make a cumulative map between histones
nn=data.frame(CHAIN_name=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'),HIST_name=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

t=ddply(prot_prot_cryst,c('PROT1_chain','PROT2_chain'),summarize,SB=table(type)['SB'],HB=table(type)['HB'],vdW=table(type)['vdW'],IM=table(type)['IM'],WM=table(type)['WM'])
t2=merge(t,nn,by.x=c('PROT1_chain'),by.y=c('CHAIN_name'))
t3=merge(t2,nn,by.x=c('PROT2_chain'),by.y=c('CHAIN_name'),suffixes=c('_1','_2'))
t4=subset(t3,select = !(colnames(t3) %in% c('PROT1_chain','PROT2_chain')))
t5=melt(t4,id=c('HIST_name_1','HIST_name_2'))
t6=transform(t5,X=1)
cm_c=rbind(t6,data.frame(HIST_name_1='H2B H',HIST_name_2='H3 A',variable='IM',value=NA,X=1))

cm_c$HIST_name_1<-factor(cm_c$HIST_name_1,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))
cm_c$HIST_name_2<-factor(cm_c$HIST_name_2,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))


a<-ggplot(data=cm_c,aes(x=X,y=value))+geom_bar(aes(x=X,y=value,color=variable,fill=variable),stat='identity',position='dodge')+
facet_grid(HIST_name_1~HIST_name_2,space='free')+
scale_fill_manual(values=c("red", "blue", "purple",'grey',  "dark green"))+
scale_color_manual(values=c("red", "blue", "purple",'grey', "dark green"))
a<-a+theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())
a<-a+ggtitle('Interactions between histones in nucleosome: crystal')
ggsave(filename="../analysis_data/int_prot_prot_sum_cryst_H2AZ.png",plot=a,width=7,height=7)


#MD---
t=ddply(prot_prot,c('PROT1_chain','PROT2_chain','type'),summarize,num=sum(av_num))
t2=merge(t,nn,by.x=c('PROT1_chain'),by.y=c('CHAIN_name'))
t3=merge(t2,nn,by.x=c('PROT2_chain'),by.y=c('CHAIN_name'),suffixes=c('_1','_2'))
t4=subset(t3,select = !(colnames(t3) %in% c('PROT1_chain','PROT2_chain')))
t6=transform(t4,X=1)
cm=rbind(t6,data.frame(HIST_name_1='H2B H',HIST_name_2='H3 A',type='IM',num=NA,X=1))

cm$HIST_name_1<-factor(cm$HIST_name_1,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))
cm$HIST_name_2<-factor(cm$HIST_name_2,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

cm$type<-factor(cm$type,levels=type_lev)

a<-ggplot(data=cm,aes(x=X,y=num))+geom_bar(aes(x=X,y=num,color=type,fill=type),stat='identity',position='dodge')+
facet_grid(HIST_name_1~HIST_name_2,space='free')+
scale_fill_manual(values=c("red", "blue", "purple",'grey',  "dark green"))+
scale_color_manual(values=c("red", "blue", "purple",'grey', "dark green"))
a<-a+theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())
a<-a+ggtitle('Interactions between histones in nucleosome: MD')
ggsave(filename="../analysis_data/int_prot_prot_sum_H2AZ.png",plot=a,width=7,height=7)



# Idea 2. Map along all sequences.

#Should correspond to cross-correlatoin maps
#Cryst ----
nn=data.frame(CHAIN_name=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'),HIST_name=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

t=ddply(prot_prot_cryst,c('PROT1_chain','PROT2_chain','PROT1_resid','PROT2_resid'),summarize,num=table(type)['SB']+table(type)['HB']+table(type)['vdW'])
t2=merge(t,nn,by.x=c('PROT1_chain'),by.y=c('CHAIN_name'))
t3=merge(t2,nn,by.x=c('PROT2_chain'),by.y=c('CHAIN_name'),suffixes=c('_1','_2'))
t4=subset(t3,select = !(colnames(t3) %in% c('PROT1_chain','PROT2_chain')))
# t5=melt(t4,id=c('HIST_name_1','HIST_name_2'))
# t6=transform(t5,X=1)
# c_c=rbind(t4,data.frame(HIST_name_1='H2B H',HIST_name_2='H3 A',PROT1_resid=60,PROT2_resid=60,num=NA))
c_c=t4

c_c$HIST_name_1<-factor(c_c$HIST_name_1,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))
c_c$HIST_name_2<-factor(c_c$HIST_name_2,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

c_cr=c_c
names(c_cr)<-c("PROT1_resid" ,"PROT2_resid", "num","HIST_name_2","HIST_name_1")
c_c=rbind(c_c,c_cr)

a<-ggplot(data=c_c,aes(x=PROT1_resid,y=PROT2_resid))+geom_point(aes(color=num))+
facet_grid(HIST_name_1~HIST_name_2)+scale_colour_gradient(low="blue", high="red")

a<-a+ggtitle('Interactions between histones in nucleosome, num = SB+HB+vdW: crystal')
ggsave(filename="../analysis_data/int_prot_prot_full_cryst_H2AZ.png",plot=a,width=7,height=7)


#MD-----

nn=data.frame(CHAIN_name=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'),HIST_name=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

t=ddply(prot_prot[prot_prot$type %in% c('vdW','HB','SB'),],c('PROT1_chain','PROT2_chain','PROT1_resid','PROT2_resid'),summarize,num=sum(av_num))
t2=merge(t,nn,by.x=c('PROT1_chain'),by.y=c('CHAIN_name'))
t3=merge(t2,nn,by.x=c('PROT2_chain'),by.y=c('CHAIN_name'),suffixes=c('_1','_2'))
t4=subset(t3,select = !(colnames(t3) %in% c('PROT1_chain','PROT2_chain')))
# t5=melt(t4,id=c('HIST_name_1','HIST_name_2'))
# t6=transform(t5,X=1)
# c_md=rbind(t4,data.frame(HIST_name_1='H2B H',HIST_name_2='H3 A',PROT1_resid=60,PROT2_resid=60,num=NA))
c_md=t4

c_md$HIST_name_1<-factor(c_md$HIST_name_1,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))
c_md$HIST_name_2<-factor(c_md$HIST_name_2,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

c_mdr=c_md
names(c_mdr)<-c("PROT1_resid" ,"PROT2_resid", "num","HIST_name_2","HIST_name_1")
c_md=rbind(c_md,c_mdr)

a<-ggplot(data=c_md,aes(x=PROT1_resid,y=PROT2_resid))+geom_point(aes(color=num))+
facet_grid(HIST_name_1~HIST_name_2)+scale_colour_gradient(low="blue", high="red")

a<-a+ggtitle('Interactions between histones in nucleosome, num = SB+HB+vdW: MD')
ggsave(filename="../analysis_data/int_prot_prot_full_H2AZ.png",plot=a,width=7,height=7)


#Maps along seqs with sequences

# # theme_set(theme_grey(base_size=30))

# # protseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
# # seqdf=data.frame(sequence=substring(protseq, seq(1,nchar(protseq)),seq(1,nchar(protseq))),X=seq(0,nchar(protseq)-1))
# # seqdf_t<-seqdf
# # seqdf_t$X<-seqdf_t$X-73.0
# # seqdf_t$Y<-rep(0,length(seqdf_t$X))
# # # #Let;s average interactions for every resid.

# ################Total number of interactions
# q=seq(73,-73,-1)

# ###Crystal
# d1=ddply(prot_prot_cryst,c("PROT1_chain","PROT1_resid","type"),function(df) c(num=nrow(df)))
# d1i=d1[d1$PROT1_chain=='CHI',]
# d1j=d1[d1$PROT1_chain=='CHJ',]
# d1j$PROT1_resid<-q[d1j$PROT1_resid+74]
# #Let's add zeros also

# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),PROT1_chain=c('CHI'),num=c(0))
# t=data.frame(type=c('VdW','WM','IP','IM','HB'))
# d1z=merge(d1zt,t)

# d1<-rbind(d1i,d1j,d1z)

# tot_cryst=ddply(d1,c("PROT1_resid","type"),summarize,number=sum(num))

# c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('VdW','WM'),]) + ggtitle("PROT1-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
# o<-ggplot(data=tot_cryst[tot_cryst$type %in% c('IP','HB'),]) + ggtitle("PROT1-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple", "grey","dark green"))

# # q<-arrangeGrob(c,o)
# # ggsave(filename="../analysis_data/int_prot_prot_prof.png",plot=q,width=15,height=5)


# ##MD----------
# d1=ddply(prot_prot,c("PROT1_chain","PROT1_resid","type"),function(df) c(num=sum(df$av_num)))
# d1i=d1[d1$PROT1_chain=='CHI',]
# d1j=d1[d1$PROT1_chain=='CHJ',]
# d1j$PROT1_resid<-q[d1j$PROT1_resid+74]
# d1<-rbind(d1i,d1j,d1z)
# tot_md=ddply(d1,c("PROT1_resid","type"),summarize,number=sum(num))


# cmd<-ggplot(data=tot_md[tot_md$type %in% c('VdW','WM'),]) + ggtitle("PROT1-protein interactions in nucleosome: MD") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
# omd<-ggplot(data=tot_md[tot_md$type %in% c('IP','HB','IM'),]) + ggtitle("PROT1-protein interactions in nucleosome: MD") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple","grey", "dark green"))


# q<-arrangeGrob(c,cmd,o,omd,ncol=1)
# ggsave(filename="../analysis_data/int_prot_prot_gen_prof.png",plot=q,width=15,height=10)


# ################Profile of interactions with bases, sugars, phosphates from MD
# theme_set(theme_gray(base_size = 18))


# d1=ddply(prot_prot,c("PROT1_chain","PROT1_resid","type","PROT1_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# d1zp=merge(d1z,t)
# d2zp=d1zp
# d2zp[,'PROT1_chain']=factor(c('CHJ'))
# d1<-rbind(d1,d1zp,d2zp)
# int_md=ddply(d1,c("PROT1_resid","type",'PROT1_part','PROT1_chain'),summarize,number=sum(num))

# #CHAIN I
# #By int type
# chi_prot_part_cwm<-ggplot(data=int_md[(int_md$type %in% c('VdW','WM')) & (int_md$PROT1_chain=='CHI'),]) + ggtitle("PROT1-protein interactions in nucleosome, MD, chain I, C & WM") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
# facet_grid(PROT1_part~.,scales='free')

# chi_prot_part_iphbim<-ggplot(data=int_md[(int_md$type %in% c('IP','HB','IM')) & (int_md$PROT1_chain=='CHI'),]) + ggtitle("PROT1-protein interactions in nucleosome, MD, chain I, IP & HB & IM") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple",'grey'))+
# facet_grid(PROT1_part~.,scales='free')

# ggsave(filename="../analysis_data/int_prot_prot_dp_chi_cwm.png",plot=chi_prot_part_cwm,width=15,height=10)
# ggsave(filename="../analysis_data/int_prot_prot_dp_chi_iphbim.png",plot=chi_prot_part_iphbim,width=15,height=10)

# #By contact type

# d1=ddply(subset(prot_prot,type=='VdW'),c("PROT1_chain","PROT1_resid","C_type","PROT1_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),PROT1_chain=c('CHI'),num=c(0))
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# d1zp=merge(d1zt,t)
# t2=data.frame(C_type=c('PN','PP','NN'))
# d1z=merge(d1zp,t2)

# d1<-rbind(d1,d1z)
# int_md=ddply(d1,c("PROT1_resid","PROT1_chain",'PROT1_part','C_type'),summarize,number=sum(num))

# chi_prot_part_c_type<-ggplot(data=int_md[(int_md$PROT1_chain=='CHI'),]) + ggtitle("PROT1-protein interactions in nucleosome, MD, chain I, contacts by type") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=C_type),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(PROT1_part~.,scales='free')
# ggsave(filename="../analysis_data/int_prot_prot_dp_chi_c_type.png",plot=chi_prot_part_c_type,width=15,height=10)


# ####Let's specifically study ARGinines

# #MD
# d1=ddply(subset(prot_prot,PROT_resname=='ARG' & type %in% c('IP','HB','VdW') ),c('PROT1_chain',"PROT1_resid","PROT1_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHI'))
# d2zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHJ'))
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('PROT1_chain',"PROT1_resid",'PROT1_part'),summarize,number=sum(num))

# chi_prot_part_arg_int<-ggplot(data=int_md) + ggtitle("PROT1-protein interactions in nucleosome, MD, interactions with ARG (C+IP+HB)") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=PROT1_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(PROT1_part~.,scales='free')
# ggsave(filename="../analysis_data/int_prot_prot_dp_arg.png",plot=chi_prot_part_arg_int,width=15,height=10)

# #Crystal
# d1=ddply(subset(prot_prot_cryst,PROT_resname=='ARG' & type %in% c('IP','HB','VdW') ),c('PROT1_chain',"PROT1_resid","PROT1_part"),function(df) c(num=length(df$param1)))
# #Add zero frames
# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHI'))
# d2zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHJ'))
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('PROT1_chain',"PROT1_resid",'PROT1_part'),summarize,number=sum(num))

# chi_prot_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("PROT1-protein interactions in nucleosome, Crystal, interactions with ARG (C+IP+HB)") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=PROT1_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(PROT1_part~.,scales='free')
# ggsave(filename="../analysis_data/int_prot_prot_dp_arg_cryst.png",plot=chi_prot_part_arg_int_cryst,width=15,height=10)


# ####Let's specifically study LYSinines


# d1=ddply(subset(prot_prot,PROT_resname=='LYS' & type %in% c('IP','HB','VdW') ),c('PROT1_chain',"PROT1_resid","PROT1_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHI'))
# d2zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHJ'))
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('PROT1_chain',"PROT1_resid",'PROT1_part'),summarize,number=sum(num))

# chi_prot_part_lys_int<-ggplot(data=int_md) + ggtitle("PROT1-protein interactions in nucleosome, MD, interactions with LYS (C+IP+HB)") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=PROT1_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(PROT1_part~.,scales='free')
# ggsave(filename="../analysis_data/int_prot_prot_dp_lys.png",plot=chi_prot_part_lys_int,width=15,height=10)



# ####Misc


# ####Let's specifically study ARGinines

# #MD
# d1=ddply(subset(prot_prot,PROT_resname=='ARG' & type %in% c('VdW') ),c('PROT1_chain',"PROT1_resid","PROT1_part"),function(df) c(num=sum(df$av_num)))
# #Add zero frames
# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHI'))
# d2zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHJ'))
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('PROT1_chain',"PROT1_resid",'PROT1_part'),summarize,number=sum(num))

# chi_prot_part_arg_int<-ggplot(data=int_md) + ggtitle("PROT1-protein interactions in nucleosome, MD, interactions with ARG (contacts only)") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=PROT1_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(PROT1_part~.,scales='free')

# #Crystal
# d1=ddply(subset(prot_prot_cryst,PROT_resname=='ARG' & type %in% c('VdW') ),c('PROT1_chain',"PROT1_resid","PROT1_part"),function(df) c(num=length(df$param1)))
# #Add zero frames
# d1zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHI'))
# d2zt=data.frame(PROT1_resid=seq(73,-73,-1),num=c(0),PROT1_chain=c('CHJ'))
# t=data.frame(PROT1_part=c('base','sugar','phosphate'))
# dzt<-rbind(d1zt,d2zt)
# dz=merge(dzt,t)
# d1=rbind(d1,dz)

# int_md=ddply(d1,c('PROT1_chain',"PROT1_resid",'PROT1_part'),summarize,number=sum(num))

# chi_prot_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("PROT1-protein interactions in nucleosome, Crystal, interactions with ARG (contacts only)") + 
# xlab("Base pair")+geom_line(aes(x=PROT1_resid,y=number,color=PROT1_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
# facet_grid(PROT1_part~.,scales='free')
# q<-arrangeGrob(chi_prot_part_arg_int,chi_prot_part_arg_int_cryst,ncol=1)

# ggsave(filename="../analysis_data/int_prot_prot_dp_arg_md_vs_cryst.png",plot=q,width=15,height=20)


