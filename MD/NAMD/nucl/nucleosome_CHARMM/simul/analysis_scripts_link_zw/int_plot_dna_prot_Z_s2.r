#R-script analysis of dna protein interactions
#Additional analysis asked by Zhurkin
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
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
dna_prot_cryst<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_avr_df.csv')

dna_prot<-subset(dna_prot,type=='Z')
dna_prot_cryst<-subset(dna_prot_cryst,type=='Z')
# dna_prot_cryst<-subset(dna_prot_cryst,param1<3.5)

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
# type_lev=c('IP','HB','VdW','IM','WM')

#Let's assign level order
# dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
# dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))

dna_prot$groove<-factor(dna_prot$groove,levels=c('minor','major'))
dna_prot_cryst$groove<-factor(dna_prot_cryst$groove,levels=c('minor','major'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)

dna_prot$PROT_part<-factor(dna_prot$PROT_part,levels=c('side chain','backbone'))
dna_prot_cryst$PROT_part<-factor(dna_prot_cryst$PROT_part,levels=c('side chain','backbone'))


#canonical 14 arginines
# resname ARG and segname CHA CHE and resid 83 63
# resname ARG and segname CHB CHF and resid 45
# resname ARG and segname CHC CHG and resid 42 77
# resname ARG and segname CHD CHH and resid 30
#(resname ARG and segname CHA CHE and resid 49)
#let's construct mask for them
arg_mask=with(dna_prot,((PROT_chain %in% c('CHA','CHE'))&(PROT_resid %in% c(83,63,49)))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid %in% c(45)))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid %in% c(42,77)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid %in% c(30))))
arg_mask_cryst=with(dna_prot_cryst,((PROT_chain %in% c('CHA','CHE'))&(PROT_resid %in% c(83,63,49)))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid %in% c(45)))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid %in% c(42,77)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid %in% c(30))))

####General statistics section:
theme_set(theme_gray(base_size = 15))

##########Histograms with DNA_part classification



####Histogram of contacts with bases
theme_set(theme_gray(base_size = 12))


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


################Profile of interactions with bases, sugars, phosphates from MD
theme_set(theme_bw(base_size = 18))


d1=ddply(dna_prot,c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))

int_md_vse=int_md
#CHAIN I
#By int type

#The same for crystal

d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_c=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_c_vse=int_c

#Comparative
int_c$data=factor('crystal')
int_c_all$data=factor('crystal')
int_md$data=factor('MD')
int_md_all$data=factor('MD')

int=rbind(int_c,int_md)
int_all=rbind(int_c_all,int_md_all)


#Let's do all the same for arginines


d1=ddply(subset(dna_prot,PROT_resname=='ARG'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_md_ARG=int_md

#The same for crystal

d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_c=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_c_ARG=int_c

#Comparative
int_c$data=factor('crystal')
int_c_all$data=factor('crystal')
int_md$data=factor('MD')
int_md_all$data=factor('MD')

int=rbind(int_c,int_md)
int_all=rbind(int_c_all,int_md_all)

####Let's subdivide all arg and canonical arg


d1=ddply(subset(dna_prot,arg_mask),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_md=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_md_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_md_ARG14=int_md

int_md_vse$data=factor('all')
int_md_ARG$data=factor('ARG')
int_md_ARG14$data=factor('ARG14')
int=rbind(int_md_vse,int_md_ARG,int_md_ARG14)

int[int$DNA_chain=='CHJ','DNA_resid']=int[int$DNA_chain=='CHJ','DNA_resid']*(-1)

q<-ggplot(data=int[(int$DNA_chain %in% c('CHI','CHJ')) & int$data=='ARG14' &int$DNA_part=='sugar',],aes(x=DNA_resid,y=number,color=DNA_chain)) + ggtitle("MD: DNA-protein interactions,  Z-contacts") + 
xlab("Base pair")+geom_line(size=1)+geom_point(subset = .(number > 0),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red','magenta'))+
facet_grid(DNA_part~.,scales='free')+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())

ggsave(filename="../analysis_data/int_dna_prot_Z_arg_comb.png",plot=q,width=15,height=4)

#The same for crystal



d1=ddply(subset(dna_prot_cryst,arg_mask_cryst),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int_c=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,number=sum(num))
int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,number=sum(num))
int_c_ARG14=int_c

int_c_vse$data=factor('all')
int_c_ARG$data=factor('ARG')
int_c_ARG14$data=factor('ARG14')
int=rbind(int_c_vse,int_c_ARG,int_c_ARG14)

int[int$DNA_chain=='CHJ','DNA_resid']=int[int$DNA_chain=='CHJ','DNA_resid']*(-1)

q<-ggplot(data=int[(int$DNA_chain %in% c('CHI','CHJ')) & int$data=='ARG14' &int$DNA_part=='sugar',],aes(x=DNA_resid,y=number,color=DNA_chain)) + ggtitle("X-ray: DNA-protein interactions,  Z-contacts 3.9A") + 
xlab("Base pair")+geom_line(size=1)+geom_point(subset = .(number > 0),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red','magenta'))+
facet_grid(DNA_part~.,scales='free')+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())



ggsave(filename="../analysis_data/int_dna_prot_Z_arg_comb_cryst.png",plot=q,width=15,height=4)

#minimum distance


d1=ddply(subset(dna_prot_cryst,arg_mask_cryst & PROT_atom=='CZ'),c("DNA_chain","DNA_resid","DNA_part"),function(df) c(min=min(df$param1)))
#Add zero frames
t=data.frame(DNA_part=c('base','sugar','phosphate'))
d1z=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),min=c(0))

d1zp=merge(d1z,t)
d2zp=d1zp
d2zp[,'DNA_chain']=factor(c('CHJ'))
d1<-rbind(d1,d1zp,d2zp)
int=ddply(d1,c("DNA_resid",'DNA_part','DNA_chain'),summarize,min_dist=sum(min))
# int_c_all=ddply(d1,c("DNA_resid",'DNA_chain'),summarize,min=sum(min))
# int_c_ARG14=int_c


int[int$DNA_chain=='CHJ','DNA_resid']=int[int$DNA_chain=='CHJ','DNA_resid']*(-1)

q<-ggplot(data=int[(int$DNA_chain %in% c('CHI','CHJ')) &(int$DNA_part %in% c('sugar','phosphate')),],aes(x=DNA_resid,y=min_dist,color=DNA_chain)) + ggtitle("X-ray: DNA-protein interactions, minimum distance to CZ of canonical ARG") + 
xlab("Base pair")+geom_line(size=1)+geom_point(subset = .(min_dist > 0),size=3)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red','magenta'))+
facet_grid(DNA_part~.,scales='free')+theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())



ggsave(filename="../analysis_data/int_dna_prot_Z_arg14_mindist_cryst.png",plot=q,width=15,height=6)


quit()
#By contact type

