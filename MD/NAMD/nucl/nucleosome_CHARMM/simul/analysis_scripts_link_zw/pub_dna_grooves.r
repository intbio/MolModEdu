#R-script analysis of dna protein interactions
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
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))
##############
###############
##Loading data frames
dna_avr<-read.csv('../analysis_data/dna2_param_df_md_avr.csv')
# dna_avr=subset(dna_avr,BPnum<168&BPnum>20)

dna_avr=rename(dna_avr, c('D12_av'='D12','D21_av'='D21','W12_av'='W12','W21_av'='W21','Shear_av'='Shear','Stretch_av'='Stretch','Stagger_av'='Stagger','Buckle_av'='Buckle','Prop.Tw_av'='Prop.Tw','Opening_av'='Opening','Shift_av'='Shift','Tilt_av'='Tilt','Rise_av'='Rise','Pairing_av'='Pairing',"x_av"='x',"y_av"='y',"z_av"='z', "Roll_av"="Roll",'Twist_av'='Twist','Slide_av'='Slide','P_1_av'='P_1','P_2_av'='P_2','chi_1_av'='chi_1','chi_2_av'='chi_2'))

dna_avr_cryst<-read.csv('../analysis_data/dna2_param_df_cryst.csv')
# dna_avr_cryst=rename(dna_avr_cryst, c('D12_av'='D12','D21_av'='D21','W12_av'='W12','W21_av'='W21','Shear_av'='Shear','Stretch_av'='Stretch','Stagger_av'='Stagger','Buckle_av'='Buckle','Prop.Tw_av'='Prop.Tw','Opening_av'='Opening','Shift_av'='Shift','Tilt_av'='Tilt','Rise_av'='Rise','Pairing_av'='Pairing',"x_av"='x',"y_av"='y',"z_av"='z', "Roll_av"="Roll",'Twist_av'='Twist','Slide_av'='Slide','P_1_av'='P_1','P_2_av'='P_2','chi_1_av'='chi_1','chi_2_av'='chi_2',))
dna_avr$BPnum=dna_avr$BPnum-74-20
dna_avr_cryst$BPnum=dna_avr_cryst$BPnum-74-20

dna_avr$turn=c('positive')
dna_avr[dna_avr$BPnum<0,'turn']=c('negative')
dna_avr[dna_avr$BPnum<0,'BPnum']=-dna_avr[dna_avr$BPnum<0,'BPnum']
t=subset(dna_avr,BPnum==0)
t$turn=c('negative')
dna_avr=rbind(dna_avr,t)

dna_avr_cryst$turn=c('positive')
dna_avr_cryst[dna_avr_cryst$BPnum<0,'turn']=c('negative')
dna_avr_cryst[dna_avr_cryst$BPnum<0,'BPnum']=-dna_avr_cryst[dna_avr_cryst$BPnum<0,'BPnum']
t=subset(dna_avr_cryst,BPnum==0)
t$turn=c('negative')
dna_avr_cryst=rbind(dna_avr_cryst,t)


d_c=melt(dna_avr_cryst,id.vars=c('BPnum','turn'),measure.vars=c('W12','W21','D12','D21'))
d_md=melt(dna_avr,id.vars=c('BPnum','turn'),measure.vars=c('W12','W21','D12','D21'))

d_c$data='X-ray'
d_md$data='MD average'
d=rbind(d_c,d_md)


d$turn<-factor(d$turn,levels=c('negative','positive'),labels=c('SHL < 0','SHL > 0'))

head(d)
q<-ggplot(data=subset(d,variable=='W12'),aes(x=BPnum/10,y=value,color=turn)) + ggtitle("DNA grooves, FN-model, average") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(data~.,scales='fixed')+ylab('Size, Å')+
scale_x_continuous(limits=c(-0.05,9.3),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Superhelix location (SHL)')+
# scale_alpha_discrete(range = c(0.3, 1.0))+
# scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
# scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
geom_line()

# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_dna_grooves.png",plot=q,width=15,height=5.5)



q()

dna_prot_max<-read.csv('../analysis_data/dna_prot_max_sc_df.csv')

#collapse by DNA symmetry and chains.
dna_prot_max[dna_prot_max$DNA_resid<0,'DNA_resid']=-dna_prot_max[dna_prot_max$DNA_resid<0,'DNA_resid']
head(dna_prot_max)
dna_prot_max$PROT_chain=revalue(dna_prot_max$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))
dna_prot_max=ddply(dna_prot_max,c('DNA_resid','PROT_chain'),summarize,max_sc=max(max_sc),core=core[which.max(max_sc)],tails=tails[which.max(max_sc)])
dna_prot_max_m=melt(dna_prot_max, id.vars=c('DNA_resid','PROT_chain'),measure.vars=c('core','tails'))
head(dna_prot_max_m)

dna_prot_max_m=arrange(dna_prot_max_m,variable,PROT_chain)


q<-ggplot(data=subset(dna_prot_max_m,value>0),aes(x=DNA_resid,y=value,fill=PROT_chain,color=PROT_chain)) + ggtitle("MAX DNA-protein interactions in nucleosome: MD, all") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_bar(stat='identity',position='stack',width=0.5)+
facet_grid(variable~.,scales='fixed')+
scale_x_continuous(limits=c(-0.5,93),breaks = round(seq(0,90, by = 10),1))+scale_alpha_discrete(range = c(0.3, 1.0))+
scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='yellow','H2B'='red'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='yellow','H2B'='red'),labels=c('H3','H4','H2A','H2B'),name='Histones')#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_max.png",plot=q,width=15,height=5)


q()
#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SC','SB','HB','vdW','IM','WM')
# head(subset(dna_prot_cryst,type=='SB'))
#Let's assign level order
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)


dna_prot_crysto=dna_prot_cryst #subset(dna_prot_cryst, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))

dna_proto=dna_prot #subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))



####Histogram of contacts with bases
theme_set(theme_gray(base_size = 12))

################Total number of interactions
q=seq(93,-93,-1)

###Crystal
d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","DNA_part",'PROT_chain',"type"),function(df) c(num=nrow(df)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
#Let's add zeros also

d1ztp=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('phosphate'),num=c(0))
d1zts=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('sugar'),num=c(0))
d1ztb=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('base'),num=c(0))
d1zt<-rbind(d1ztp,d1zts,d1ztb)


t=data.frame(type=c('vdW','WM','SB','IM','HB','SC'))
tch=data.frame(PROT_chain=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'))
d1z=merge(d1zt,t)
d1z=merge(d1z,tch)

d1<-rbind(d1i,d1j,d1z)

tot_cryst=ddply(d1,c("DNA_resid","DNA_part",'PROT_chain',"type"),summarize,number=sum(num))
c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('SC'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal, all") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(DNA_part~.,scales='free')+
scale_x_continuous(limits=c(-93,93),breaks = round(seq(-90,90, by = 10),1))+scale_color_manual(values=c("blue",'green','yellow','red',"dark blue",'dark green','orange','dark red'))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)

# q<-arrangeGrob(c,o)
# ggsave(filename="../analysis_data/int_dna_prot_prof.png",plot=q,width=15,height=5)

##MD----------
d1=ddply(dna_prot,c("DNA_chain","DNA_resid","DNA_part",'PROT_chain',"type"),function(df) c(num=sum(df$av_num)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
d1<-rbind(d1i,d1j,d1z)
tot_md=ddply(d1,c("DNA_resid","DNA_part",'PROT_chain',"type"),summarize,number=sum(num))


cmd<-ggplot(data=tot_md[tot_md$type %in% c('SC'),]) + ggtitle("DNA-protein interactions in nucleosome: MD, all") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(DNA_part~.,scales='free')+
scale_x_continuous(limits=c(-93,93),breaks = round(seq(-90,90, by = 10),1))+scale_color_manual(values=c("blue",'green','yellow','red',"dark blue",'dark green','orange','dark red'))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)


q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_all.png",plot=q,width=15,height=10)

#####Core


dna_prot_cryst=subset(dna_prot_crysto, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))

dna_prot=subset(dna_proto, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))



################Total number of interactions
q=seq(93,-93,-1)

###Crystal
d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","DNA_part",'PROT_chain',"type"),function(df) c(num=nrow(df)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
#Let's add zeros also

d1ztp=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('phosphate'),num=c(0))
d1zts=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('sugar'),num=c(0))
d1ztb=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('base'),num=c(0))
d1zt<-rbind(d1ztp,d1zts,d1ztb)


t=data.frame(type=c('vdW','WM','SB','IM','HB','SC'))
tch=data.frame(PROT_chain=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'))
d1z=merge(d1zt,t)
d1z=merge(d1z,tch)

d1<-rbind(d1i,d1j,d1z)

tot_cryst=ddply(d1,c("DNA_resid","DNA_part",'PROT_chain',"type"),summarize,number=sum(num))
c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('SC'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal, core") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(DNA_part~.,scales='free')+
scale_x_continuous(limits=c(-93,93),breaks = round(seq(-90,90, by = 10),1))+scale_color_manual(values=c("blue",'green','yellow','red',"dark blue",'dark green','orange','dark red'))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)

# q<-arrangeGrob(c,o)
# ggsave(filename="../analysis_data/int_dna_prot_prof.png",plot=q,width=15,height=5)

##MD----------
d1=ddply(dna_prot,c("DNA_chain","DNA_resid","DNA_part",'PROT_chain',"type"),function(df) c(num=sum(df$av_num)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
d1<-rbind(d1i,d1j,d1z)
tot_md=ddply(d1,c("DNA_resid","DNA_part",'PROT_chain',"type"),summarize,number=sum(num))


cmd<-ggplot(data=tot_md[tot_md$type %in% c('SC'),]) + ggtitle("DNA-protein interactions in nucleosome: MD, core") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(DNA_part~.,scales='free')+
scale_x_continuous(limits=c(-93,93),breaks = round(seq(-90,90, by = 10),1))+scale_color_manual(values=c("blue",'green','yellow','red',"dark blue",'dark green','orange','dark red'))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)


q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_core.png",plot=q,width=15,height=10)

###Tails



dna_prot_cryst=subset(dna_prot_crysto, !(((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))

dna_prot=subset(dna_proto, !(((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))



################Total number of interactions
q=seq(93,-93,-1)

###Crystal
d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","DNA_part",'PROT_chain',"type"),function(df) c(num=nrow(df)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
#Let's add zeros also

d1ztp=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('phosphate'),num=c(0))
d1zts=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('sugar'),num=c(0))
d1ztb=data.frame(DNA_resid=seq(93,-93,-1),DNA_chain=c('CHI'),DNA_part=c('base'),num=c(0))
d1zt<-rbind(d1ztp,d1zts,d1ztb)


t=data.frame(type=c('vdW','WM','SB','IM','HB','SC'))
tch=data.frame(PROT_chain=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'))
d1z=merge(d1zt,t)
d1z=merge(d1z,tch)

d1<-rbind(d1i,d1j,d1z)

tot_cryst=ddply(d1,c("DNA_resid","DNA_part",'PROT_chain',"type"),summarize,number=sum(num))
c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('SC'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal, tails") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(DNA_part~.,scales='free')+
scale_x_continuous(limits=c(-93,93),breaks = round(seq(-90,90, by = 10),1))+scale_color_manual(values=c("blue",'green','yellow','red',"dark blue",'dark green','orange','dark red'))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)

# q<-arrangeGrob(c,o)
# ggsave(filename="../analysis_data/int_dna_prot_prof.png",plot=q,width=15,height=5)

##MD----------
d1=ddply(dna_prot,c("DNA_chain","DNA_resid","DNA_part",'PROT_chain',"type"),function(df) c(num=sum(df$av_num)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
d1<-rbind(d1i,d1j,d1z)
tot_md=ddply(d1,c("DNA_resid","DNA_part",'PROT_chain',"type"),summarize,number=sum(num))


cmd<-ggplot(data=tot_md[tot_md$type %in% c('SC'),]) + ggtitle("DNA-protein interactions in nucleosome: MD, tails") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
facet_grid(DNA_part~.,scales='free')+
scale_x_continuous(limits=c(-93,93),breaks = round(seq(-90,90, by = 10),1))+scale_color_manual(values=c("blue",'green','yellow','red',"dark blue",'dark green','orange','dark red'))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)


q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_tails.png",plot=q,width=15,height=10)


quit()


