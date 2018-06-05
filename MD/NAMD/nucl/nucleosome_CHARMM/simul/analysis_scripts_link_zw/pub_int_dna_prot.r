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
#core H3 > 43 < 132
#H4 > 23 < 99
#H2A > 15 < 118
#H2B > 29 < 124
dna_prot_cryst<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_avr_df.csv')

dna_prot<-subset(dna_prot,type=='SC')
dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')


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


