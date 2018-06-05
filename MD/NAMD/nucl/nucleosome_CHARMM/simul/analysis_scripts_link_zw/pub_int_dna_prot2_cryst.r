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
# dna_prot_cryst<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')
dna_prot_avr<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')



# q()


dna_prot<-subset(dna_prot_avr,type=='SC')
# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))
dna_prot_tails=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid <= 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid <= 23))|((PROT_chain %in% c('CHC','CHG'))&((PROT_resid <= 15)|(PROT_resid >= 119)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid <= 29))))

dna_prot_core_temp=ddply(dna_prot_core,c('DNA_resid','PROT_chain'),summarize,tot_num_core=length(param1))
dna_prot_tails_temp=ddply(dna_prot_tails,c('DNA_resid','PROT_chain'),summarize,tot_num_tails=length(param1))
dna_prot_temp=merge(dna_prot_core_temp,dna_prot_tails_temp,all=TRUE)
dna_prot_temp[is.na(dna_prot_temp)] <- 0
# dna_prot_max_sc=ddply(dna_prot_temp,c('DNA_resid','PROT_chain'),summarize,max_sc=max(tot_num_core+tot_num_tails),core=tot_num_core[which.max(tot_num_core+tot_num_tails)],tails=tot_num_tails[which.max(tot_num_core+tot_num_tails)])
dna_prot_avr=ddply(dna_prot_temp,c('DNA_resid','PROT_chain'),summarize,avr_sc=mean(tot_num_core)+mean(tot_num_tails),core=mean(tot_num_core),tails=mean(tot_num_tails))

dna_prot_avr$turn=c('positive')
dna_prot_avr[dna_prot_avr$DNA_resid<0,'turn']=c('negative')
dna_prot_avr[dna_prot_avr$DNA_resid<0,'DNA_resid']=-dna_prot_avr[dna_prot_avr$DNA_resid<0,'DNA_resid']
# t=subset(dna_prot_avr,DNA_resid==0)
# t$turn=c('negative')
# dna_prot_avr=rbind(dna_prot_avr,t)
head(dna_prot_avr)
dna_prot_avr$turn<-factor(dna_prot_avr$turn,levels=c('negative','positive'),labels=c('SHL < 0','SHL > 0'))

dna_prot_avr$PROT_chain=revalue(dna_prot_avr$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))
dna_prot_avr=ddply(dna_prot_avr,c('DNA_resid','PROT_chain','turn'),summarize,avr_sc=mean(avr_sc),core=mean(core),tails=mean(tails))
dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('DNA_resid','PROT_chain','turn'),measure.vars=c('core','tails'))
head(dna_prot_avr_m)

dna_prot_avr_m=arrange(dna_prot_avr_m,variable,PROT_chain)


dna_prot_avr_m$variable<-factor(dna_prot_avr_m$variable,levels=c('core','tails'),labels=c('CORE','TAILS'))


ann_text <- data.frame(x = c(8.3-0.15,0.8,4.85), y = c(36.5,36.5,36.5),lab = c('Linker DNA','Inner DNA turn','Outer DNA turn'),
                       variable = factor('CORE'),turn=factor('negative',levels=c('negative','positive'),labels=c('SHL < 0','SHL > 0')),PROT_chain='H3',
                       xmin=c(7.65-0.15,0,4.85-0.85),xmax=c(8.3+0.65-0.15,1.0+0.6,4.85+0.85),ymin=c(30,30,30),ymax=c(43,43,43))

theme_set(theme_grey(base_size = 18))


ann_num_inner=ddply(subset(dna_prot_avr_m,DNA_resid<39.5),c('PROT_chain','turn','variable'),summarize,avr_sum=round(sum(value),0))
t=ddply(subset(dna_prot_avr_m,DNA_resid<39.5),c('turn','variable'),summarize,avr_sum=round(sum(value),0))
t$PROT_chain='ALL'
ann_num_inner=rbind(ann_num_inner,t)
ann_num_inner$x=3.5
ann_num_inner$y=35
ann_num_inner[ann_num_inner$PROT_chain=='H3','x']=2.625-0.25
ann_num_inner[ann_num_inner$PROT_chain=='H3','y']=35
ann_num_inner[ann_num_inner$PROT_chain=='H4','x']=2.875-0.25
ann_num_inner[ann_num_inner$PROT_chain=='H4','y']=35
ann_num_inner[ann_num_inner$PROT_chain=='H2A','x']=3.125-0.25
ann_num_inner[ann_num_inner$PROT_chain=='H2A','y']=35
ann_num_inner[ann_num_inner$PROT_chain=='H2B','x']=3.375-0.25
ann_num_inner[ann_num_inner$PROT_chain=='H2B','y']=35
ann_num_inner[ann_num_inner$PROT_chain=='ALL','x']=3.375+0.25-0.25
ann_num_inner[ann_num_inner$PROT_chain=='ALL','y']=35


ann_num_inner[(ann_num_inner$PROT_chain=='H3')&(ann_num_inner$variable=='CORE')&(ann_num_inner$turn=='SHL > 0'),'y']=27
ann_num_inner[(ann_num_inner$PROT_chain=='H4')&(ann_num_inner$variable=='CORE')&(ann_num_inner$turn=='SHL > 0'),'y']=27
ann_num_inner[(ann_num_inner$PROT_chain=='H2A')&(ann_num_inner$variable=='CORE')&(ann_num_inner$turn=='SHL > 0'),'y']=27
ann_num_inner[(ann_num_inner$PROT_chain=='H2B')&(ann_num_inner$variable=='CORE')&(ann_num_inner$turn=='SHL > 0'),'y']=27
ann_num_inner[(ann_num_inner$PROT_chain=='ALL')&(ann_num_inner$variable=='CORE')&(ann_num_inner$turn=='SHL > 0'),'y']=27

ann_num_inner[(ann_num_inner$PROT_chain=='H3')&(ann_num_inner$variable=='TAILS')&(ann_num_inner$turn=='SHL < 0'),'y']=40
ann_num_inner[(ann_num_inner$PROT_chain=='H4')&(ann_num_inner$variable=='TAILS')&(ann_num_inner$turn=='SHL < 0'),'y']=40
ann_num_inner[(ann_num_inner$PROT_chain=='H2A')&(ann_num_inner$variable=='TAILS')&(ann_num_inner$turn=='SHL < 0'),'y']=40
ann_num_inner[(ann_num_inner$PROT_chain=='H2B')&(ann_num_inner$variable=='TAILS')&(ann_num_inner$turn=='SHL < 0'),'y']=40
ann_num_inner[(ann_num_inner$PROT_chain=='ALL')&(ann_num_inner$variable=='TAILS')&(ann_num_inner$turn=='SHL < 0'),'y']=40


ann_num_outer=ddply(subset(dna_prot_avr_m,(DNA_resid>39.5)&(DNA_resid<73.5)),c('PROT_chain','turn','variable'),summarize,avr_sum=round(sum(value),0))
t=ddply(subset(dna_prot_avr_m,(DNA_resid>39.5)&(DNA_resid<73.5)),c('turn','variable'),summarize,avr_sum=round(sum(value),0))
t$PROT_chain='ALL'
ann_num_outer=rbind(ann_num_outer,t)

ann_num_outer$x=7
ann_num_outer$y=35
ann_num_outer[ann_num_outer$PROT_chain=='H3','x']=6.25-0.125
ann_num_outer[ann_num_outer$PROT_chain=='H3','y']=35
ann_num_outer[ann_num_outer$PROT_chain=='H4','x']=6.375-0.25
ann_num_outer[ann_num_outer$PROT_chain=='H4','y']=35
ann_num_outer[ann_num_outer$PROT_chain=='H2A','x']=6.625-0.25
ann_num_outer[ann_num_outer$PROT_chain=='H2A','y']=35
ann_num_outer[ann_num_outer$PROT_chain=='H2B','x']=6.875-0.25
ann_num_outer[ann_num_outer$PROT_chain=='H2B','y']=35
ann_num_outer[ann_num_outer$PROT_chain=='ALL','x']=6.875+0.25-0.25
ann_num_outer[ann_num_outer$PROT_chain=='ALL','y']=35

ann_num_outer[(ann_num_outer$PROT_chain=='H3')&(ann_num_outer$variable=='CORE')&(ann_num_outer$turn=='SHL > 0'),'y']=27
ann_num_outer[(ann_num_outer$PROT_chain=='H4')&(ann_num_outer$variable=='CORE')&(ann_num_outer$turn=='SHL > 0'),'y']=27
ann_num_outer[(ann_num_outer$PROT_chain=='H2A')&(ann_num_outer$variable=='CORE')&(ann_num_outer$turn=='SHL > 0'),'y']=27
ann_num_outer[(ann_num_outer$PROT_chain=='H2B')&(ann_num_outer$variable=='CORE')&(ann_num_outer$turn=='SHL > 0'),'y']=27
ann_num_outer[(ann_num_outer$PROT_chain=='ALL')&(ann_num_outer$variable=='CORE')&(ann_num_outer$turn=='SHL > 0'),'y']=27

ann_num_outer[(ann_num_outer$PROT_chain=='H3')&(ann_num_outer$variable=='TAILS')&(ann_num_outer$turn=='SHL < 0'),'y']=40
ann_num_outer[(ann_num_outer$PROT_chain=='H4')&(ann_num_outer$variable=='TAILS')&(ann_num_outer$turn=='SHL < 0'),'y']=40
ann_num_outer[(ann_num_outer$PROT_chain=='H2A')&(ann_num_outer$variable=='TAILS')&(ann_num_outer$turn=='SHL < 0'),'y']=40
ann_num_outer[(ann_num_outer$PROT_chain=='H2B')&(ann_num_outer$variable=='TAILS')&(ann_num_outer$turn=='SHL < 0'),'y']=40
ann_num_outer[(ann_num_outer$PROT_chain=='ALL')&(ann_num_outer$variable=='TAILS')&(ann_num_outer$turn=='SHL < 0'),'y']=40


ann_sites <- data.frame(x = c(0.3,0.7,0.9,1.3,1.8,2.4,2.8,3.4,3.8,4.3,4.8,5.4,5.8,6.7), y = c(42),lab = c('L * 2','L * 1','alpha * N','alpha * 1','alpha * 1','L * 1','L * 2','L * 2','L * 1','alpha * 1','alpha * 1','L * 1','L * 2','alpha * N'),
                       variable = factor('CORE'),turn=factor('positive',levels=c('negative','positive'),labels=c('SHL < 0','SHL > 0')),PROT_chain=c('H3','H4','H3','H4','H3','H3','H4','H2B','H2A','H2A','H2B','H2B','H2A','H3'),
                       xmin=c(1.0),xmax=c(2.0),ymin=c(30),ymax=c(43))


q<-ggplot(data=subset(dna_prot_avr_m,value>0),aes(x=DNA_resid/10,y=value,fill=PROT_chain,color=PROT_chain)) + ggtitle("X-ray") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_bar(stat='identity',position='stack',width=0.05)+
facet_grid(variable+turn~.,scales='fixed')+ylab('Protein-DNA contacts')+
scale_x_continuous(limits=c(-0.05,9.3),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+scale_alpha_discrete(range = c(0.3, 1.0))+
scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(0,46),breaks = c(0,10,20,30,40),labels=c('','','20','','40'),name='Protein-DNA contacts',expand=c(0,0))+
geom_rect(data=ann_num_inner,aes(xmin=x-0.125,ymin=y-5,xmax=x+0.125,ymax=y+5,x=x,y=y),color='pink',fill='dark grey')+
geom_rect(data=ann_num_outer,aes(xmin=x-0.125,ymin=y-5,xmax=x+0.125,ymax=y+5,x=x,y=y),color='pink',fill='dark grey')+
# geom_rect(data=ann_num_linker,aes(xmin=x-0.125,ymin=y-5,xmax=x+0.125,ymax=y+5,x=x,y=y),color='pink',fill='dark grey')+

geom_text(data=ann_num_inner,aes(x=x,y=y,label=avr_sum))+
geom_text(data=ann_num_outer,aes(x=x,y=y,label=avr_sum))+
# geom_text(data=ann_num_linker,aes(x=x,y=y,label=avr_sum))+
geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=37,xmax=x+0.12,ymax=46,x=x,y=y),color='white',fill='grey')+
geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_avr_cryst.png",plot=q,width=15,height=5.5)



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


