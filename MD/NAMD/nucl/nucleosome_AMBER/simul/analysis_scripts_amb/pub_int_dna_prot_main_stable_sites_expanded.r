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
dna_prot_raw<-read.csv('../analysis_data/dna_prot_avr_df.csv')



dna_prot<-subset(dna_prot_raw,type=='SC')
dna_prot=subset(dna_prot,av_num>0.8)

# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))

#symmetrization/splitting
dna_prot_core$turnI='SHL > 0'
dna_prot_core[dna_prot_core$DNA_resid<0,'turnI']='SHL < 0'

dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']=(-1)*dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']

dna_prot_core$PROT_chain=revalue(dna_prot_core$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))

#remove duplicates but leave it inclusive!!!! might do exclusive as well
dna_prot_core=ddply(dna_prot_core,c('turnI','DNA_resid','PROT_chain','PROT_resid','PROT_atom','DNA_atom','type','DNA_part','PROT_part'),summarize,av_num=length(av_num))
dna_prot_core$av_num=1

dna_prot_avr=ddply(dna_prot_core,c('turnI','DNA_resid','PROT_chain'),summarize,core=sum(av_num))
dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('turnI','DNA_resid','PROT_chain'),measure.vars=c('core'))
head(dna_prot_avr_m)
dna_prot_avr_m=arrange(dna_prot_avr_m,variable,PROT_chain)

# dna_prot_avr_m$variable<-factor(dna_prot_avr_m$variable,levels=c('core'),labels=c('CORE'))

ann_text <- data.frame(x = c(0.8,4.85), y = c(22.5,22.5),lab = c('Inner DNA turn','Outer DNA turn'),
                       variable = factor('CORE'),turnI=factor('negative',levels=c('negative','positive'),labels=c('SHL < 0','SHL > 0')),PROT_chain='H3',
                       xmin=c(0,4.85-0.85),xmax=c(1.0+0.6,4.85+0.85),ymin=c(20,20),ymax=c(25,25))

theme_set(theme_grey(base_size = 18))


ann_sites <- data.frame(x = c(0.3,0.7,0.9,1.3,1.8,2.4,2.8,3.4,3.8,4.3,4.8,5.4,5.8,6.7), y = c(19),lab = c('L * 2','L * 1','alpha * N','alpha * 1','alpha * 1','L * 1','L * 2','L * 2','L * 1','alpha * 1','alpha * 1','L * 1','L * 2','alpha * N'),
                       variable = factor('CORE'),turnI=factor('positive',levels=c('negative','positive'),labels=c('SHL < 0','SHL > 0')),PROT_chain=c('H3','H4','H3','H4','H3','H3','H4','H2B','H2A','H2A','H2B','H2B','H2A','H3'),
                       xmin=c(1.0),xmax=c(2.0),ymin=c(22),ymax=c(35))

#add dnd---


dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))

dnaseq2<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
dnaseq2<-gsub("([AG])", "R", dnaseq2)
dnaseq2<-gsub("([CT])", "Y", dnaseq2)

seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))
seqdf2=data.frame(sequence=substring(dnaseq2, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))


#dfcryst$X<-as.factor(dfcryst$X)
df_dna<-read.csv("../analysis_data/dna_rot_df_avr.csv",header=TRUE,check.name=FALSE)
# df_cryst<-read.csv("../analysis_data/dna_rot_df_cryst.csv",header=TRUE,check.name=FALSE)

# df_avr$DATA='Average'
# head(df_avr)
# df_cryst$DATA='X-ray'
# head(df_cryst)
# df=rbind(df_avr,df_cryst)
dfd=df_dna

dfp=subset(dfd,Basepair>=0)
dfm=subset(dfd,Basepair<=0)
dfm$Basepair=(-1)*dfm$Basepair
dfp$turn='positive'
dfm$turn='negative'
dfd=rbind(dfp,dfm)

dfd[dfd$Angle>0,'Angle']=dfd[dfd$Angle>0,'Angle']-180

q<-ggplot(data=subset(dna_prot_avr_m,value>0),aes(x=DNA_resid/10,y=value)) + ggtitle("MD simulations, NCPamb-model, average, interactions with > 80% stability") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_bar(aes(fill=PROT_chain,color=PROT_chain),stat='identity',position='stack',width=0.05)+
# facet_grid(turn~.,scales='fixed')+
ylab('Stable protein-DNA contacts: total number')+
scale_x_continuous(limits=c(-0.05,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+
scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(0,25),breaks = c(0,5,10,15),labels=c('0','5','10','15'),name='Contacts or angle/18',expand=c(0,0))+

geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=16,xmax=x+0.12,ymax=21,x=x,y=y),color='white',fill='grey')+
geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)+
geom_line(data=dfd,aes(x=Basepair/10,y=-Angle/18,alpha=turn),fill='black',color='black',size=1)+
geom_point(data=dfd,aes(x=Basepair/10,y=-Angle/18,alpha=turn),fill='black',color='black',size=2)+
scale_alpha_discrete(range = c(0.3, 0.6),name='SHL')+
facet_grid(turnI~.)
# scale_y_continuous(limits=c(0,25),breaks = c(0,5,10,15),labels=c('0','90','180','270'),name='Angle, deg',expand=c(0,0))


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_main_sites_expanded.png",plot=q,width=15,height=3.5)
q()

#Draw int types####################################

dna_prot<-subset(dna_prot_raw,type%in%c('SB','HB','vdW'))
dna_prot=subset(dna_prot,av_num>0.8)

# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))

#symmetrization
dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']=(-1)*dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']

dna_prot_core$PROT_chain=revalue(dna_prot_core$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))

#remove duplicates but leave it inclusive!!!! might do exclusive as well
dna_prot_core=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid','PROT_atom','DNA_atom','type','DNA_part','PROT_part'),summarize,av_num=length(av_num))
dna_prot_core$av_num=1





dna_prot_avr=ddply(dna_prot_core,c('DNA_resid','type'),summarize,core=sum(av_num))
dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('DNA_resid','type'),measure.vars=c('core'))
head(dna_prot_avr_m)
type_lev=c('SB','HB','vdW','IM','WM','SC')

dna_prot_avr_m$type<-factor(dna_prot_avr_m$type,levels=type_lev)

# dna_prot_avr_m$type=revalue(dna_prot_avr_m$type, c("SB"="1SB","HB"="2HB","vdW"="3vdW"))
# dna_prot_avr_m=arrange(dna_prot_avr_m,variable,type)



q<-ggplot(data=subset(dna_prot_avr_m,value>0),aes(x=DNA_resid/10,y=value,fill=type,color=type,order=type)) + ggtitle("MD simulations, FN-system, average, interactions with > 80% stability") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_bar(stat='identity',position='stack',width=0.05)+
# facet_grid(turn~.,scales='fixed')+
ylab('Stable protein-DNA contacts: total number')+
scale_x_continuous(limits=c(-0.05,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+scale_alpha_discrete(range = c(0.3, 1.0))+
# scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
# scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
# geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(0,15),breaks = c(0,5,10,15),labels=c('0','5','10','15'),name='Stable contacts, total',expand=c(0,0))+

# geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=14,xmax=x+0.12,ymax=19,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),name='Type intr')+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),guide=FALSE)


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_main_sites_type.png",plot=q,width=15,height=2.5)

################DNA part

dna_prot<-subset(dna_prot_raw,type%in%c('SC'))
dna_prot=subset(dna_prot,av_num>0.8)

# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))

#symmetrization
dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']=(-1)*dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']

dna_prot_core$PROT_chain=revalue(dna_prot_core$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))

#remove duplicates but leave it inclusive!!!! might do exclusive as well
dna_prot_core=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid','PROT_atom','DNA_atom','type','DNA_part','PROT_part'),summarize,av_num=length(av_num))
dna_prot_core$av_num=1


dna_prot_avr=ddply(dna_prot_core,c('DNA_resid','DNA_part'),summarize,core=sum(av_num))
dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('DNA_resid','DNA_part'),measure.vars=c('core'))
head(dna_prot_avr_m)
dna_prot_avr_m$DNA_part<-factor(dna_prot_avr_m$DNA_part,levels=c('phosphate','sugar','base'))

# type_lev=c('SB','HB','vdW','IM','WM','SC')

# dna_prot_avr_m$type<-factor(dna_prot_avr_m$type,levels=type_lev)

# dna_prot_avr_m$type=revalue(dna_prot_avr_m$type, c("SB"="1SB","HB"="2HB","vdW"="3vdW"))
# dna_prot_avr_m=arrange(dna_prot_avr_m,variable,type)



d<-ggplot(data=subset(dna_prot_avr_m,value>0),aes(x=DNA_resid/10,y=value,fill=DNA_part,color=DNA_part)) +# ggtitle("MD simulations, FN-system, average, interactions with > 80% stability,DNA part") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_bar(stat='identity',position='stack',width=0.05)+
# facet_grid(turn~.,scales='fixed')+
ylab('Stable protein-DNA contacts: total number')+
scale_x_continuous(limits=c(-0.05,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+scale_alpha_discrete(range = c(0.3, 1.0))+
# scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
# scale_fill_manual(name='DNA part')+

# scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
# geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(0,15),breaks = c(0,5,10,15),labels=c('0','5','10','15'),name='Stable contacts, total',expand=c(0,0))

# geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=14,xmax=x+0.12,ymax=19,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)+
# scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),name='Type intr')+
# scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),guide=FALSE)


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_main_sites_dna_part.png",plot=d,width=17,height=2.5)



####Protein part


dna_prot<-subset(dna_prot_raw,type%in%c('SC'))
dna_prot=subset(dna_prot,av_num>0.8)

# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))

#symmetrization
dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']=(-1)*dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']

dna_prot_core$PROT_chain=revalue(dna_prot_core$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))

#remove duplicates but leave it inclusive!!!! might do exclusive as well
dna_prot_core=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid','PROT_atom','DNA_atom','type','DNA_part','PROT_part'),summarize,av_num=length(av_num))
dna_prot_core$av_num=1


dna_prot_avr=ddply(dna_prot_core,c('DNA_resid','PROT_part'),summarize,core=sum(av_num))
dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('DNA_resid','PROT_part'),measure.vars=c('core'))
head(dna_prot_avr_m)
# dna_prot_avr_m$DNA_part<-factor(dna_prot_avr_m$DNA_part,levels=c('phosphate','sugar','base'))

# type_lev=c('SB','HB','vdW','IM','WM','SC')

# dna_prot_avr_m$type<-factor(dna_prot_avr_m$type,levels=type_lev)

# dna_prot_avr_m$type=revalue(dna_prot_avr_m$type, c("SB"="1SB","HB"="2HB","vdW"="3vdW"))
# dna_prot_avr_m=arrange(dna_prot_avr_m,variable,type)



d<-ggplot(data=subset(dna_prot_avr_m,value>0),aes(x=DNA_resid/10,y=value,fill=PROT_part,color=PROT_part)) + #ggtitle("MD simulations, FN-system, average, interactions with > 80% stability,DNA part") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
geom_bar(stat='identity',position='stack',width=0.05)+
# facet_grid(turn~.,scales='fixed')+
ylab('Stable protein-DNA contacts: total number')+
scale_x_continuous(limits=c(-0.05,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+scale_alpha_discrete(range = c(0.3, 1.0))+
scale_fill_manual(breaks=c('backbone','side chain'),values=c('backbone'='magenta','side chain'='dark green'),labels=c('backbone','sidechain'),name='Protein part')+
# scale_fill_manual(name='DNA part')+
scale_color_manual(breaks=c('backbone','side chain'),values=c('backbone'='magenta','side chain'='dark green'),labels=c('backbone','sidechain'),name='Protein part',guide=FALSE)+

# scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
# geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(0,15),breaks = c(0,5,10,15),labels=c('0','5','10','15'),name='Stable contacts, total',expand=c(0,0))

# geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=14,xmax=x+0.12,ymax=19,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)+
# scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),name='Type intr')+
# scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),guide=FALSE)


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_main_sites_prot_part.png",plot=d,width=17,height=2.5)


###Residue numbers

dna_prot<-subset(dna_prot_raw,type%in%c('SC'))
dna_prot=subset(dna_prot,av_num>0.8)

# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))

#symmetrization
dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']=(-1)*dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']

dna_prot_core$PROT_chain=revalue(dna_prot_core$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))

#remove duplicates but leave it inclusive!!!! might do exclusive as well
dna_prot_core=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid','PROT_atom','DNA_atom','type','DNA_part','PROT_part'),summarize,av_num=length(av_num))
dna_prot_core$av_num=1


dna_prot_avr=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid'),summarize,core=sum(av_num))
dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('DNA_resid','PROT_chain','PROT_resid'),measure.vars=c('core'))
head(dna_prot_avr_m)
# dna_prot_avr_m$DNA_part<-factor(dna_prot_avr_m$DNA_part,levels=c('phosphate','sugar','base'))

# type_lev=c('SB','HB','vdW','IM','WM','SC')

# dna_prot_avr_m$type<-factor(dna_prot_avr_m$type,levels=type_lev)

# dna_prot_avr_m$type=revalue(dna_prot_avr_m$type, c("SB"="1SB","HB"="2HB","vdW"="3vdW"))
# dna_prot_avr_m=arrange(dna_prot_avr_m,variable,type)



d<-ggplot(data=subset(dna_prot_avr_m,value>0),aes(x=DNA_resid/10,y=2.2,color=PROT_chain,label=PROT_resid)) + ggtitle("MD simulations, FN-system, average, interactions with > 80% stability,Resnumbers") + 
xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
# geom_bar(stat='identity',position='stack',width=0.05)+
# facet_grid(turn~.,scales='fixed')+
ylab('Stable protein-DNA contacts: total number')+
scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
scale_x_continuous(limits=c(-0.05,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+scale_alpha_discrete(range = c(0.3, 1.0))+
# scale_fill_manual(breaks=c('backbone','side chain'),values=c('backbone'='magenta','side chain'='dark green'),labels=c('backbone','sidechain'),name='Protein part')+
# scale_fill_manual(name='DNA part')+
# scale_color_manual(breaks=c('backbone','side chain'),values=c('backbone'='magenta','side chain'='dark green'),labels=c('backbone','sidechain'),name='Protein part',guide=FALSE)+
geom_text(position = 'stack',size=5)+

# scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
# geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(0,15),breaks = c(0,5,10,15),labels=c('0','5','10','15'),name='Stable contacts, total',expand=c(0,0))

# geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=14,xmax=x+0.12,ymax=19,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)+
# scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),name='Type intr')+
# scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),guide=FALSE)


# q<-arrangeGrob(c,cmd,ncol=1)
ggsave(filename="../analysis_data/pub_int_dna_prot_main_sites_resnum.png",plot=d,width=16,height=2.5)



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


# 