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

dna_prot<-subset(dna_prot_raw,type%in%c('SC'))
dna_prot=subset(dna_prot,av_num>0.8)

# dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']


# dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))
dna_prot_core=dna_prot
#symmetrization
dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']=(-1)*dna_prot_core[dna_prot_core$DNA_resid<0,'DNA_resid']

dna_prot_core$PROT_chain=revalue(dna_prot_core$PROT_chain, c("CHA"="H3","CHE"="H3","CHB"="H4","CHF"="H4","CHC"="H2A","CHG"="H2A","CHD"="H2B","CHH"="H2B"))

#remove duplicates but leave it inclusive!!!! might do exclusive as well
dna_prot_core=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid','PROT_atom','DNA_atom','type','DNA_part','PROT_part'),summarize,av_num=length(av_num))
dna_prot_core[dna_prot_core$av_num<=1,'av_num']=1
dna_prot_core[dna_prot_core$av_num>1,'av_num']=2



# dna_prot_avr=ddply(dna_prot_core,c('DNA_resid','PROT_chain','PROT_resid'),summarize,core=max(av_num))
dna_prot_avr=ddply(dna_prot_core,c('PROT_chain','PROT_resid'),summarize,Chain=max(av_num))

dna_prot_avr[dna_prot_avr$Chain==1,'Chain']='Single'
dna_prot_avr[dna_prot_avr$Chain==2,'Chain']='Both'
dna_prot_avr_m=dna_prot_avr
# dna_prot_avr_m=melt(dna_prot_avr, id.vars=c('PROT_chain','PROT_resid'),measure.vars=c('core'))
head(dna_prot_avr_m)
# dna_prot_avr_m$DNA_part<-factor(dna_prot_avr_m$DNA_part,levels=c('phosphate','sugar','base'))

# type_lev=c('SB','HB','vdW','IM','WM','SC')

# dna_prot_avr_m$type<-factor(dna_prot_avr_m$type,levels=type_lev)

# dna_prot_avr_m$type=revalue(dna_prot_avr_m$type, c("SB"="1SB","HB"="2HB","vdW"="3vdW"))
# dna_prot_avr_m=arrange(dna_prot_avr_m,variable,type)


img <- readPNG(paste("seq_img/",'H3full_new',".png",sep=''))
h3 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4full_new',".png",sep=''))
h4 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Afull_new',".png",sep=''))
h2a <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Bfull_new',".png",sep=''))
h2b <- rasterGrob(img, interpolate=TRUE,width=1)

a<-ggplot(data=subset(dna_prot_avr_m,PROT_chain=='H3'),aes(x=PROT_resid,y=(-50),fill=Chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-101,58),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CHA','CHE'),labels=c('A','E'),name='Chain')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Super stable")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h3, ymin=0, ymax=53, xmin=0.5,xmax=135.5)
ggsave("../analysis_data/pub_int_prot_map_superstable_h3.png",plot=a,height=2,width=12)


b<-ggplot(data=subset(dna_prot_avr_m,PROT_chain=='H4'),aes(x=PROT_resid,y=(-50),fill=Chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-101,63),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,103),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CHA','CHE'),labels=c('A','E'),name='Chain')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Super stable")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h4, ymin=0, ymax=60, xmin=0.5,xmax=102.5)

ggsave("../analysis_data/pub_int_prot_map_superstable_h4.png",plot=b,height=2,width=10)



c<-ggplot(data=subset(dna_prot_avr_m,PROT_chain=='H2A'),aes(x=PROT_resid,y=(-50),fill=Chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-101,53),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,129),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CHA','CHE'),labels=c('A','E'),name='Chain')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Super stable")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2a, ymin=0, ymax=45, xmin=0.5,xmax=128.5)

ggsave("../analysis_data/pub_int_prot_map_superstable_h2a.png",plot=c,height=2,width=10)



d<-ggplot(data=subset(dna_prot_avr_m,PROT_chain=='H2B'),aes(x=PROT_resid,y=(-50),fill=Chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-101,60),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,123),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CHA','CHE'),labels=c('A','E'),name='Chain')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Super stable")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2b, ymin=0, ymax=50, xmin=0.5,xmax=122.5)

ggsave("../analysis_data/pub_int_prot_map_superstable_h2b.png",plot=d,height=2,width=10)



q<-arrangeGrob(a,b,c,d,ncol=1)


ggsave("../analysis_data/pub_int_prot_map_superstable.png",plot=q,height=8,width=12)



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