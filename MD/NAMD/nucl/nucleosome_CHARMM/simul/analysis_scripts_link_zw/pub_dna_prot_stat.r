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
##Loading data frames
#DNA-protein interactions
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
dna_prot_cryst<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_avr_df.csv')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
prot_rn_lev=c('ARG','LYS','THR','SER','VAL','GLY','ILE','PHE','GLN','ALA','LEU','PRO','GLU','TYR','MET','HSE','ASP','CYS','ASN')
prot_rn_levn=c('ARG','LYS','THR','SER','VAL','GLY','ILE','PHE','GLN','ALA','LEU','PRO','GLU','TYR','MET','HIS','ASP','CYS','ASN')

type_lev=c('SC','SB','HB','vdW','IM','WM')
# head(subset(dna_prot_cryst,type=='SB'))
#Let's assign level order
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)


dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))
dna_prot_tails=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid <= 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid <= 23))|((PROT_chain %in% c('CHC','CHG'))&((PROT_resid <= 15)|(PROT_resid >= 119)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid <= 29))))
dna_prot_core$hist_part='Core'
dna_prot_tails$hist_part='Tails'
dna_prot=rbind(dna_prot_core,dna_prot_tails)

dna_prot_core_cryst=subset(dna_prot_cryst, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))
dna_prot_tails_cryst=subset(dna_prot_cryst, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid <= 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid <= 23))|((PROT_chain %in% c('CHC','CHG'))&((PROT_resid <= 15)|(PROT_resid >= 119)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid <= 29))))
dna_prot_core_cryst$hist_part='Core'
dna_prot_tails_cryst$hist_part='Tails'
dna_prot_cryst=rbind(dna_prot_core_cryst,dna_prot_tails_cryst)

dna_prot$DATA='MD'
dna_prot_cryst$DATA='X-ray'
cn=colnames(dna_prot_cryst)
avl=cn[!(cn=='param1' | cn=='param2' | cn=='param3')]
dna_prot_cryst=dna_prot_cryst[,avl]
dna_prot_cryst$av_num=1

dna_prot_comb=rbind(dna_prot_cryst,dna_prot)

dna_prot_comb=subset(dna_prot_comb,type%in%c('SB','vdW','HB'))
####General statistics section:
theme_set(theme_bw(base_size = 15)+theme(panel.grid.minor=element_blank()))

dna_prot_comb_sum=ddply(dna_prot_comb,c('hist_part','PROT_resname','DATA','type'),summarize,tot_num=sum(av_num))

# dna_prot_comb_sum$PROT_resname=revalue(dna_prot_comb_sum$PROT_resname, 
	# c('ARG'=1,'LYS'=2,'THR'=3,'SER'=4,'VAL'=5,'GLY'=6,'ILE'=7,'PHE'=8,'GLN'=9,'TYR'=10,'HSE'=11,'ASN'=12,'ALA'=13,'LEU'=14,'PRO'=15,'GLU'=16,'ASP'=17,'MET'=18,'CYS'=19))
# head(dna_prot_comb_sum,n=30)
dna_prot_comb_sum$PROT_resname=as.numeric(dna_prot_comb_sum$PROT_resname)
# head(dna_prot_comb_sum,n=30)

# dna_prot_comb_sum=arrange(dna_prot_comb_sum,PROT_resname)


a<-ggplot(data=subset(dna_prot_comb_sum,hist_part=='Core'),aes(x=PROT_resname-0.2,alpha=hist_part,fill=type,y=tot_num))+
geom_bar(aes(color=type),stat='identity',width=0.4,position='stack')+
geom_bar(data=subset(dna_prot_comb_sum,hist_part=='Tails'),aes(x=PROT_resname+0.2,color=type),stat='identity',width=0.4,position='stack')+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),label=c('SB'='Salt bridges','HB'='H-bonds','vdW'='vdW contacts'),name='Type')+
scale_alpha_manual(values=c(1.0,0.3,0.1),name='Region')+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"),guide=FALSE)+
ggtitle('Interactions between DNA and protein')+
scale_x_continuous(limits=c(0.5,17.5),breaks = seq(1, 19),labels=prot_rn_levn,expand=c(0,0),name='Amino acid residues')+
facet_grid(DATA~.,scales = "free_x")+ylab('Number of interactions')


ggsave(filename="../analysis_data/pub_dna_prot_stat.png",plot=a,width=10,height=5)



a<-ggplot(data=subset(dna_prot_comb_sum,(hist_part=='Core') & !(PROT_resname%in%c(1,2))),aes(x=PROT_resname-0.2,alpha=hist_part,fill=type,y=tot_num))+
geom_bar(aes(color=type),stat='identity',width=0.4,position='stack')+
geom_bar(data=subset(dna_prot_comb_sum,(hist_part=='Tails') & !(PROT_resname%in%c(1,2))),aes(x=PROT_resname+0.2,color=type),stat='identity',width=0.4,position='stack')+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.3,0.1),name='')+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
ggtitle('Interactions between DNA and protein')+
scale_x_continuous(limits=c(2.5,17.5),breaks = seq(1, 19),labels=prot_rn_levn,expand=c(0,0),name='Amino acid residues')+
facet_grid(DATA~.,scales = "free_x")+ylab('Number of interactions')



ggsave(filename="../analysis_data/pub_dna_prot_stat_zoom.png",plot=a,width=9,height=5)



# q()

#Now draw table
theme_set(theme_bw(base_size = 15))


dna_prot_comb=rbind(dna_prot_cryst,dna_prot)

dna_prot_comb=subset(dna_prot_comb,type%in%c('SB','vdW','HB'))

ann_num=ddply(dna_prot_comb,c('type','DATA','PROT_part','DNA_part','hist_part'),summarize,avr_sum=round(sum(av_num),0))
ann_num$DNA_part=revalue(ann_num$DNA_part,c('phosphate'='phosph')) 


q<-ggplot(data=ann_num,aes(x=as.numeric(type),y=0,label=avr_sum)) + #ggtitle("MD simulations, FN-system, average") + 
# xlab("Base pair")+#geom_line(aes(x=DNA_resid,y=number,color=PROT_chain),size=1)+
# geom_bar(stat='identity',position='stack',width=0.05)+
# facet_grid(variable+turn~.,scales='fixed')+ylab('Protein-DNA contacts: average number')+
scale_x_continuous(limits=c(1.5,4.5),breaks = c(2.5,3.5),labels=c(),expand=c(0,0),name='')+
# scale_alpha_discrete(range = c(0.3, 1.0))+
# scale_fill_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+
# scale_color_manual(breaks=c('H3','H4','H2A','H2B'),values=c('H3'='blue','H4'='green','H2A'='gold','H2B'='red','ALL'='black'),labels=c('H3','H4','H2A','H2B'),name='Histones')+#+#+xlim(-70,70)
# geom_vline(xintercept = 7.35)+geom_vline(xintercept = 3.95)+
# geom_rect(data=ann_text,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_text,aes(x=x,y=y,label=lab),color='black',fill='black',size=8)+
scale_y_continuous(limits=c(-0.6,0.6),breaks = c(),labels=c(),name='',expand=c(0,0))+
# geom_rect(data=ann_num_inner,aes(xmin=x-0.125,ymin=y-5,xmax=x+0.125,ymax=y+5,x=x,y=y),color='pink',fill='dark grey')+
# geom_rect(data=ann_num_outer,aes(xmin=x-0.125,ymin=y-5,xmax=x+0.125,ymax=y+5,x=x,y=y),color='pink',fill='dark grey')+
# geom_rect(data=ann_num_linker,aes(xmin=x-0.125,ymin=y-5,xmax=x+0.125,ymax=y+5,x=x,y=y),color='pink',fill='dark grey')+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+

# geom_text(data=ann_num_inner,aes(x=x,y=y,label=avr_sum))+
# geom_text(data=ann_num_outer,aes(x=x,y=y,label=avr_sum))+
# geom_text(data=ann_num_linker,aes(x=x,y=y,label=avr_sum))+
# geom_rect(data=ann_sites,aes(xmin=x-0.12,ymin=37,xmax=x+0.12,ymax=46,x=x,y=y),color='white',fill='grey')+
# geom_text(data=ann_sites,aes(x=x,y=y,label=lab,color=PROT_chain),fill='black',size=6,parse=TRUE)
geom_text(aes(color=type),size=6)+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1),name='')+
facet_grid(DATA+DNA_part~hist_part+PROT_part,scales = "fixed")+
theme( panel.grid.minor = element_blank())


ggsave(filename="../analysis_data/pub_dna_prot_table.png",plot=q,width=10,height=5)






dna_prot_comb=rbind(dna_prot_cryst,dna_prot)

dna_prot_comb=subset(dna_prot_comb,type%in%c('SC'))
ann_num=ddply(dna_prot_comb,c('type','DATA','PROT_part','DNA_part','hist_part'),summarize,avr_sum=round(sum(av_num),0))
ann_num$DNA_part=revalue(ann_num$DNA_part,c('phosphate'='phosph')) 
head(ann_num)


q<-ggplot(data=ann_num,aes(x=0,y=0,label=avr_sum)) + #ggtitle("MD simulations, FN-system, average") + 
scale_x_continuous(limits=c(-0.6,0.6),breaks = c(),labels=c(),name='',expand=c(0,0))+
scale_y_continuous(limits=c(-0.6,0.6),breaks = c(),labels=c(),name='',expand=c(0,0))+
scale_color_manual(values=c('SC'="black"))+

geom_text(aes(color=type),size=6)+
# scale_fill_manual(values=c('SC'="black"))+
scale_alpha_manual(values=c(1.0,0.5,0.1),name='')+
facet_grid(DATA+DNA_part~hist_part+PROT_part,scales = "fixed")+
theme( panel.grid.minor = element_blank())


ggsave(filename="../analysis_data/pub_dna_prot_table_SC.png",plot=q,width=10,height=5)




q()


##########Histograms with DNA_part classification

#---Cryst
a2<-ggplot(data=dna_prot_cryst[dna_prot_cryst$PROT_resname!='LYS' & dna_prot_cryst$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=DNA_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_y_continuous(limits=c(0,120),minor_breaks = seq(0, 600, 5),breaks = seq(0, 600, 10))+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+theme(legend.position="none")


a1<-ggplot(data=dna_prot_cryst,aes(x=PROT_resname,alpha=DNA_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_y_continuous(limits=c(0,590),minor_breaks = seq(0, 600, 25),breaks = seq(0, 600, 50))

a1<-a1+ggtitle('Interactions between DNA and protein in crystal')
c1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=120,ymax=620)
# ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart_cryst.png",plot=a1,width=10,height=5)

##------MD

a2<-ggplot(data=dna_prot[dna_prot$PROT_resname!='LYS' & dna_prot$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=DNA_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_y_continuous(limits=c(0,120),minor_breaks = seq(0, 600, 5),breaks = seq(0, 600, 10))+
scale_fill_manual(values=c( 'SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot,aes(x=PROT_resname,alpha=DNA_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in MD simulations (average count)')+
scale_y_continuous(limits=c(0,590),minor_breaks = seq(0, 600, 25),breaks = seq(0, 600, 50))

# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
m1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=120,ymax=620)
# ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart.png",plot=a1,width=10,height=5)

q<-arrangeGrob(c1,m1,ncol=1)
ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart_mdcryst.png",plot=q,width=10,height=10)

###Histograms with PROT_part classification


#---Cryst
a2<-ggplot(data=dna_prot_cryst[dna_prot_cryst$PROT_resname!='LYS' & dna_prot_cryst$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_y_continuous(limits=c(0,120),minor_breaks = seq(0, 600, 5),breaks = seq(0, 600, 10))+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+theme(legend.position="none")


a1<-ggplot(data=dna_prot_cryst,aes(x=PROT_resname,alpha=PROT_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_y_continuous(limits=c(0,590),minor_breaks = seq(0, 600, 25),breaks = seq(0, 600, 50))

a1<-a1+ggtitle('Interactions between DNA and protein in crystal')
c1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=120,ymax=620)
# ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart_cryst.png",plot=a1,width=10,height=5)

##------MD

a2<-ggplot(data=dna_prot[dna_prot$PROT_resname!='LYS' & dna_prot$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_y_continuous(limits=c(0,120),minor_breaks = seq(0, 600, 5),breaks = seq(0, 600, 10))+
scale_fill_manual(values=c( 'SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1))+
scale_color_manual(values=c( 'SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot,aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
a1<-a1+ggtitle('Interactions between DNA and protein in MD simulations (average count)')+
scale_y_continuous(limits=c(0,590),minor_breaks = seq(0, 600, 25),breaks = seq(0, 600, 50))

# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
m1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=120,ymax=620)
# ggsave(filename="../analysis_data/int_dna_prot_hist_dnapart.png",plot=a1,width=10,height=5)

q<-arrangeGrob(c1,m1,ncol=1)
ggsave(filename="../analysis_data/int_dna_prot_hist_protpart_mdcryst.png",plot=q,width=10,height=10)

####Histogram of contacts with bases
theme_set(theme_gray(base_size = 12))

#---Cryst
c<-ggplot(data=dna_prot_cryst[dna_prot_cryst$DNA_part=='base',],aes(x=PROT_resname,alpha=PROT_part,fill=type))+
geom_bar(aes(color=type,y=..count..),position='stack')+ylim(0,65)+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c('backbone'=1.0,'side chain'=0.5))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))
c<-c+ggtitle('Interactions between DNA bases and protein in crystal')

##------MD

md<-ggplot(data=dna_prot[dna_prot$DNA_part=='base',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
geom_bar(aes(color=type,y=..count..),position='stack')+ylim(0,65)+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c('backbone'=1.0,'side chain'=0.5),breaks=c('backbone','side chain'))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))
md<-md+ggtitle('Interactions between DNA bases and protein in MD simulations (average count)')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
q<-arrangeGrob(c,md)
ggsave(filename="../analysis_data/int_dna_base_prot_hist.png",plot=q,width=10,height=5)




###Histograms of contact types
#---Cryst

a1<-ggplot(data=dna_prot_cryst[dna_prot_cryst$type=='vdW',],aes(x=PROT_resname,fill=C_type))+
geom_bar(aes(color=C_type,y=..count..),position='stack')+
scale_y_continuous(limits=c(0,470),minor_breaks = seq(0, 600, 25),breaks = seq(0, 600, 50))+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1),breaks=c('backbone','side chain'))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))
c1<-a1+ggtitle('Interactions between DNA and protein in crystal: contact types')
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
# ggsave(filename="../analysis_data/int_dna_prot_hist_ctype_cryst.png",plot=a1,width=10,height=5)

##------MD

a1<-ggplot(data=dna_prot[dna_prot$type=='vdW',],aes(x=PROT_resname,fill=C_type,weight=av_num))+
geom_bar(aes(color=C_type,y=..count..),position='stack')+
scale_y_continuous(limits=c(0,470),minor_breaks = seq(0, 600, 25),breaks = seq(0, 600, 50))+
scale_fill_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))+
scale_alpha_manual(values=c(1.0,0.5,0.1),breaks=c('backbone','side chain'))+
scale_color_manual(values=c('SB'="red",'HB'="blue",'vdW'="purple", 'IM'='grey', 'WM'="dark green"))
m1<-a1+ggtitle('Interactions between DNA and protein in MD: contact types')
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=3,xmax=18,ymin=110,ymax=800)
# ggsave(filename="../analysis_data/int_dna_prot_hist_ctype.png",plot=a1,width=10,height=5)
q<-arrangeGrob(c1,m1)
ggsave(filename="../analysis_data/int_dna_prot_hist_ctype_mdcryst.png",plot=q,width=10,height=10)

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
d1j$DNA_resid<-q[d1j$DNA_resid+74]
#Let's add zeros also

d1zt=data.frame(DNA_resid=seq(73,-73,-1),DNA_chain=c('CHI'),num=c(0))
t=data.frame(type=c('vdW','WM','SB','IM','HB'))
d1z=merge(d1zt,t)

d1<-rbind(d1i,d1j,d1z)

tot_cryst=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))

c<-ggplot(data=tot_cryst[tot_cryst$type %in% c('vdW','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
o<-ggplot(data=tot_cryst[tot_cryst$type %in% c('SB','HB'),]) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple", "grey","dark green"))

# q<-arrangeGrob(c,o)
# ggsave(filename="../analysis_data/int_dna_prot_prof.png",plot=q,width=15,height=5)

##MD----------
d1=ddply(dna_prot,c("DNA_chain","DNA_resid","type"),function(df) c(num=sum(df$av_num)))
d1i=d1[d1$DNA_chain=='CHI',]
d1j=d1[d1$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+74]
d1<-rbind(d1i,d1j,d1z)
tot_md=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))


cmd<-ggplot(data=tot_md[tot_md$type %in% c('vdW','WM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
omd<-ggplot(data=tot_md[tot_md$type %in% c('SB','HB','IM'),]) + ggtitle("DNA-protein interactions in nucleosome: MD") + 
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
chi_dna_part_cwm<-ggplot(data=int_md[(int_md$type %in% c('vdW','WM')) & (int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, vdW & WM") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))+
facet_grid(DNA_part~.,scales='free')

chi_dna_part_iphbim<-ggplot(data=int_md[(int_md$type %in% c('SB','HB','IM')) & (int_md$DNA_chain=='CHI'),]) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, SB & HB & IM") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("red", "purple",'grey'))+
facet_grid(DNA_part~.,scales='free')

ggsave(filename="../analysis_data/int_dna_prot_dp_chi_cwm.png",plot=chi_dna_part_cwm,width=15,height=10)
ggsave(filename="../analysis_data/int_dna_prot_dp_chi_iphbim.png",plot=chi_dna_part_iphbim,width=15,height=10)

#By contact type

d1=ddply(subset(dna_prot,type=='vdW'),c("DNA_chain","DNA_resid","C_type","DNA_part"),function(df) c(num=sum(df$av_num)))
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
d1=ddply(subset(dna_prot,PROT_resname=='ARG' & type %in% c('SB','HB','vdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with ARG (vdW+SB+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_arg.png",plot=chi_dna_part_arg_int,width=15,height=10)

#Crystal
d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG' & type %in% c('SB','HB','vdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_arg_int_cryst<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, Crystal, interactions with ARG (vdW+SB+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_arg_cryst.png",plot=chi_dna_part_arg_int_cryst,width=15,height=10)


####Let's specifically study LYSinines


d1=ddply(subset(dna_prot,PROT_resname=='LYS' & type %in% c('SB','HB','vdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
#Add zero frames
d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHI'))
d2zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0),DNA_chain=c('CHJ'))
t=data.frame(DNA_part=c('base','sugar','phosphate'))
dzt<-rbind(d1zt,d2zt)
dz=merge(dzt,t)
d1=rbind(d1,dz)

int_md=ddply(d1,c('DNA_chain',"DNA_resid",'DNA_part'),summarize,number=sum(num))

chi_dna_part_lys_int<-ggplot(data=int_md) + ggtitle("DNA-protein interactions in nucleosome, MD, interactions with LYS (vdW+SB+HB)") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=DNA_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('red',"blue", "purple"))+
facet_grid(DNA_part~.,scales='free')
ggsave(filename="../analysis_data/int_dna_prot_dp_lys.png",plot=chi_dna_part_lys_int,width=15,height=10)



####Misc


####Let's specifically study ARGinines

#MD
d1=ddply(subset(dna_prot,PROT_resname=='ARG' & type %in% c('vdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=sum(df$av_num)))
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
d1=ddply(subset(dna_prot_cryst,PROT_resname=='ARG' & type %in% c('vdW') ),c('DNA_chain',"DNA_resid","DNA_part"),function(df) c(num=length(df$param1)))
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


