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

prot_prot_cryst=subset(prot_prot_cryst,type!='Z')
prot_prot=subset(prot_prot,type=='WM')

type_lev=c('IP','HB','VdW','IM','WM')


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

##########Histograms with PROT_part classification

#---Cryst
fh1=prot_prot_cryst[,c('PROT1_resname','PROT1_part','type')]
names(fh1)<-c('PROT_resname','PROT_part','type')
fh2=prot_prot_cryst[,c('PROT2_resname','PROT2_part','type')]
names(fh2)<-c('PROT_resname','PROT_part','type')
fh=rbind(fh1,fh2)

##------MD

fh1=prot_prot[,c('PROT1_resname','PROT1_part','type','av_num')]
names(fh1)<-c('PROT_resname','PROT_part','type','av_num')
fh2=prot_prot[,c('PROT2_resname','PROT2_part','type','av_num')]
names(fh2)<-c('PROT_resname','PROT_part','type','av_num')
fh=rbind(fh1,fh2)

a1<-ggplot(data=fh,aes(x=PROT_resname,fill=PROT_part,weight=av_num))+
geom_bar(aes(y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "green", "purple", "grey", "dark green"),name='Protein part')+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
a1<-a1+ggtitle('Water-mediated interactions between protein  chains in MD simulations')+xlab('Residue name')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=prot_prot_cryst,aes(x=PROT_resname,fill=PROT1_part,alpha=type))+geom_bar(aes(color=PROT1_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
ggsave(filename="../analysis_data/int_prot_prot_hist_protpart_anya_talk.png",plot=a1,width=10,height=5)

#Also we need to output statistics of WM interactions between histones

a=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH')
wm_prot_prot=data.frame()
for (i in a) {
for (j in a){

	print(i)
	print(j)
	s=sum(subset(prot_prot,(PROT1_chain == i)&(PROT2_chain == j))$av_num)
	print(s)
	wm_prot_prot=rbind(wm_prot_prot,data.frame(PROT1_chain=i,PROT2_chain=j,sum=s))
}

}

print(dcast(wm_prot_prot,PROT1_chain~PROT2_chain))

print('Total')
print(sum(wm_prot_prot$sum))