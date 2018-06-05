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

dna_prot<-subset(dna_prot,type=='WM')
dna_prot_cryst<-subset(dna_prot_cryst,type!='Z')
#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('IP','HB','VdW','IM','WM')

#Let's assign level order
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)




####General statistics section:
theme_set(theme_gray(base_size = 15))

##########Histograms with DNA_part classification


##------MD

# a2<-ggplot(data=dna_prot[dna_prot$PROT_resname!='LYS' & dna_prot$PROT_resname!='ARG',],aes(x=PROT_resname,alpha=PROT_part,fill=type,weight=av_num))+
# geom_bar(aes(color=type,y=..count..),position='stack')+
# scale_fill_manual(values=c( "blue", "purple", "grey", "dark green"))+
# scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c( "blue", "purple", "grey", "dark green"))+theme(legend.position="none")

a1<-ggplot(data=dna_prot,aes(x=PROT_resname,fill=PROT_part,weight=av_num))+
geom_bar(aes(y=..count..),position='stack')+
scale_fill_manual(values=c("blue", "green", "purple", "grey", "dark green"),name='Protein part')+scale_alpha_manual(values=c(1.0,0.5,0.1))+scale_color_manual(values=c("red", "blue", "purple", "grey", "dark green"))
a1<-a1+ggtitle('Water-mediated inteactions between DNA and protein in MD simulations')+xlab('Residue name')
# png(file="../analysis_data/int_dan_prot_cryst.png",width=2000,height=800)
# ggplot(data=dna_prot_cryst,aes(x=PROT_resname,fill=DNA_part,alpha=type))+geom_bar(aes(color=DNA_part,y=..count..),position='dodged')+scale_fill_manual(values=c("red", "blue", "green", "grey", "purple"))+scale_alpha_manual(values=c(1,0.5,0.2,0.1))+scale_color_manual(values=c("red", "blue", "green", "grey", "purple"))
# a1<-a1+annotation_custom(grob=ggplotGrob(a2),xmin=2.5,xmax=19,ymin=100,ymax=620)
ggsave(filename="../analysis_data/int_dna_prot_hist_protpart_anya_talk.png",plot=a1,width=10,height=5)

print('Total number of interactions')
print(sum(dna_prot$av_num))
