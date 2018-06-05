#This script will produce average for dna md data frame

library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
# library(gridExtra)
##############
INDEP=10 # frequency of independent points


errest<-function(x,numf){

return(sd(x)/sqrt(numf/INDEP))
}

myshift<-function(x){ #fits the vector to the left by 1, add 1 on end
return(c(x[-1],1))

}

dna<-read.csv('../../analysis_data/dna2_param_df_md.csv',na.strings=c("NA",'---'))
nf1=length(table(dna$Time))
dna$BPnum=dna$BPnum-20
dna2<-subset(dna,Time %in%seq(10,7500,10))
dna2$Time=as.integer(dna2$Time/10)
dna_dz=subset(dna2,BPnum%in%c(-39+74,39+74))[c('Time','BPnum','z')]
dna2_dz=dcast(dna_dz,Time~BPnum)
delta_z=dna2_dz['35'] - dna2_dz['113']
names(delta_z)='delta_z'
dna_prot<-read.csv('../../analysis_data/dna_prot_raw_df.csv')

#We need to get contacts
#between 24-28 'KKRK' H2B tails
#let's get contacts
dna_prot2=subset(dna_prot,type=='C'&PROT_chain%in%c('CHD','CHH')&PROT_resid%in%c(24,25,26,27,28))
dna_g1=subset(dna_prot2,DNA_resid<0)
dna_g2=subset(dna_prot2,DNA_resid>0)

dna_g1_cont=ddply(dna_g1,c('Time','PROT_chain'),summarize,num=length(param1))
dna_g2_cont=ddply(dna_g2,c('Time','PROT_chain'),summarize,num=length(param1))


dna_g1_c=dcast(dna_g1_cont,Time~PROT_chain)
dna_g2_c=dcast(dna_g2_cont,Time~PROT_chain)

sumdf1=merge(dna_g1_c,dna_g2_c,by=c('Time'),suffixes=c('gyre1','gyre2'))


findf=cbind(sumdf1,delta_z)

findf$all=findf$CHDgyre1+findf$CHDgyre2+findf$CHHgyre1+findf$CHHgyre2
findf$gyre1=findf$CHDgyre1+findf$CHHgyre1
findf$gyre2=findf$CHDgyre2+findf$CHHgyre2
findf$CHD=findf$CHDgyre1+findf$CHDgyre2
findf$CHH=findf$CHHgyre1+findf$CHHgyre2

write.csv(findf, file="contacts_and_deltaZ.csv")


fdf_m=melt(findf,id.vars='Time')
quit()


#let's do averaging and combination
#########
nf=length(table(dna_prot$Time))
#nf=750

####
cn=colnames(dna_prot)
avl=cn[!(cn=='Time' | cn=='param1' | cn=='param2' | cn=='param3')]

dna_prot=ddply(dna_prot,avl,summarize,av_num=length(param1)/nf)


#############
#adding important name columns, resnames, etc.

rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


sugar<-c('C1\'','C2\'','C3\'','C4\'','C5\'','O4\'')
phosphate<-c('O3\'','O5\'','O1P','O2P','P')
protbb<-c('C','O','N','CA')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('IP','HB','VdW','IM','WM','Z')

#We want to know for every interaction type get number of interactions between different types
# of amino acids and DNA groups
dp <- function(name){

if(name %in% sugar) { return('sugar') }
else if(name %in% phosphate){return('phosphate')}
else {return('base')}

}

pp <- function(name){

if(name %in% protbb) { return('backbone') }
else {return('side chain')}

}

#Let's add DNA parts classification
dna_prot_cryst$DNA_part<-as.factor(sapply(dna_prot_cryst$DNA_atom,dp))
dna_prot$DNA_part<-as.factor(sapply(dna_prot$DNA_atom,dp))

#Let's add prtoein parts classification
dna_prot_cryst$PROT_part<-as.factor(sapply(dna_prot_cryst$PROT_atom,pp))
dna_prot$PROT_part<-as.factor(sapply(dna_prot$PROT_atom,pp))


#Let's add DNA bases names and protein residies names
colnames(rr)<-c('DNA_chain','DNA_resid','DNA_resname')
dna_prot_cryst=merge(dna_prot_cryst,rr)
dna_prot=merge(dna_prot,rr)

colnames(rr)<-c('PROT_chain','PROT_resid','PROT_resname')
dna_prot_cryst=merge(dna_prot_cryst,rr)
dna_prot=merge(dna_prot,rr)

#Let's assign level order
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot$type<-factor(dna_prot$type,levels=type_lev)

dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))


dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)

#Now let's add contact classification
#PP, PN, NP

pol_nonpol <- function(df){

if(!(df['type']%in%c('C','VdW','Z'))) { return('NA') }
else {
if(substring(df['DNA_atom'],1,1) %in% c('O','P','N')){
if(substring(df['PROT_atom'],1,1) %in% c('C')) {return('PN')}
else{return('PP')}
}
else{
if(substring(df['PROT_atom'],1,1) %in% c('C')) {return('NN')}
else{return('PN')}

}
}
}

#Let's add pol-nonplo classification
dna_prot_cryst$C_type<-as.factor(apply(dna_prot_cryst[,c('DNA_atom','PROT_atom','type')],1,pol_nonpol))
dna_prot$C_type<-as.factor(apply(dna_prot[,c('DNA_atom','PROT_atom','type')],1,pol_nonpol))

#Let's add groove classification

# A: major: N6, C6, C5, N7, C8, minor:  C2,N3,C4
# G: major: O6, C6, C5, N7, C8, minor: N2, C2,N3,C4
# T: major: C6, C5, C4, C5M, O4 minor: N3, C2,O2
# C: major: C6, C5, C4, N4 minor: N3, C2,O2

m_m <- function(df){


if(df['DNA_resname']=='ADE') {
	if(df['DNA_atom'] %in% c('N6','C6','C5','N7','C8')){ return('major')}
	else if(df['DNA_atom'] %in% c('C2','N3','C4')){ return('minor')}
	else return('NA')}
if(df['DNA_resname']=='GUA') {
	if(df['DNA_atom'] %in% c('O6','C6','C5','N7','C8')){ return('major')}
	else if(df['DNA_atom'] %in% c('N2','C2','N3','C4')){ return('minor')}
	else return('NA')}
if(df['DNA_resname']=='THY') {
	if(df['DNA_atom'] %in% c('C6','C5','C4','C5M','O4')){ return('major')}
	else if(df['DNA_atom'] %in% c('N3','C2','O2')){ return('minor')}
	else return('NA')}

if(df['DNA_resname']=='CYT') {
	if(df['DNA_atom'] %in% c('C6','C5','C4','N4')){ return('major')}
	else if(df['DNA_atom'] %in% c('N3','C2','O2')){ return('minor')}
	else return('NA')}

}


dna_prot_cryst$groove<-as.factor(apply(dna_prot_cryst[,c('DNA_atom','DNA_resname')],1,m_m))
dna_prot$groove<-as.factor(apply(dna_prot[,c('DNA_atom','DNA_resname')],1,m_m))



write.csv(dna_prot, file="../analysis_data/dna_prot_avr_df.csv")
write.csv(dna_prot_cryst, file="../analysis_data/dna_prot_avr_df_cryst.csv")


