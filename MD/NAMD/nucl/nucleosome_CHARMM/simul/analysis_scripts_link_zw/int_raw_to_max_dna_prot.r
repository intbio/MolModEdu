#Here we max a frame of maximum interactions

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


#########Starting with DNA-protein interactions#########
###############
##Loading data frames
#DNA-protein interactions
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')

###Drop 250 nsec of equilibration
dna_prot=subset(dna_prot,Time>max(dna_prot$Time)/4)
print('kuku')

################
#Here we need to filter out HB and SB from contacts, and recalcualte WM interactions.
############

dna_prot=subset(dna_prot,type=='SC')
dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']=-dna_prot[dna_prot$DNA_chain=='CHJ','DNA_resid']

dna_prot_core=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 119))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29))))
dna_prot_tails=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid <= 43))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid <= 23))|((PROT_chain %in% c('CHC','CHG'))&((PROT_resid <= 15)|(PROT_resid >= 119)))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid <= 29))))


#we want to get maximum number of SC per nucleotide
dna_prot_core_temp=ddply(dna_prot_core,c('Time','DNA_resid','PROT_chain'),summarize,tot_num_core=length(param1))
dna_prot_tails_temp=ddply(dna_prot_tails,c('Time','DNA_resid','PROT_chain'),summarize,tot_num_tails=length(param1))
dna_prot_temp=merge(dna_prot_core_temp,dna_prot_tails_temp,all=TRUE)
dna_prot_temp[is.na(dna_prot_temp)] <- 0
# dna_prot_max_sc=ddply(dna_prot_temp,c('DNA_resid','PROT_chain'),summarize,max_sc=max(tot_num_core+tot_num_tails),core=tot_num_core[which.max(tot_num_core+tot_num_tails)],tails=tot_num_tails[which.max(tot_num_core+tot_num_tails)])
dna_prot_max_sc=ddply(dna_prot_temp,c('DNA_resid','PROT_chain'),summarize,max_sc=max(tot_num_core+tot_num_tails),core=max(tot_num_core),tails=max(tot_num_tails))

write.csv(dna_prot_max_sc, file="../analysis_data/dna_prot_max_sc_df.csv")
q()

#############
#adding important name columns, resnames, etc.

rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


sugar<-c('C1\'','C2\'','C3\'','C4\'','C5\'','O4\'')
phosphate<-c('O3\'','O5\'','O1P','O2P','P')
protbb<-c('C','O','N','CA','O','OT1','OT2')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SB','HB','vdW','IM','WM','SC')

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

if(!(df['type']%in%c('SC','vdW'))) { return('NA') }
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


