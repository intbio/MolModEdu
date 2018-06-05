#R-script  conversts raw data frames to avarage data frames
# makes important filtering - eg. deletes duplicates of HB and SB from contacts
# and adds many useful columns for analysis:
# residue names, atomic groups etc.
# introdues some order to factors

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
dna_prot_cryst<-read.csv('../analysis_data/dna_prot_raw_df_cryst.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')

###Drop 250 nsec of equilibration
dna_prot=subset(dna_prot,Time>max(dna_prot$Time)/4)


################
#Here we need to filter out HB and SB from contacts, and recalcualte WM interactions.
############
#data frame SQL-like magic
a=split(dna_prot_cryst,dna_prot_cryst$type)
#Important in raw df C-means all contacts
#in average we wnated C -to be contacts without HB and SB
#now we understand that we still need original C, but all downstream scripts already 
# assume C is without HB and SB
#we want to transition to new nomenclature
# VdW - contacts without HB and SB
# Z - all contacts
# to transition we retain all C type, but add Z and VdW
#get rid of duplicated in contacts via tricky merge and join
m=merge(a$SC,a$SB,by.x=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','z'))
a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(11:20)]
a$vdW$type<-factor('vdW')

#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

dna_prot_cryst=do.call(rbind,a)

# head(subset(dna_prot_cryst,type=='SB'))

a=split(dna_prot,dna_prot$type)
#get rid of duplicated in contacts via tricky merge and join
m=merge(a$SC,a$SB,by.x=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','z'))
a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(12:20)]
# a$VdW=a$C
a$vdW$type<-factor('vdW')

#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

dna_prot=do.call(rbind,a)

#-----------


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
phosphate<-c('O3\'','O5\'','O1P','O2P','P','OP1','OP2')
protbb<-c('C','O','N','CA','O','OT1','OT2','OXT')

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


if(df['DNA_resname']%in%c('ADE','DA','DA3','DA5')) {
	if(df['DNA_atom'] %in% c('N6','C6','C5','N7','C8')){ return('major')}
	else if(df['DNA_atom'] %in% c('C2','N3','C4')){ return('minor')}
	else return('NA')}
if(df['DNA_resname']%in%c('GUA','DG','DG3','DG5')) {
	if(df['DNA_atom'] %in% c('O6','C6','C5','N7','C8')){ return('major')}
	else if(df['DNA_atom'] %in% c('N2','C2','N3','C4')){ return('minor')}
	else return('NA')}
if(df['DNA_resname']%in%c('THY','DT','DT3','DT5')) {
	if(df['DNA_atom'] %in% c('C6','C5','C4','C5M','O4','C7')){ return('major')}
	else if(df['DNA_atom'] %in% c('N3','C2','O2')){ return('minor')}
	else return('NA')}

if(df['DNA_resname']%in%c('CYT','DC','DC3','DC5')) {
	if(df['DNA_atom'] %in% c('C6','C5','C4','N4')){ return('major')}
	else if(df['DNA_atom'] %in% c('N3','C2','O2')){ return('minor')}
	else return('NA')}

}


dna_prot_cryst$groove<-as.factor(apply(dna_prot_cryst[,c('DNA_atom','DNA_resname')],1,m_m))
dna_prot$groove<-as.factor(apply(dna_prot[,c('DNA_atom','DNA_resname')],1,m_m))



write.csv(dna_prot, file="../analysis_data/dna_prot_avr_df.csv")
write.csv(dna_prot_cryst, file="../analysis_data/dna_prot_avr_df_cryst.csv")


