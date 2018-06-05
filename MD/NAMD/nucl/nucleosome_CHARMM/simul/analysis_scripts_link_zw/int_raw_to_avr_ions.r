#R-script  conversts raw data frames to avarage data frames
# makes important filtering - eg. deletes duplicates of HB and IP from contacts
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


#########Starting with ION-intein interactions#########
###############
##Loading data frames
#ION-intein interactions
# ion_int<-read.csv('../analysis_data/ion_int_raw_df.csv')
ion_int_cryst<-read.csv('../analysis_data/ions_int_raw_df_cryst.csv')
ion_int<-read.csv('../analysis_data/ions_int_raw_df.csv')

################
#Here we need to filter out HB and IP from contacts, and recalcualte WM interactions.
############
#data frame SQL-like magic
# a=split(ion_int_cryst,ion_int_cryst$type)

# a$Z=a$C
# a$Z$type<-factor('Z')
#get rid of duplicated in contacts via tricky merge and join
# m=merge(a$C,a$IP,by.x=c('ION_chain','ION_resid','ION_atom',"INT_chain","INT_resid"),by.y=c('ION_chain','ION_resid','ION_atom',"INT_chain","INT_resid"),all.x=TRUE,suffix=c('','y'))
# m2=merge(m,a$HB,by.x=c('ION_chain','ION_resid','ION_atom',"INT_chain","INT_resid"),by.y=c('ION_chain','ION_resid','ION_atom',"INT_chain","INT_resid"),all.x=TRUE,suffix=c('','z'))
# a$C=m2[is.na(m2$typez) & is.na(m2$typey),-(10:20)]
# a$C$type<-factor('VdW')

#Let's simlify now the water mediated interactions
# a$WM=a$WM[!duplicated(a$WM),]

# ion_int_cryst=do.call(rbind,a)

# a=split(ion_int,ion_int$type)
# a$Z=a$C
# a$Z$type<-factor('Z')
#get rid of duplicated in contacts via tricky merge and join
# m=merge(a$C,a$IP,by.x=c("Time","ION_atom","ION_ind","INT_atom","INT_chain","INT_resid"),by.y=c("Time","ION_atom","ION_ind","INT_atom","INT_chain","INT_resid"),all.x=TRUE,suffix=c('','y'))
# m2=merge(m,a$HB,by.x=c("Time","ION_atom","ION_ind","INT_atom","INT_chain","INT_resid"),by.y=c("Time","ION_atom","ION_ind","INT_atom","INT_chain","INT_resid"),all.x=TRUE,suffix=c('','z'))
# a$C=m2[is.na(m2$typez) & is.na(m2$typey),-(11:20)]
# a$C$type<-factor('VdW')

#Let's simlify now the water mediated interactions
# a$WM=a$WM[!duplicated(a$WM),]

# ion_int=do.call(rbind,a)

#-----------


#let's do averaging and combination
#########
nf=length(table(ion_int$Time))
#nf=750

####
cn=colnames(ion_int)
avl=cn[!(cn=='Time' | cn=='param1' | cn=='param2' | cn=='param3')]

ion_int=ddply(ion_int,avl,summarize,av_num=length(param1)/nf)


#############
#adding important name columns, resnames, etc.

rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


sugar<-c('C1\'','C2\'','C3\'','C4\'','C5\'','O4\'')
phosphate<-c('O3\'','O5\'','O1P','O2P','P')
intbb<-c('C','O','N','CA')

#Let's define levels order for better graphics
int_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS','ADE','THY','GUA','CYT')
type_lev=c('IP','HB','VdW','IM','WM','Z')



colnames(rr)<-c('INT_chain','INT_resid','INT_resname')
ion_int_cryst=merge(ion_int_cryst,rr)
ion_int=merge(ion_int,rr)


#We want to know for every interaction type get number of interactions between different types
# of amino acids and ION groups
dp <- function(resname,atname){
if(resname %in% c('ADE','THY','GUA','CYT'))
{
if(atname %in% sugar) { return('sugar') }
else if(atname %in% phosphate){return('phosphate')}
else {return('base')}

}
else{
if(atname %in% intbb) { return('backbone') }
else {return('side chain')}
}

}

#Let's add prtoein parts classification
ion_int_cryst$INT_part<-as.factor(mapply(dp,ion_int_cryst$INT_resname,ion_int_cryst$INT_atom))
ion_int$INT_part<-as.factor(mapply(dp,ion_int$INT_resname,ion_int$INT_atom))


#Let's assign level order
ion_int_cryst$type<-factor(ion_int_cryst$type,levels=type_lev)
ion_int$type<-factor(ion_int$type,levels=type_lev)

ion_int$INT_part<-factor(ion_int$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))
ion_int_cryst$INT_part<-factor(ion_int_cryst$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))


ion_int_cryst$INT_resname<-factor(ion_int_cryst$INT_resname,levels=int_rn_lev)
ion_int$INT_resname<-factor(ion_int$INT_resname,levels=int_rn_lev)

#Now let's add contact classification
#PP, PN, NP

pol_nonpol <- function(df){

if(!(df['type']%in%c('C','VdW','Z'))) { return('NA') }
else {
if(substring(df['INT_atom'],1,1) %in% c('C','S')) {return('NONPOL')}
else{return('POL')}
}

}


#Let's add pol-nonplo classification
ion_int_cryst$C_type<-as.factor(apply(ion_int_cryst[,c('INT_atom','type')],1,pol_nonpol))
ion_int$C_type<-as.factor(apply(ion_int[,c('INT_atom','type')],1,pol_nonpol))


#Let's add groove classification

# A: major: N6, C6, C5, N7, C8, minor:  C2,N3,C4
# G: major: O6, C6, C5, N7, C8, minor: N2, C2,N3,C4
# T: major: C6, C5, C4, C5M, O4 minor: N3, C2,O2
# C: major: C6, C5, C4, N4 minor: N3, C2,O2

m_m <- function(df){

if(!(df['INT_resname']%in%c('ADE','GUA','CYT','THY'))) {return('NA')}
if(df['INT_resname']=='ADE') {
	if(df['INT_atom'] %in% c('N6','C6','C5','N7','C8')){ return('major')}
	else if(df['INT_atom'] %in% c('C2','N3','C4')){ return('minor')}
	else return('NA')}
if(df['INT_resname']=='GUA') {
	if(df['INT_atom'] %in% c('O6','C6','C5','N7','C8')){ return('major')}
	else if(df['INT_atom'] %in% c('N2','C2','N3','C4')){ return('minor')}
	else return('NA')}
if(df['INT_resname']=='THY') {
	if(df['INT_atom'] %in% c('C6','C5','C4','C5M','O4')){ return('major')}
	else if(df['INT_atom'] %in% c('N3','C2','O2')){ return('minor')}
	else return('NA')}

if(df['INT_resname']=='CYT') {
	if(df['INT_atom'] %in% c('C6','C5','C4','N4')){ return('major')}
	else if(df['INT_atom'] %in% c('N3','C2','O2')){ return('minor')}
	else return('NA')}

}


ion_int_cryst$groove<-as.factor(apply(ion_int_cryst[,c('INT_atom','INT_resname')],1,m_m))
ion_int$groove<-as.factor(apply(ion_int[,c('INT_atom','INT_resname')],1,m_m))





write.csv(ion_int, file="../analysis_data/ion_int_avr_df.csv")
write.csv(ion_int_cryst, file="../analysis_data/ion_int_avr_df_cryst.csv")


