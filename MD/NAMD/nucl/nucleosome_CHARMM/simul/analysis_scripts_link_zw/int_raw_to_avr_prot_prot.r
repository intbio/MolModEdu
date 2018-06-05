#R-script  conversts raw data frames to avarage data frames for protein-protein interactions
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


#########Starting with Protein-protein interactions#########
###############
##Loading data frames
prot_prot_cryst<-read.csv('../analysis_data/prot_prot_raw_df_cryst.csv')
prot_prot<-read.csv('../analysis_data/prot_prot_raw_df.csv')

################
#Here we need to filter out HB and IP from contacts, and recalcualte WM interactions.
############
#data frame SQL-like magic
a=split(prot_prot_cryst,prot_prot_cryst$type)

# a$Z=a$SC
# a$Z$type<-factor('Z')


#get rid of duplicated in contacts via tricky merge and join

m=merge(a$SC,a$SB,by.x=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','z'))
a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(11:20)]
a$vdW$type<-factor('vdW')

#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

prot_prot_cryst=do.call(rbind,a)

#MD data
# print('KUKU')

a=split(prot_prot,prot_prot$type)
# a$Z=a$SC
# a$Z$type<-factor('Z')

#get rid of duplicated in contacts via tricky merge and join
m=merge(a$SC,a$SB,by.x=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),by.y=c("Time","PROT1_atom","PROT1_chain","PROT1_resid","PROT2_atom","PROT2_chain","PROT2_resid"),all.x=TRUE,suffix=c('','z'))
a$vdW=m2[is.na(m2$typez) & is.na(m2$typey),-(12:20)]
a$vdW$type<-factor('vdW')

#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

prot_prot=do.call(rbind,a)

#-----------


#let's do averaging and combination
#########
nf=length(table(prot_prot$Time))
#nf=750

####
cn=colnames(prot_prot)
avl=cn[!(cn=='Time' | cn=='param1' | cn=='param2' | cn=='param3')]

prot_prot=ddply(prot_prot,avl,summarize,av_num=length(param1)/nf)


#############
#adding important name columns, resnames, etc.

rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


# sugar<-c('C1\'','C2\'','C3\'','C4\'','C5\'','O4\'')
# phosphate<-c('O3\'','O5\'','O1P','O2P','P')
protbb<-c('C','O','N','CA','OT1','OT2')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SB','HB','vdW','IM','WM','SC')

#We want to know for every interaction type get number of interactions between different types
# of amino acids and PROT1 groups
# dp <- function(name){

# if(name %in% sugar) { return('sugar') }
# else if(name %in% phosphate){return('phosphate')}
# else {return('base')}

# }

pp <- function(name){

if(name %in% protbb) { return('backbone') }
else {return('side chain')}

}

#Let's add PROT1 parts classification
prot_prot_cryst$PROT1_part<-as.factor(sapply(prot_prot_cryst$PROT1_atom,pp))
prot_prot$PROT1_part<-as.factor(sapply(prot_prot$PROT1_atom,pp))

#Let's add prtoein parts classification
prot_prot_cryst$PROT2_part<-as.factor(sapply(prot_prot_cryst$PROT2_atom,pp))
prot_prot$PROT2_part<-as.factor(sapply(prot_prot$PROT2_atom,pp))


#Let's add PROT1 bases names and protein residies names
colnames(rr)<-c('PROT1_chain','PROT1_resid','PROT1_resname')
prot_prot_cryst=merge(prot_prot_cryst,rr)
prot_prot=merge(prot_prot,rr)

colnames(rr)<-c('PROT2_chain','PROT2_resid','PROT2_resname')
prot_prot_cryst=merge(prot_prot_cryst,rr)
prot_prot=merge(prot_prot,rr)

#Let's assign level order
prot_prot_cryst$type<-factor(prot_prot_cryst$type,levels=type_lev)
prot_prot$type<-factor(prot_prot$type,levels=type_lev)

prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
prot_prot_cryst$PROT1_part<-factor(prot_prot_cryst$PROT1_part,levels=c('backbone','side chain'))


prot_prot_cryst$PROT2_resname<-factor(prot_prot_cryst$PROT2_resname,levels=prot_rn_lev)
prot_prot$PROT2_resname<-factor(prot_prot$PROT2_resname,levels=prot_rn_lev)

#Now let's add contact classification
#PP, PN, NP

pol_nonpol <- function(df){

if(!(df['type']%in%c('SC','vdW'))) { return('NA') }
else {
if(substring(df['PROT1_atom'],1,1) %in% c('O','P','N')){
if(substring(df['PROT2_atom'],1,1) %in% c('C','S')) {return('PN')}
else{return('PP')}
}
else{
if(substring(df['PROT2_atom'],1,1) %in% c('C','S')) {return('NN')}
else{return('PN')}

}
}
}

#Let's add pol-nonplo classification
prot_prot_cryst$C_type<-as.factor(apply(prot_prot_cryst[,c('PROT1_atom','PROT2_atom','type')],1,pol_nonpol))
prot_prot$C_type<-as.factor(apply(prot_prot[,c('PROT1_atom','PROT2_atom','type')],1,pol_nonpol))


write.csv(prot_prot, file="../analysis_data/prot_prot_avr_df.csv")
write.csv(prot_prot_cryst, file="../analysis_data/prot_prot_avr_df_cryst.csv")


