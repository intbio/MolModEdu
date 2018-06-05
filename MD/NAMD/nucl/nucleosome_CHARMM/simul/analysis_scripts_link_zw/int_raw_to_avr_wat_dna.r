#R-script  conversts raw data frames to avarage data frames
# in this case we will not do averaging but will assign names

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
# wat_dna<-read.csv('../analysis_data/wat_dna_raw_df.csv')
wat_dna_cryst<-read.csv('../analysis_data/wat_dna_int_raw_df_cryst.csv')
wat_dna<-read.csv('../analysis_data/wat_dna_int_raw_df.csv')

################
#Here we need to filter out HB and IP from contacts, and recalcualte WM interactions.
############

#-----------


#let's do averaging and combination
#########
# nf=length(table(wat_dna$Time))
#nf=750

####
# cn=colnames(wat_dna)
# avl=cn[!(cn=='Time' | cn=='param1' | cn=='param2' | cn=='param3')]

# wat_dna=ddply(wat_dna,avl,summarize,av_num=length(param1)/nf)


#############
#adding important name columns, resnames, etc.

rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


sugar<-c('C1\'','C2\'','C3\'','C4\'','C5\'','O4\'')
phosphate<-c('O3\'','O5\'','O1P','O2P','P')
protbb<-c('C','O','N','CA','OT1','OT2')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('HB','SC')

#We want to know for every interaction type get number of interactions between different types
# of amino acids and DNA groups
dp <- function(name){

if(name %in% sugar) { return('sugar') }
else if(name %in% phosphate){return('phosphate')}
else {return('base')}

}


#Let's add DNA parts classification
wat_dna_cryst$DNA_part<-as.factor(sapply(wat_dna_cryst$DNA_atom,dp))
wat_dna$DNA_part<-as.factor(sapply(wat_dna$DNA_atom,dp))


#Let's add DNA bases names and protein residies names
colnames(rr)<-c('DNA_chain','DNA_resid','DNA_resname')
wat_dna_cryst=merge(wat_dna_cryst,rr)
wat_dna=merge(wat_dna,rr)


#Major-minor groove

#Let's try to discriminate between major and minor groove

# A: major: N6, C6, C5, N7, C8, minor:  C2,N3,C4
# G: major: O6, C6, C5, N7, C8, minor: N2, C2,N3,C4
# T: major: C6, C5, C4, C5M, O4 minor: N3, C2,O2
# C: major: C6, C5, C4, N4 minor: N3, C2,O2

m_m <- function(df){

if(df['type']!='SC') { return('NA') }
else {
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
}

#Let's add pol-nonplo classification
wat_dna_cryst$groove<-as.factor(apply(wat_dna_cryst[,c('DNA_atom','DNA_resname','type')],1,m_m))
wat_dna$groove<-as.factor(apply(wat_dna[,c('DNA_atom','DNA_resname','type')],1,m_m))




#Let's assign level order
wat_dna_cryst$type<-factor(wat_dna_cryst$type,levels=type_lev)
wat_dna$type<-factor(wat_dna$type,levels=type_lev)

wat_dna$DNA_part<-factor(wat_dna$DNA_part,levels=c('phosphate','sugar','base'))
wat_dna_cryst$DNA_part<-factor(wat_dna_cryst$DNA_part,levels=c('phosphate','sugar','base'))



#Now let's add contact classification
#PP, PN, NP


write.csv(wat_dna, file="../analysis_data/wat_dna_avr_df.csv")
write.csv(wat_dna_cryst, file="../analysis_data/wat_dna_avr_df_cryst.csv")


