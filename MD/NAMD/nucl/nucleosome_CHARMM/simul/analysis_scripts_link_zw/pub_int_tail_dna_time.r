#Qunatifies the evolution of interactions of tails with DNA
# % of amino acid residues in contact with DNA

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



# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')

# dna_prot=subset(dna_prot,type %in%c('SC','WM'))

# #tails
# # PROT=atomsel("noh and ((segname CHA CHE and resid < 44 ) or (segname CHB CHF and resid < 24) or (segname CHC CHG and (resid <16 or resid > 118)) or (segname CHD CHH and resid < 30))")
# #Let's get only tails
# #New def
# #set core_ca "((segname CHA CHE and (resid 37 to 135)) or (segname CHB CHF and (resid 16 to 102))  or (segname CHC CHG and (resid 12 to 118)) or (segname CHD CHH and (resid 21 to 122))) and name CA"

# # dna_prot=subset(dna_prot, ((PROT_chain %in% c('CHA','CHE')) & (PROT_resid < 37)) | ((PROT_chain %in% c('CHB','CHF')) & (PROT_resid < 16)) | ((PROT_chain %in% c('CHC','CHG')) & ((PROT_resid < 12 | PROT_resid > 118)) ) | ((PROT_chain %in% c('CHD','CHH')) & (PROT_resid < 21))) 
# dna_prot=subset(dna_prot, ((PROT_chain %in% c('CHA','CHE')) & (PROT_resid < 44)) | ((PROT_chain %in% c('CHB','CHF')) & (PROT_resid < 24)) | ((PROT_chain %in% c('CHC','CHG')) & ((PROT_resid < 16 | PROT_resid > 118)) ) | ((PROT_chain %in% c('CHD','CHH')) & (PROT_resid < 30))) 

# #now let's reduce to individual residues
# tail_len=data.frame(PROT_chain=c('CHA','CHE','CHB','CHF','CHC','CHG','CHD','CHH'),resnum=c(43,43,23,23,25,25,29,29))
# tailh_len=data.frame(Histones=c('H3','H4','H2A','H2B'),resnum=c(43,23,25,29))

# # tail_len=data.frame(PROT_chain=c('CHA','CHE','CHB','CHF','CHC','CHG','CHD','CHH'),resnum=c(36,36,15,15,21,21,20,20))
# # tailh_len=data.frame(Histones=c('H3','H4','H2A','H2B'),resnum=c(36,15,21,20))


# dna_prot=subset(dna_prot,Time%in%seq(0,1000,10))

# cn=colnames(dna_prot)
# avl=c('Time','PROT_chain','PROT_resid')

# dna_prot=ddply(dna_prot,avl,summarize,num=length(param1))
# avl=c('Time','PROT_chain')
# dna_prot=ddply(dna_prot,avl,summarize,num=length(PROT_resid))

# ch_types=data.frame(PROT_chain=c('CHA','CHE','CHB','CHF','CHC','CHG','CHD','CHH'),Histones=c('H3','H3','H4','H4','H2A','H2A','H2B','H2B'))
# dna_prot=merge(dna_prot,ch_types)
# dna_prot=ddply(dna_prot,c('Time','Histones'),summarize,num=sum(num))
# dna_tail_int=merge(dna_prot,tailh_len)

# write.csv(dna_tail_int, file="../analysis_data/dna_tail_int.csv")

# theme_set(theme_bw()+theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank()))
theme_set(theme_bw())


dna_prot<-read.csv('../analysis_data/dna_prot_raw_df_cryst.csv')

dna_prot=subset(dna_prot,type %in%c('SC','WM'))

#tails
# PROT=atomsel("noh and ((segname CHA CHE and resid < 44 ) or (segname CHB CHF and resid < 24) or (segname CHC CHG and (resid <16 or resid > 118)) or (segname CHD CHH and resid < 30))")
#Let's get only tails
#New def
#set core_ca "((segname CHA CHE and (resid 37 to 135)) or (segname CHB CHF and (resid 16 to 102))  or (segname CHC CHG and (resid 12 to 118)) or (segname CHD CHH and (resid 21 to 122))) and name CA"

# dna_prot=subset(dna_prot, ((PROT_chain %in% c('CHA','CHE')) & (PROT_resid < 37)) | ((PROT_chain %in% c('CHB','CHF')) & (PROT_resid < 16)) | ((PROT_chain %in% c('CHC','CHG')) & ((PROT_resid < 12 | PROT_resid > 118)) ) | ((PROT_chain %in% c('CHD','CHH')) & (PROT_resid < 21))) 
dna_prot=subset(dna_prot, ((PROT_chain %in% c('CHA','CHE')) & (PROT_resid < 44)) | ((PROT_chain %in% c('CHB','CHF')) & (PROT_resid < 24)) | ((PROT_chain %in% c('CHC','CHG')) & ((PROT_resid < 16 | PROT_resid > 118)) ) | ((PROT_chain %in% c('CHD','CHH')) & (PROT_resid < 30))) 

# #now let's reduce to individual residues
tail_len=data.frame(PROT_chain=c('CHA','CHE','CHB','CHF','CHC','CHG','CHD','CHH'),resnum=c(43,43,23,23,25,25,29,29))
tailh_len=data.frame(Histones=c('H3','H4','H2A','H2B'),resnum=c(43,23,25,29))

# tail_len=data.frame(PROT_chain=c('CHA','CHE','CHB','CHF','CHC','CHG','CHD','CHH'),resnum=c(36,36,15,15,21,21,20,20))
# tailh_len=data.frame(Histones=c('H3','H4','H2A','H2B'),resnum=c(36,15,21,20))


# dna_prot=subset(dna_prot,Time%in%seq(0,1000,10))

cn=colnames(dna_prot)
avl=c('PROT_chain','PROT_resid')

dna_prot=ddply(dna_prot,avl,summarize,num=length(param1))
avl=c('PROT_chain')
dna_prot=ddply(dna_prot,avl,summarize,num=length(PROT_resid))

ch_types=data.frame(PROT_chain=c('CHA','CHE','CHB','CHF','CHC','CHG','CHD','CHH'),Histones=c('H3','H3','H4','H4','H2A','H2A','H2B','H2B'))
dna_prot=merge(dna_prot,ch_types)
dna_prot=ddply(dna_prot,c('Histones'),summarize,num=sum(num))
dna_tail_int_cryst=merge(dna_prot,tailh_len)
dna_tail_int_cryst$Time=0
write.csv(dna_tail_int_cryst, file="../analysis_data/dna_tail_int_cryst.csv")



dna_tail_int<-read.csv('../analysis_data/dna_tail_int.csv')
dna_tail_int_cryst<-read.csv('../analysis_data/dna_tail_int_cryst.csv')

dna_tail_int=rbind(subset(dna_tail_int,Time>0),dna_tail_int_cryst)

dna_tail_int=arrange(dna_tail_int,Histones)

q<-ggplot(data=dna_tail_int,aes(x=Time,y=num/resnum/2*100,color=Histones)) + geom_line(size=0.75)+
xlab('Time, ns')+ylab('Interacting residues, %')+xlim(0,1000)+scale_y_continuous(limits=c(10,90),breaks = seq(0,100, by = 10))+
scale_color_manual(values=c('H3'="blue",'H4'='green','H2A'='gold','H2B'='red'),breaks=c('H3','H4','H2A','H2B'),name='Histone tails')+
geom_point(data=dna_tail_int_cryst)

ggsave(filename="../analysis_data/pub_dna_tails_time.png",plot=q,width=7,height=3)
q()

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
#-----------


#let's do averaging and combination
#########
nf=length(table(dna_prot$Time))
#nf=750

####
# cn=colnames(dna_prot)
# avl=cn[!(cn=='Time' | cn=='param1' | cn=='param2' | cn=='param3')]

# dna_prot=ddply(dna_prot,avl,summarize,av_num=length(param1)/nf)


#############
#adding important name columns, resnames, etc.

# rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


#Now let's add contact classification
#PP, PN, NP


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


