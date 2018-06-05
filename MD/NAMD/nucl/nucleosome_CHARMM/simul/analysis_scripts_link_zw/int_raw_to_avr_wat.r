#R-script  conversts raw data frames to avarage data frames
# water df is already averaged, so be add here some additional info

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


#########Starting with DNA-intein interactions#########
###############
##Loading data frames
#DNA-intein interactions
# wat_int<-read.csv('../analysis_data/wat_int_raw_df.csv')
wat_int_cryst<-read.csv('../analysis_data/wat_int_raw_df_cryst.csv')
wat_int<-read.csv('../analysis_data/wat_int_raw_df.csv')

################


#let's do and combination
#########

# ####
# cn=colnames(wat_int)
# avl=cn[!(cn=='Time' | cn=='param1' | cn=='param2')]

# wat_int=ddply(wat_int,avl,summarize,av_num=length(param1)/nf)


#############
#adding important name columns, resnames, etc.

rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)


sugar<-c('C1\'','C2\'','C3\'','C4\'','C5\'','O4\'')
phosphate<-c('O3\'','O5\'','O1P','O2P','P')
intbb<-c('C','O','N','CA')

#Let's define levels order for better graphics
int_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS','ADE','THY','GUA','CYT')

#Let's add DNA bases names and intein residies names
colnames(rr)<-c('INT_chain','INT_resid','INT_resname')
wat_int_cryst=merge(wat_int_cryst,rr)
wat_int=merge(wat_int,rr)

wat_int_cryst$INT_resname<-factor(wat_int_cryst$INT_resname,levels=int_rn_lev)
wat_int$INT_resname<-factor(wat_int$INT_resname,levels=int_rn_lev)


#We want to know for every interaction type get number of interactions between different types
# of amino acids and DNA groups
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
#Let's add  parts classification
wat_int_cryst$INT_part<-as.factor(mapply(dp,wat_int_cryst$INT_resname,wat_int_cryst$INT_atom))
wat_int$INT_part<-as.factor(mapply(dp,wat_int$INT_resname,wat_int$INT_atom))


wat_int$INT_part<-factor(wat_int$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))
wat_int_cryst$INT_part<-factor(wat_int_cryst$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))


write.csv(wat_int, file="../analysis_data/wat_int_avr_df.csv")
write.csv(wat_int_cryst, file="../analysis_data/wat_int_avr_df_cryst.csv")


