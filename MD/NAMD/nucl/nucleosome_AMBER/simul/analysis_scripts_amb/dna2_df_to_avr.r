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

dna<-read.csv('../analysis_data/dna2_param_df_md.csv',na.strings=c("NA",'---'))
nf=length(table(dna$Time))

#get pairing
dna_pair=ddply(dna,'BPnum',summarize,Pairing_av=mean(Pairing))

#Important - we should remove values were pairing is 0, and adjacent points
dna=subset(dna,Pairing==1&myshift(Pairing==1))
#Magic averaging
dna_aver=ddply(dna,'BPnum',function(df) apply(df[,lapply(df,class)%in%c('numeric')],2,mean))
dna_var=ddply(dna,'BPnum',function(df) apply(df[,lapply(df,class)%in%c('numeric')],2,sd))
dna_err=ddply(dna,'BPnum',function(df) apply(df[,lapply(df,class)%in%c('numeric')],2,errest,nf))
dna1=merge(dna_aver,dna_var,by=c('BPnum'),suffixes = c("","_sd"))
dna2=merge(dna1,dna_err,by=c('BPnum'),suffixes = c("_av","_er"))
dna_all=merge(dna2,dna_pair,by=c('BPnum'),suffixes = c("",""))
# dna_all=dna2
#add naming
rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)
#let's get the sequence
#get sequence for CHAIN I
chI_seq=sapply(rr[rr$chain=='CHI','resname'],substring,1,1)
chJ_seq=rev(sapply(rr[rr$chain=='CHJ','resname'],substring,1,1))

dna_all$SEQ_1=chI_seq
dna_all$SEQ_2=chJ_seq

write.csv(dna_all, file="../analysis_data/dna2_param_df_md_avr.csv")

quit()
