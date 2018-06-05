#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)


#dfcryst$X<-as.factor(dfcryst$X)
df<-read.csv("../analysis_data/ss_nucl.csv",skip=0,header=TRUE,check.name=FALSE)
df=subset(df,Time>250)
nf=length(table(df$Time))
df_avr=ddply(df,c('Chain','resid'),summarize,ss_prop=table(SS),ss_type=names(table(SS)))
colnames(df_avr)<-c('Chain','Resid','ss_prop','ss_type')
df_avr$ss_prop=df_avr$ss_prop/nf

df_avr[df_avr$ss_type=='C','ss_prop']=df_avr[df_avr$ss_type=='C','ss_prop']+df_avr[df_avr$ss_type=='T','ss_prop']
df_avr[df_avr$ss_type=='B','ss_prop']=df_avr[df_avr$ss_type=='B','ss_prop']+df_avr[df_avr$ss_type=='E','ss_prop']
df_avr=subset(df_avr,ss_type%in%c('H','G','B','C'))

write.csv(df_avr, file="../analysis_data/ss_nucl_avr.csv")

