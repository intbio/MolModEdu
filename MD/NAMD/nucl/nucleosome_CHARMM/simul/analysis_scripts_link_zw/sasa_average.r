#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)

df_ref<-read.csv("../analysis_data/sasa_ref.csv",skip=0,header=TRUE,check.name=FALSE)
df_ref$Resname=revalue(df_ref$Resname, c("HIS"="HSE"))
#dfcryst$X<-as.factor(dfcryst$X)
# print(df_ref)
df<-read.csv("../analysis_data/sasa_nucl.csv",skip=0,header=TRUE,check.name=FALSE)
df=subset(df,Time>(max(Time)/4))
nf=length(table(df$Time))
print(nf)
df_avr=ddply(df,c('Chain','Resid','Resname'),summarize,sasa_avr=mean(SASA))
df_avr=merge(df_avr,df_ref)
df_avr$sasa_fract=df_avr$sasa_avr/df_avr$SASA_ref

head(df_avr,n=100)

write.csv(df_avr, file="../analysis_data/sasa_nucl_avr.csv")

