#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)

df<-read.csv("../analysis_data/dna_rot_df_all.csv",skip=0,header=TRUE,check.name=FALSE)

df_cryst=subset(df,Time==0)[,c('Basepair','Angle')]
df_avr=subset(df,Time==1)[,c('Basepair','Angle')]

# df=subset(df,Time>(max(Time)/4))
# nf=length(table(df$Time))
# print(nf)
# df_avr=ddply(df,c('Basepair'),summarize,Angle=mean(Angle))

# head(df_avr,n=100)




write.csv(df_avr, file="../analysis_data/dna_rot_df_avr.csv")
write.csv(df_cryst, file="../analysis_data/dna_rot_df_cryst.csv")

