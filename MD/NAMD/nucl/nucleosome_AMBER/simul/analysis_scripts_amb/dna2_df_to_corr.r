#This script will produce average 
#correlation analysis for parameter data

library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(foreach)
library(doParallel)
# library(gridExtra)
##############
FREQ=10 #number of data points in 1 ns
errest<-function(x,numf){

return(sd(x)/sqrt(numf))
}

myshift<-function(x){ #fits the vector to the left by 1, add 1 on end
return(c(x[-1],1))

}

dna<-read.csv('../analysis_data/dna2_param_df_md.csv',na.strings=c("NA",'---'))
nf=length(table(dna$Time))

#Important - we should remove values were pairing is 0, and adjacent points
#We can not do this here, since ACF won't work!!!
# dna=subset(dna,Pairing==1&myshift(Pairing==1))


#Now we are ready to calculate correlations	
# <dTW(i), dTW(j)> 
# <dRL(i), dRL(j)> 
# <dSL(i), dSL(j)> 
# For all pairs (i,j), i.e., matrices 146 x 146. 
# [here d is for delta; I hope we understand each other]
# If possible, I would like to have access to numeric data.
# I expect that TW(i) and TW(i+1) should be anti-correlated, but I don’t know what to expect for TW(i) and TW(i+5)… the same for ROLL & SLIDE.

# Also, please calculate
# <dTW(i), dRL(i)> 
# <dTW(i), dSL(i)> 
# <dSL(i), dRL(i)> 

###############################################
#But we need to start with auto correlations

#We will have a data frame with 147*N parameters and their autocor times.


get_act<-function(x){ #fits the vector to the left by 1, add 1 on end
if(NA%in%x){return(NA)}
else{
xt=ts(x,frequency=FREQ)
a=acf(xt,plot=FALSE)$acf
t=seq(length(a))/FREQ
res=nls(a~exp(-t/tau),start=list(tau=1),control=list(minFactor=1e-12,warnOnly = TRUE))
# print(coef(res)[[1]])
return(coef(res)[[1]])
}
}


#Magic calculation
dna_act=ddply(dna,'BPnum',function(df) apply(df[,lapply(df,class)%in%c('numeric')&!(names(df)%in%c('Pairing','Level','Time'))],2,failwith(NA,get_act)))

dna_all=dna_act



#add naming
rr<-read.csv('../analysis_data/resid_resname_df.csv',stringsAsFactors=FALSE)
#let's get the sequence
#get sequence for CHAIN I
chI_seq=sapply(rr[rr$chain=='CHI','resname'],substring,1,1)
chJ_seq=rev(sapply(rr[rr$chain=='CHJ','resname'],substring,1,1))

dna_all$SEQ_1=chI_seq
dna_all$SEQ_2=chJ_seq

write.csv(dna_all, file="../analysis_data/dna2_param_df_md_act.csv")

################################
#Now cross-correlations
#This should be a very big matrix
#We want to have cross-correlation coefficients and probably some parameters of cross-correlation functions.
#Columns BPnum1, BPnum2, Param1, Param2, cross-corr


get_cor<-function(df1,df2){ #fits the vector to the left by 1, add 1 on end
	# print(names(df))
	# print(nrow(df1))
	# return(ddply(df2,c('BPnum2','variable2'),summarize,corr=cor(df1$value,value)))
	return(ddply(df2,c('BPnum2','variable2'),function(df) data.frame(corr=cor(df1$value,df$value,use='na.or.complete'))))
# return(data.frame(corr=c(mean(df1[,'value']),10)))
}
# function(df) apply(df[,lapply(df,class)%in%c('numeric')],2,mean)
#get rid of 
#where Pairing in not 1 we need to put NA
# dna=subset(dna,Pairing==1&myshift(Pairing==1))
dna[dna$Pairing==0|myshift(dna$Pairing==0),!names(dna)%in%c('BPnum','Time','Pairing')]=NA


d2=dna[,(!(names(dna)%in%c('Pairing','Level')))]
# d2=subset(d2,Time<100&BPnum<40)
dna1=melt(d2,id.vars=c('BPnum','Time'),measure.vars=c('x','y','z','Roll','Twist','Slide','chi_1','chi_2','P_1','P_2'))

# dna1=melt(d2,id.vars=c('BPnum','Time'),measure.vars=c('Roll','Twist','Slide'))
dna2=rename(dna1, c('BPnum'='BPnum2','variable'='variable2'))
dna1=rename(dna1, c('BPnum'='BPnum1','variable'='variable1'))

workers <- makeCluster(4) # My computer has 2 cores
registerDoParallel(workers)

dna_cor=ddply(dna1,c('BPnum1','variable1'),get_cor,dna2,.progress='text',.parallel=FALSE)
write.csv(dna_cor, file="../analysis_data/dna2_param_df_md_cor.csv")


quit()
