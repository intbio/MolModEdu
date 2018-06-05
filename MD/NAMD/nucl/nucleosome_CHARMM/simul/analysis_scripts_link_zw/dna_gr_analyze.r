#R-script to analyze DNA-param data of grooves from Curves+
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
##############
###############

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm),
          min   = min     (xx[[col]], na.rm=na.rm),
          max   = max     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

#############


###############
#####
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
#seqnumbering starts from 0 here
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
#Get crystal and simulation data
dfcryst<-read.csv("../analysis_data/dna_gr_param_cryst.csv")
#dfcryst$X<-as.factor(dfcryst$X)
df<-read.csv("../analysis_data/dna_gr_param_big_df.csv")
#df$X<-as.factor(df$X)

dfc<-melt(dfcryst,id.vars="Level",measure.vars=c('W12','W21'))

# ggplot(data=dfc,aes(x=Level,y=value,color=variable))+geom_line()+ggtitle("Width of DNA grooves")


data<-summarySE(df,measurevar="W12", groupvars=c("Level"),na.rm=TRUE)
data2<-summarySE(df,measurevar="W21", groupvars=c("Level"),na.rm=TRUE)

ggplot(data=data,aes(x=Level,y=W12))+geom_line(color="blue",size=2)+
geom_line(data=data2,aes(x=Level,y=W21),color="green",size=2)+
#geom_line(data=data2,aes(x=Level,y=min),color="green",size=1)+
#geom_line(data=data2,aes(x=Level,y=max),color="green",size=1)+
#geom_line(data=data,aes(x=Level,y=max),color="blue",size=1)+
#geom_line(data=data,aes(x=Level,y=min),color="blue",size=1)+
geom_line(data=dfcryst,aes(x=Level,y=W12),color="red")+
geom_line(data=dfcryst,aes(x=Level,y=W21),color="red")+ggtitle("Width of DNA grooves")
# geom_line(data=dfcryst[dfcryst$Strand==1,],aes(x=Resid,y=B),color="red",size=2)+
# ggtitle("Epsilon-Zeta distribution in strand I")
#geom_errorbar(ymin=data$min,ymax=data$max, color="blue",size=1)
ggsave("../analysis_data/dna_grooves.png",height=6,width=14)

# plot.new()
# z <- locator(1)


