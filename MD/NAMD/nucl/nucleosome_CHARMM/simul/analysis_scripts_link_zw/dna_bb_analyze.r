#R-script to analyze DNA-param data of backbone from Curves+
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
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
dfcryst<-read.csv("../analysis_data/dna_bb_param_cryst.csv")
#dfcryst$X<-as.factor(dfcryst$X)
df<-read.csv("../analysis_data/dna_bb_param_big_df.csv")
#df$X<-as.factor(df$X)

#Let's determine the BI-BII conformation
#see in manual of Canal
#BI is epsilon-Zeta<0
#BII is epsilon-Zeta>0
dfcryst$B<-dfcryst$Epsil-dfcryst$Zeta
dfcryst$B[dfcryst$B>180 & !is.na(dfcryst$B)]=dfcryst$B[dfcryst$B>180 & !is.na(dfcryst$B)]-180
dfcryst$B[dfcryst$B<(-180) & !is.na(dfcryst$B)]=dfcryst$B[dfcryst$B<(-180) & !is.na(dfcryst$B)]+180

#ggplot(data=dfcryst[dfcryst$Strand==1,],aes(x=Resid,y=B))+geom_point()

df$B<-df$Epsil-df$Zeta
df$B[df$B>180 & !is.na(df$B)]=df$B[df$B>180 & !is.na(df$B)]-180
df$B[df$B<(-180) & !is.na(df$B)]=df$B[df$B<(-180) & !is.na(df$B)]+180

data<-summarySE(df,measurevar="B", groupvars=c("Resid","Strand"))

ggplot(data=df[data$Strand==1,],aes(x=Resid,y=B))+geom_point(color="blue")+
geom_line(data=data[data$Strand==1,],aes(x=Resid,y=B),color="green",size=2)+
geom_line(data=dfcryst[dfcryst$Strand==1,],aes(x=Resid,y=B),color="red",size=2)+
ggtitle("Epsilon-Zeta distribution in strand I")
#geom_errorbar(ymin=data[data$Strand==1,'min'],ymax=data[data$Strand==1,'max'], color="blue",size=1)
ggsave("../analysis_data/dna_b1_b2.png",height=6,width=14)

# plot.new()
# z <- locator(1)


