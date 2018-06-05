#R-script to analyze DNA-param data
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
##############
#Some useful functions from cookbook-r.com
#############
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#############
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
#############

#Let's start the code,
#it is sequence and nucleosome dependent, so look through the code before applying it!
#############
#set bigger font
theme_set(theme_grey(base_size=40))


#####
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
#seqnumbering starts from 0 here
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
#Get crystal and simulation data
dfcryst<-read.csv("../analysis_data/dna_param_cryst.csv")
#dfcryst$X<-as.factor(dfcryst$X)
df<-read.csv("../analysis_data/dna_param_big_df.csv")
#df$X<-as.factor(df$X)
#Let's add a sequence column


#Now let's analyze where BP become unpaired - plot the df$BP - column
theme_set(theme_grey())
#png(file="../analysis_data/dna_pairing.png",width=1000,height=800)


####Now we'll build all the profiles with SD
pl=list()
theme_set(theme_grey(base_size=40))

for (i in c("Shear","Stagger","Stretch","Buckle","Opening","Prop.Tw","Rise","Shift","Slide","Roll","Tilt","Twist")) {


if(i %in% c("Buckle","Opening","Prop.Tw","Stagger","Stretch","Shear")) {
shift=73.0
} else {
shift=73.5
}
print(shift)
data<-summarySE(df,measurevar=i, groupvars=c("X"))
data<-data[4:143,]
data$X=data$X-shift
dfcryst_t<-dfcryst
dfcryst_t$X=dfcryst_t$X-shift
seqdf_t<-seqdf
seqdf_t$X<-seqdf_t$X-73.0
seqdf_t$Y<-min(dfcryst[[i]],data[[i]]-data$sd,na.rm=TRUE)

cor=cor(x=data[,i],y=dfcryst[4:143,i], method = c("pearson"))
str2<-paste("Correlation coef. =",round(cor,digits=2), collapse = '') 
print(cor)

pl[[i]]<-ggplot(data=data,aes_string(x='X',y=i))+
#geom_point(size=3)+
geom_line(size=3)+
#geom_errorbar(ymin=data[[i]]-data$sd,ymax=data[[i]]+data$sd, color="blue",size=3)+
#geom_errorbar(ymin=data[[i]]-data$se,ymax=data[[i]]+data$se, color="black",size=3)+
xlab("Base pair / base pair step")+geom_line(data=dfcryst_t,aes_string(x='X',y=i),color="red",size=3)+
ggtitle(paste(i,"profile"))+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+#+xlim(-70,70)
geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)+
annotate("text", x=Inf, y=Inf, label = str2,hjust=1,vjust=1,size=30)

}
png(file="../analysis_data/dna_param_profiles_aver.png",width=8000,height=4000)
multiplot(plotlist=pl,cols=2)
dev.off()


#plot.new()
# z <- locator(1)


