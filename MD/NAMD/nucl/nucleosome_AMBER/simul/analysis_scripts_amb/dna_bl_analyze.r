#R-script to analyze DNA-param data in a broken line approximation
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(plyr)
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


#Function to compute angles

ang_comp <- function(DNA_bp_centers,step) {

      # vectors=DNA_bp_centers[2:nrow(DNA_bp_centers),c('x','y','z')]-DNA_bp_centers[1:nrow(DNA_bp_centers)-1,c('x','y','z')]
      vectors=DNA_bp_centers[seq(step+1,nrow(DNA_bp_centers),1),c('x','y','z')]-DNA_bp_centers[seq(1,nrow(DNA_bp_centers)-step,1),c('x','y','z')]

      v1=vectors[seq(1,nrow(vectors)-step,1),]
      v2=vectors[seq(step+1,nrow(vectors),1),]
      ip=apply(v1*v2,1,sum)
      n1=apply(v1*v1,1,sum)
      n2=apply(v2*v2,1,sum)
      norm=sqrt(n1*n2)
      angles <- acos( ip / norm )*180/3.14159
      anglesdf <- data.frame(angles)
      #SHL should be the center of the points connected by vector
      anglesdf['SHL']=((seq(1,nrow(DNA_bp_centers)-step-step,1)+step)-74)/10
    return(anglesdf)
}



#############
#############

#Let's start the code,
#it is sequence and nucleosome dependent, so look through the code before applying it!
#############
#set bigger font
# theme_set(theme_grey(base_size=40))
theme_set(theme_grey())


#####
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
#seqnumbering starts from 0 here
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
#Get crystal and simulation data
dfcryst<-read.csv("../analysis_data/dna_bl_cryst.csv")

df<-read.csv("../analysis_data/dna_bl_big_df.csv")
data<-df[df['Time']>250,] 
#dfcryst$X<-as.factor(dfcryst$X)
# df<-read.csv("../analysis_data/dna_bl_big_df.csv")
#df$X<-as.factor(df$X)
#Let's add a sequence column
####################
#Now we will make plot at overall parameter distributions



ang_cryst<-ang_comp(dfcryst,10)

# ggplot(data=ang_cryst,aes(x=SHL,y=angles))+
# xlab("SHL of angle vertex")+ylab("Angle, deg")+
# geom_line() + 
# scale_x_continuous(breaks = round(seq(-7,7, by = 0.5),0.1)) + ggtitle("Angle between segments, 10 bp step")
# ggsave("../analysis_data/dna_bl_cryst_10.png",height=6,width=14)    
# dev.off()

#Now let's deal with simulation data

ang_sim<-ddply(data, .(Time), ang_comp,step=10)
ang_aver<-summarySE(ang_sim,measurevar="angles", groupvars=c("SHL"))
#plot.new()
# z <- locator(1)
ggplot(data=ang_aver,aes_string(x='SHL',y='angles'))+
geom_point(size=2)+
#geom_line(size=3)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$sd,ymax=ang_aver$angles+ang_aver$sd, color="blue",size=1)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$se,ymax=ang_aver$angles+ang_aver$se, color="black",size=1)+
xlab("SHL of angle vertex")+ylab('Angle,degrees')+geom_line(data=ang_cryst,aes_string(x='SHL',y='angles'),color="red",size=1)+
ggtitle('Angle between segments, 10 bp step')+
scale_x_continuous(breaks = round(seq(-7,7, by = 0.5),0.1))
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)
ggsave("../analysis_data/dna_bl_10.png",height=6,width=14)    


ang_cryst<-ang_comp(dfcryst,5)

ang_sim<-ddply(data, .(Time), ang_comp,step=5)
ang_aver<-summarySE(ang_sim,measurevar="angles", groupvars=c("SHL"))
#plot.new()
# z <- locator(1)
ggplot(data=ang_aver,aes_string(x='SHL',y='angles'))+
geom_point(size=2)+
#geom_line(size=3)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$sd,ymax=ang_aver$angles+ang_aver$sd, color="blue",size=1)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$se,ymax=ang_aver$angles+ang_aver$se, color="black",size=1)+
xlab("SHL of angle vertex")+ylab('Angle,degrees')+geom_line(data=ang_cryst,aes_string(x='SHL',y='angles'),color="red",size=1)+
ggtitle('Angle between segments, 5 bp step')+
scale_x_continuous(breaks = round(seq(-7,7, by = 0.5),0.1))
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)
ggsave("../analysis_data/dna_bl_5.png",height=6,width=14)    


######20
ang_cryst<-ang_comp(dfcryst,20)

ang_sim<-ddply(data, .(Time), ang_comp,step=20)
ang_aver<-summarySE(ang_sim,measurevar="angles", groupvars=c("SHL"))
#plot.new()
# z <- locator(1)
ggplot(data=ang_aver,aes_string(x='SHL',y='angles'))+
geom_point(size=2)+
#geom_line(size=3)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$sd,ymax=ang_aver$angles+ang_aver$sd, color="blue",size=1)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$se,ymax=ang_aver$angles+ang_aver$se, color="black",size=1)+
xlab("SHL of angle vertex")+ylab('Angle,degrees')+geom_line(data=ang_cryst,aes_string(x='SHL',y='angles'),color="red",size=1)+
ggtitle('Angle between segments, 20 bp step')+
scale_x_continuous(breaks = round(seq(-7,7, by = 0.5),0.1))
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)
ggsave("../analysis_data/dna_bl_20.png",height=6,width=14)   



ang_cryst<-ang_comp(dfcryst,1)

ang_sim<-ddply(data, .(Time), ang_comp,step=1)
ang_aver<-summarySE(ang_sim,measurevar="angles", groupvars=c("SHL"))
#plot.new()
# z <- locator(1)
ggplot(data=ang_aver,aes_string(x='SHL',y='angles'))+
geom_point(size=2)+
#geom_line(size=3)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$sd,ymax=ang_aver$angles+ang_aver$sd, color="blue",size=1)+
geom_errorbar(ymin=ang_aver$angles-ang_aver$se,ymax=ang_aver$angles+ang_aver$se, color="black",size=1)+
xlab("SHL of angle vertex")+ylab('Angle,degrees')+geom_line(data=ang_cryst,aes_string(x='SHL',y='angles'),color="red",size=1)+
ggtitle('Angle between segments, 1 bp step')+
scale_x_continuous(breaks = round(seq(-7,7, by = 0.5),0.1))
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)
ggsave("../analysis_data/dna_bl_1.png",height=6,width=14)  

# ggplot(data=ang_cryst,aes(x=SHL,y=angles))+
# xlab("SHL of angle vertex")+ylab("Angle, deg")+
# geom_line() + 
# scale_x_continuous(breaks = round(seq(-7,7, by = 0.5),0.1)) + ggtitle("Angle between segments, 5 bp step")
# ggsave("../analysis_data/dna_bl_cryst_5.png",height=6,width=14)    
# dev.off()
