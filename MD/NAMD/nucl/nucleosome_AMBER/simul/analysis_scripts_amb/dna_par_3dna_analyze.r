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
####################
#Now we will make plot at overall parameter distributions

#make graph list
pl=list()
meand=c()
for (i in c("Shear","Stagger","Stretch","Buckle","Opening","Prop.Tw","Rise","Shift","Slide","Roll","Tilt","Twist")) {
# get rid of BP which are around those that are considered not paired by X3DNA, and also first BP
mask<-(df['BP']==1) & c(TRUE,df[1:nrow(df)-1,'BP'])
data<-df[!is.na(df[i]) & mask,i] 

#We assume that crystal has all base-pairs
datacryst<-dfcryst[!is.na(dfcryst[i]),i]

#Fit distrib params to MD simulation data
dp<-fitdist(data,"norm",method="mle")
str2<-paste(capture.output(dp), collapse = '\n')[1] 
meand<-as.numeric(dp$estimate[1]) 
stdev<-as.numeric(dp$estimate[2])
print(i)
#load pictures
img <- readPNG(paste("dna_par_labels/",i,".png",sep=''))
g <- rasterGrob(img, interpolate=TRUE)

#Do the actual plotting in object
pl[[i]]<-ggplot(df, aes_string(x=i)) + ggtitle(paste(i,"distribution")) + 
geom_density(aes_string(x=i),size=3,na.rm=TRUE,kernel=c("rectangular"),adjust=0.2) + geom_density(data=dfcryst, aes_string(x=i),size=3,color="red",na.rm=TRUE,kernel=c("rectangular"),adjust=0.2) +
#geom_histogram(aes_string(x=i,y="..density.."),size=3,na.rm=TRUE) +
#geom_histogram(data=dfcryst, aes_string(x=i,y="..density.."),size=3,color="red",na.rm=TRUE) +
geom_vline(aes_string(xintercept=meand), color="green", linetype="dashed", size=3) +
annotate("text", x=Inf, y=Inf, label = str2,hjust=1,vjust=1) + annotation_custom(g, ymin=0, xmax=-stdev+meand) + xlim(-3*stdev+meand,3*stdev+meand)
}

#collect graphs to one page and print them
png(file="../analysis_data/dna_param_distr.png",width=4000,height=2000)
multiplot(plotlist=pl,cols=4)
#ggsave("../analysis_data/dna_param_distr.png")
dev.off()


#Now let's analyze where BP become unpaired - plot the df$BP - column
theme_set(theme_grey())
#png(file="../analysis_data/dna_pairing.png",width=1000,height=800)
bp<-summarySE(df,measurevar="BP", groupvars=c("X"))

ggplot(data=bp[,],aes(X-74,BP*100))+
xlab("Base pair number")+ylab("Stability, %")+
geom_point() + geom_point(data=bp[bp$BP<1,],color="red") + 
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1)) + ggtitle("Stability of base pairing during simulation")
ggsave("../analysis_data/dna_pairing.png",height=6,width=14)
dev.off()

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


pl[[i]]<-ggplot(data=data,aes_string(x='X',y=i))+
geom_point(size=3)+
#geom_line(size=3)+
geom_errorbar(ymin=data[[i]]-data$sd,ymax=data[[i]]+data$sd, color="blue",size=3)+
geom_errorbar(ymin=data[[i]]-data$se,ymax=data[[i]]+data$se, color="black",size=3)+
xlab("Base pair / base pair step")+geom_line(data=dfcryst_t,aes_string(x='X',y=i),color="red",size=3)+
ggtitle(paste(i,"profile"))+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+#+xlim(-70,70)
geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)
}
png(file="../analysis_data/dna_param_profiles_sd.png",width=8000,height=4000)
multiplot(plotlist=pl,cols=2)
dev.off()

####Let's do the same with minmax
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
seqdf_t$Y<-min(data$min,dfcryst[[i]],na.rm=TRUE)


pl[[i]]<-ggplot(data=data,aes_string(x='X',y=i))+geom_point(size=3)+
geom_errorbar(ymin=data$min,ymax=data$max, color="blue",size=3)+
geom_errorbar(ymin=data[[i]]-data$se,ymax=data[[i]]+data$se, color="black",size=3)+
xlab("Base pair / base pair step")+geom_line(data=dfcryst_t,aes_string(x='X',y=i),color="red",size=3)+
ggtitle(paste(i,"profile"))+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+#+xlim(-70,70)
geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)
}
png(file="../analysis_data/dna_param_profiles_minmax.png",width=8000,height=4000)
multiplot(plotlist=pl,cols=2)
dev.off()

### Now let's make video
system('mkdir -p dat')
system('rm -f ../analysis_data/dna_param_profiles.mov')

#Let's get the time vector
time<-sort(unlist(unique(df['Time'])))
time<-time[seq(1,length(time),25)]
numf=0
for(t in time)
{
print(t)
####Let's loop to make a series of images
pl=list()
theme_set(theme_grey(base_size=40))

for (i in c("Shear","Stagger","Stretch","Buckle","Opening","Prop.Tw","Rise","Shift","Slide","Roll","Tilt","Twist")) {


if(i %in% c("Buckle","Opening","Prop.Tw","Stagger","Stretch","Shear")) {
shift=73.0
} else {
shift=73.5
}

data<-df[df$Time==t,]
data<-data[4:143,]
data$X=data$X-shift
dfcryst_t<-dfcryst
dfcryst_t$X=dfcryst_t$X-shift
seqdf_t<-seqdf
seqdf_t$X<-seqdf_t$X-73.0
seqdf_t$Y<-min(dfcryst[[i]],data[[i]],na.rm=TRUE)


pl[[i]]<-ggplot(data=data,aes_string(x='X',y=i))+geom_line(size=3)+
#geom_errorbar(ymin=data[[i]]-data$sd,ymax=data[[i]]+data$sd, color="blue",size=3)+
#geom_errorbar(ymin=data[[i]]-data$se,ymax=data[[i]]+data$se, color="black",size=3)+
xlab("Base pair / base pair step")+geom_line(data=dfcryst_t,aes_string(x='X',y=i),color="red",size=3)+
ggtitle(paste(i,"profile"))+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+#+xlim(-70,70)
geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)

}
png(file=paste("dat/",numf,".png",sep=""),width=8000,height=4000)
multiplot(plotlist=pl,cols=2)
dev.off()
numf<-numf+1
}

system('ffmpeg -r 5 -i dat/%d.png -s 2000x1000 -q:v 0 -pix_fmt yuv420p  ../analysis_data/dna_param_profiles.mov')

system('rm -r dat')

#plot.new()
# z <- locator(1)


