#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)



df_cryst<-read.table("1kx5_R.pdb",skip=0,header=TRUE,check.name=FALSE)
df_cryst<-df_cryst[(df_cryst$Aname %in% c('CA','P')),c('Chname','Resid','Bfactor')]
chname_conv=data.frame(Chname=c('A','B','C','D','E','F','G','H','I','J'),Chain=c('H3_1','H4_1','H2A_1','H2B_1','H3_2','H4_2','H2A_2','H2B_2','DNA_I','DNA_J'))
dfc=merge(df_cryst,chname_conv)[,c('Chain','Resid','Bfactor')]
dfc$RMSF=sqrt(dfc$Bfactor/8/3.1415/3.1415)
dfcm=dfc[,c('Chain','Resid','RMSF')]

head(dfcm)
#dfcryst$X<-as.factor(dfcryst$X)
df<-read.table("../analysis_data/rmsf_chains.dat",skip=4,header=TRUE,check.name=FALSE)

df<-df[,c('Resid','DNA_I','DNA_J')]
dfm=melt(df,id.var=c("Resid"))
dfm=subset(dfm, value>0)
dfm=rename(dfm,c("value"="RMSF","variable"="Chain"))
head(dfm)

dfcm$data='X-ray'
dfm$data='MD'

df=rbind(dfm,dfcm)
head(df)


# ggplot(data=df,aes(x=as.numeric(as.character(variable)),y=Time))+

a<-ggplot(data=df[df$Chain %in% c('DNA_I','DNA_J') ,],aes(x=Resid,y=RMSF,linetype=data,color=Chain))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,5.0)+
#scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
#scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
scale_color_manual(values=c("red","blue"),name='',breaks=c('DNA_I','DNA_J'),labels=c('Chain I','Chain J'))+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, P-atoms, nucleosome DNA.")
# geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
# annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)


ggsave("../analysis_data/pub_dna_rmsf.png",plot=a,height=6,width=12)

#plot.new()
# z <- locator(1)
quit()
q<-arrangeGrob(t,s,ncol=1)

img <- readPNG(paste("dna_par_labels/",i,".png",sep=''))
g <- rasterGrob(img, interpolate=TRUE)
+ annotation_custom(g, ymin=0,7xmax=-st2e0+meand)
