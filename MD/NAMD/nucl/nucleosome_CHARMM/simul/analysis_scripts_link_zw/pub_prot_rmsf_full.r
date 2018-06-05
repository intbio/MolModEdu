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

df<-df[,1:9]
dfm=melt(df,id.var=c("Resid"))
dfm=subset(dfm, value>0)
dfm=rename(dfm,c("value"="RMSF","variable"="Chain"))
head(dfm)

dfcm$data='X-ray'
dfm$data='MD'

df=rbind(dfm,dfcm)
head(df)

img <- readPNG(paste("seq_img/",'H3',".png",sep=''))
h3 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4',".png",sep=''))
h4 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2A',".png",sep=''))
h2a <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2B',".png",sep=''))
h2b <- rasterGrob(img, interpolate=TRUE,width=1)



# ggplot(data=df,aes(x=as.numeric(as.character(variable)),y=Time))+

a<-ggplot(data=df[df$Chain=='H3_1' ,],aes(x=Resid,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,10.0)+
#scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
#scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H3, A")+
geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)


e<-ggplot(data=df[df$Chain=='H3_2' ,],aes(x=Resid,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,10.0)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H3, E")+
geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)


b<-ggplot(data=df[df$Chain=='H4_1' ,],aes(x=Resid,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,10.0)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H4, B")+
geom_vline(xintercept = c(24,29,30,41,49,75,82,93), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h4, ymin=1.75, ymax=2.0, xmin=23.5,xmax=98.5)


f<-ggplot(data=df[df$Chain=='H4_2' ,],aes(x=Resid,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,2.0)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H4, F")+
geom_vline(xintercept = c(24,29,30,41,49,75,82,93), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h4, ymin=1.75, ymax=2.0, xmin=23.5,xmax=98.5)


c<-ggplot(data=df[df$Chain=='H2A_1' & df$Resid >=16 & df$Resid <=117,],aes(x=Resid,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,3.3)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H2A, C")+
geom_vline(xintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h2a, ymin=2.9, ymax=3.25, xmin=15.5,xmax=117.5)


g<-ggplot(data=df[df$Chain=='H2A_2' & df$Resid >=16 & df$Resid <=117,],aes(x=Resid,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,3.3)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H2A, G")+
geom_vline(xintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h2a, ymin=2.9, ymax=3.25, xmin=15.5,xmax=117.5)


d<-ggplot(data=df[df$Chain=='H2B_1' & df$Resid >=30 & df$Resid <=120,],aes(x=Resid+3,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,2.0)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H2B, D")+
geom_vline(xintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h2b, ymin=1.75, ymax=2.0, xmin=32.5,xmax=123.5)


h<-ggplot(data=df[df$Chain=='H2B_2' & df$Resid >=30 & df$Resid <=120,],aes(x=Resid+3,y=RMSF,color=data))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,2.0)+
xlab("Residue number")+ylab("RMSF, A")+ggtitle("RMSF, C-alpha, nucleosome core. Histone H2B, H")+
geom_vline(xintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
annotation_custom(h2b, ymin=1.75, ymax=2.0, xmin=32.5,xmax=123.5)



#facet_grid(Chain~.)
q<-arrangeGrob(a,b,e,f,c,d,g,h,ncol=2)

ggsave("../analysis_data/pub_prot_rmsf_full.png",plot=q,height=12,width=12)

#plot.new()
# z <- locator(1)
quit()
q<-arrangeGrob(t,s,ncol=1)

img <- readPNG(paste("dna_par_labels/",i,".png",sep=''))
g <- rasterGrob(img, interpolate=TRUE)
+ annotation_custom(g, ymin=0,7xmax=-st2e0+meand)
