#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)


#dfcryst$X<-as.factor(dfcryst$X)
df<-read.table("../analysis_data/pub_prot_rmsd.csv",skip=0,header=TRUE,check.name=FALSE)

# dfm=melt(df,id.var=c("Time"))

img <- readPNG(paste("seq_img/",'H3core_new',".png",sep=''))
h3 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4core_new',".png",sep=''))
h4 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Acore_new',".png",sep=''))
h2a <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Bcore_new',".png",sep=''))
h2b <- rasterGrob(img, interpolate=TRUE,width=1)


# ggplot(data=df,aes(x=as.numeric(as.character(variable)),y=Time))+

a<-ggplot(data=df[df$Chain=='CHA'&df$Resid>43,],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H3, A")+
annotation_custom(h3, ymin=1000, ymax=1270, xmin=43.5,xmax=135.5)

e<-ggplot(data=df[df$Chain=='CHE'&df$Resid>43,],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H3, E")+
annotation_custom(h3, ymin=1000, ymax=1270, xmin=43.5,xmax=135.5)


b<-ggplot(data=df[df$Chain=='CHB'&df$Resid>33,],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1260),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H4, B")+
annotation_custom(h4, ymin=1000, ymax=1335, xmin=33.5,xmax=102.5)

f<-ggplot(data=df[df$Chain=='CHF'&df$Resid>33,],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1260),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H4, F")+
annotation_custom(h4, ymin=1000, ymax=1335, xmin=33.5,xmax=102.5)



c<-ggplot(data=df[df$Chain=='CHC'&df$Resid>15&df$Resid<119,],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1200),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H2A, C")+
annotation_custom(h2a, ymin=1000, ymax=1270, xmin=15.5,xmax=118.5)

g<-ggplot(data=df[df$Chain=='CHG'&df$Resid>15&df$Resid<119,],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1200),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H2A, G")+
annotation_custom(h2a, ymin=1000, ymax=1270, xmin=15.5,xmax=118.5)


d<-ggplot(data=df[df$Chain=='CHD'&df$Resid>29,],aes(x=Resid+3,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1220),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H2B, D")+
annotation_custom(h2b, ymin=1000, ymax=1290, xmin=32.5,xmax=125.5)

h<-ggplot(data=df[df$Chain=='CHH'&df$Resid>29,],aes(x=Resid+3,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"),na.value='orange')+
scale_y_continuous(limits=c(0,1220),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, NCPamb-model. Histone H2B, H")+
annotation_custom(h2b, ymin=1000, ymax=1290, xmin=32.5,xmax=125.5)




#facet_grid(Chain~.)
q<-arrangeGrob(a,b,e,f,c,d,g,h,ncol=2)


ggsave("../analysis_data/pub_prot_rmsd.png",plot=q,height=12,width=12)
#plot.new()
# z <- locator(1)
quit()
q<-arrangeGrob(t,s,ncol=1)

img <- readPNG(paste("dna_par_labels/",i,".png",sep=''))
g <- rasterGrob(img, interpolate=TRUE)
+ annotation_custom(g, ymin=0, xmax=-stdev+meand)
