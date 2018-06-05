#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)


#dfcryst$X<-as.factor(dfcryst$X)
df<-read.csv("../analysis_data/ss_nucl_avr.csv",skip=0,header=TRUE,check.name=FALSE)

# dfm=melt(df,id.var=c("Time"))

img <- readPNG(paste("seq_img/",'H3full_new',".png",sep=''))
h3 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4full_new',".png",sep=''))
h4 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Afull_new',".png",sep=''))
h2a <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Bfull_new',".png",sep=''))
h2b <- rasterGrob(img, interpolate=TRUE,width=1)


# ggplot(data=df,aes(x=as.numeric(as.character(variable)),y=Time))+
h3data=df[df$Chain=='CHA' | df$Chain=='CHE',]
h3data=ddply(h3data,c('Resid','ss_type'),summarize,ss_prop=mean(ss_prop))
h3data$ss_type=factor(h3data$ss_type, levels=c('B','H','G','C'))
h3data=arrange(h3data,ss_type)
head(h3data)

a<-ggplot(data=h3data,aes(x=Resid,y=(ss_prop)*(-100),fill=ss_type))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='stack')+scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('B','H','G','C'),values=c('B'='green','H'='orange','G'='black','C'='blue'),labels=c(expression(paste(beta,'-strand')),expression(paste(alpha,'-helix')),'3-10 helix','Coil/Turn'),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Propensity, %")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h3, ymin=0, ymax=50, xmin=0.5,xmax=135.5)
print(subset(h3data,Resid==119))
ggsave("../analysis_data/pub_prot_ss_h3.png",plot=a,height=2,width=12)


####H4

h4data=df[df$Chain=='CHB' | df$Chain=='CHF',]
h4data=ddply(h4data,c('Resid','ss_type'),summarize,ss_prop=mean(ss_prop))
h4data$ss_type=factor(h4data$ss_type, levels=c('B','H','G','C'))
h4data=arrange(h4data,ss_type)
head(h4data)

b<-ggplot(data=h4data,aes(x=Resid,y=(ss_prop)*(-100),fill=ss_type))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='stack')+scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,103),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('B','H','G','C'),values=c('B'='green','H'='orange','G'='black','C'='blue'),labels=c(expression(paste(beta,'-strand')),expression(paste(alpha,'-helix')),'3-10 helix','Coil/Turn'),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Propensity, %")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h4, ymin=0, ymax=50, xmin=0.5,xmax=102.5)

ggsave("../analysis_data/pub_prot_ss_h4.png",plot=b,height=2,width=10)


###H2A

h2adata=df[df$Chain=='CHC' | df$Chain=='CHG',]
h2adata=ddply(h2adata,c('Resid','ss_type'),summarize,ss_prop=mean(ss_prop))
h2adata$ss_type=factor(h2adata$ss_type, levels=c('B','H','G','C'))
h2adata=arrange(h2adata,ss_type)
head(h2adata)

c<-ggplot(data=h2adata,aes(x=Resid,y=(ss_prop)*(-100),fill=ss_type))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='stack')+scale_y_continuous(limits=c(-101,60),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,129),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('B','H','G','C'),values=c('B'='green','H'='orange','G'='black','C'='blue'),labels=c(expression(paste(beta,'-strand')),expression(paste(alpha,'-helix')),'3-10 helix','Coil/Turn'),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Propensity, %")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2a, ymin=0, ymax=60, xmin=0.5,xmax=128.5)

ggsave("../analysis_data/pub_prot_ss_h2a.png",plot=c,height=2,width=13)


####H2B
h2bdata=df[df$Chain=='CHD' | df$Chain=='CHH',]
h2bdata=ddply(h2bdata,c('Resid','ss_type'),summarize,ss_prop=mean(ss_prop))
h2bdata$ss_type=factor(h2bdata$ss_type, levels=c('B','H','G','C'))
h2bdata=arrange(h2bdata,ss_type)
head(h2bdata)

d<-ggplot(data=h2bdata,aes(x=Resid,y=(ss_prop)*(-100),fill=ss_type))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='stack')+scale_y_continuous(limits=c(-101,60),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,123),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('B','H','G','C'),values=c('B'='green','H'='orange','G'='black','C'='blue'),labels=c(expression(paste(beta,'-strand')),expression(paste(alpha,'-helix')),'3-10 helix','Coil/Turn'),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Propensity, %")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2b, ymin=0, ymax=60, xmin=0.5,xmax=122.5)

ggsave("../analysis_data/pub_prot_ss_h2b.png",plot=d,height=2,width=12)

q<-arrangeGrob(a,b,c,d,ncol=1)


ggsave("../analysis_data/pub_prot_ss.png",plot=q,height=8,width=12)

q()

e<-ggplot(data=df[df$Chain=='CHE',],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, E")+
annotation_custom(h3, ymin=1000, ymax=1270, xmin=43.5,xmax=131.5)


b<-ggplot(data=df[df$Chain=='CHB',],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1260),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H4, B")+
annotation_custom(h4, ymin=1000, ymax=1335, xmin=23.5,xmax=98.5)

f<-ggplot(data=df[df$Chain=='CHF',],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1260),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H4, F")+
annotation_custom(h4, ymin=1000, ymax=1335, xmin=23.5,xmax=98.5)



c<-ggplot(data=df[df$Chain=='CHC',],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1200),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H2A, C")+
annotation_custom(h2a, ymin=1000, ymax=1270, xmin=15.5,xmax=117.5)

g<-ggplot(data=df[df$Chain=='CHG',],aes(x=Resid,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1200),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H2A, G")+
annotation_custom(h2a, ymin=1000, ymax=1270, xmin=15.5,xmax=117.5)


d<-ggplot(data=df[df$Chain=='CHD',],aes(x=Resid+3,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1220),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H2B, D")+
annotation_custom(h2b, ymin=1000, ymax=1290, xmin=32.5,xmax=123.5)

h<-ggplot(data=df[df$Chain=='CHH',],aes(x=Resid+3,y=Time))+
geom_tile(aes(fill=RMSD)) + 
scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
scale_y_continuous(limits=c(0,1220),breaks = round(seq(0,1000, by = 100),1))+
xlab("Residue number")+ylab("Time, ns")+ggtitle("RMSD, C-alpha, nucleosome core. Histone H2B, H")+
annotation_custom(h2b, ymin=1000, ymax=1290, xmin=32.5,xmax=123.5)




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
