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
df<-read.table("../analysis_data/pub_prot_rmsd.csv",skip=0,header=TRUE,check.name=FALSE)
dfsch<-read.table("../analysis_data/pub_prot_rmsd_sch.csv",skip=0,header=TRUE,check.name=FALSE)

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
h3data=ddply(h3data,c('Resid'),summarize,max_rmsd=max(RMSD))
h3data$color='CA'
h3data[h3data$max_rmsd>6,'color']='amore'
h3data[h3data$max_rmsd>6,'max_rmsd']=6

h3data_sch=dfsch[dfsch$Chain=='CHA' | dfsch$Chain=='CHE',]
h3data_sch=ddply(h3data_sch,c('Resid'),summarize,max_rmsd=max(RMSD))
h3data_sch$color='bSCH'
h3data_sch[h3data_sch$max_rmsd>6,'color']='amore'
h3data_sch[h3data_sch$max_rmsd>6,'max_rmsd']=6
h3data=rbind(h3data,h3data_sch)
h3data=arrange(h3data,color,Resid)
head(h3data)

a<-ggplot(data=h3data,aes(x=Resid,y=(max_rmsd)*(-1),fill=color))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='identity')+scale_y_continuous(limits=c(-6,0.5),breaks=seq(-6,0),labels=seq(6,0,by=-1))+
scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("RMSD, Å")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h3, ymin=-1, ymax=0, xmin=0.5,xmax=135.5)

ggsave("../analysis_data/pub_prot_max_rmsd_h3.png",plot=a,height=2,width=12)

####H4

h4data=df[df$Chain=='CHB' | df$Chain=='CHF',]
h4data=ddply(h4data,c('Resid'),summarize,max_rmsd=max(RMSD))
h4data$color='CA'
h4data[h4data$max_rmsd>6,'color']='amore'
h4data[h4data$max_rmsd>6,'max_rmsd']=6

h4data_sch=dfsch[dfsch$Chain=='CHB' | dfsch$Chain=='CHF',]
h4data_sch=ddply(h4data_sch,c('Resid'),summarize,max_rmsd=max(RMSD))
h4data_sch$color='bSCH'
h4data_sch[h4data_sch$max_rmsd>6,'color']='amore'
h4data_sch[h4data_sch$max_rmsd>6,'max_rmsd']=6
h4data=rbind(h4data,h4data_sch)
h4data=arrange(h4data,color,Resid)
head(h4data)

a<-ggplot(data=h4data,aes(x=Resid,y=(max_rmsd)*(-1),fill=color))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='identity')+scale_y_continuous(limits=c(-6,0.5),breaks=seq(-6,0),labels=seq(6,0,by=-1))+
scale_x_continuous(limits=c(0,103),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("RMSD, Å")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h4, ymin=-1, ymax=0, xmin=0.5,xmax=102.5)

ggsave("../analysis_data/pub_prot_max_rmsd_h4.png",plot=a,height=2,width=10)

###H2A

h2adata=df[df$Chain=='CHC' | df$Chain=='CHG',]
h2adata=ddply(h2adata,c('Resid'),summarize,max_rmsd=max(RMSD))
h2adata$color='CA'
h2adata[h2adata$max_rmsd>6,'color']='amore'
h2adata[h2adata$max_rmsd>6,'max_rmsd']=6

h2adata_sch=dfsch[dfsch$Chain=='CHC' | dfsch$Chain=='CHG',]
h2adata_sch=ddply(h2adata_sch,c('Resid'),summarize,max_rmsd=max(RMSD))
h2adata_sch$color='bSCH'
h2adata_sch[h2adata_sch$max_rmsd>6,'color']='amore'
h2adata_sch[h2adata_sch$max_rmsd>6,'max_rmsd']=6
h2adata=rbind(h2adata,h2adata_sch)
h2adata=arrange(h2adata,color,Resid)
head(h2adata)

a<-ggplot(data=h2adata,aes(x=Resid,y=(max_rmsd)*(-1),fill=color))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='identity')+scale_y_continuous(limits=c(-6,0.5),breaks=seq(-6,0),labels=seq(6,0,by=-1))+
scale_x_continuous(limits=c(0,129),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("RMSD, Å")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2a, ymin=-1, ymax=0, xmin=0.5,xmax=128.5)

ggsave("../analysis_data/pub_prot_max_rmsd_h2a.png",plot=a,height=2,width=13)

####H2B
h2bdata=df[df$Chain=='CHD' | df$Chain=='CHH',]
h2bdata=ddply(h2bdata,c('Resid'),summarize,max_rmsd=max(RMSD))
h2bdata$color='CA'
h2bdata[h2bdata$max_rmsd>6,'color']='amore'
h2bdata[h2bdata$max_rmsd>6,'max_rmsd']=6

h2bdata_sch=dfsch[dfsch$Chain=='CHD' | dfsch$Chain=='CHH',]
h2bdata_sch=ddply(h2bdata_sch,c('Resid'),summarize,max_rmsd=max(RMSD))
h2bdata_sch$color='bSCH'
h2bdata_sch[h2bdata_sch$max_rmsd>6,'color']='amore'
h2bdata_sch[h2bdata_sch$max_rmsd>6,'max_rmsd']=6
h2bdata=rbind(h2bdata,h2bdata_sch)
h2bdata=arrange(h2bdata,color,Resid)
head(h2bdata)

a<-ggplot(data=h2bdata,aes(x=Resid,y=(max_rmsd)*(-1),fill=color))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(stat='identity',position='identity')+scale_y_continuous(limits=c(-6,0.5),breaks=seq(-6,0),labels=seq(6,0,by=-1))+
scale_x_continuous(limits=c(0,123),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("RMSD, Å")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2b, ymin=-1, ymax=0, xmin=0.5,xmax=122.5)

ggsave("../analysis_data/pub_prot_max_rmsd_h2b.png",plot=a,height=2,width=12)


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
