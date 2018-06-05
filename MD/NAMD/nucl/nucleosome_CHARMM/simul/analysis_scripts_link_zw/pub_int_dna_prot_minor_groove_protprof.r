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
dna_prot_avr<-read.csv('../analysis_data/dna_prot_avr_df.csv')


dna_prot<-subset(dna_prot_avr,type=='SC' & groove=='minor' & PROT_part=='side chain')
# dfm=melt(df,id.var=c("Time"))
head(dna_prot)
# q()
img <- readPNG(paste("seq_img/",'H3_full_inv',".png",sep=''))
h3 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4_full_inv',".png",sep=''))
h4 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2A_full_inv',".png",sep=''))
h2a <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2B_full_inv',".png",sep=''))
h2b <- rasterGrob(img, interpolate=TRUE,width=1)

theme_set(theme_bw()+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))

dna_prot$PROT_chain=revalue(dna_prot$PROT_chain, c("CHA"="A","CHB"="B","CHC"="C","CHD"="D","CHE"="E","CHF"="F","CHG"="G","CHH"="H"))
# ggplot(data=df,aes(x=as.numeric(as.character(variable)),y=Time))+
h3data=dna_prot[dna_prot$PROT_chain=='A' | dna_prot$PROT_chain=='E',]
# h3data=ddply(h3data,c('Resid'),summarize,sasa_fract=mean(sasa_fract))
# h3data$ss_type=factor(h3data$ss_type, levels=c('B','H','G','C'))
h3data=ddply(h3data,c('PROT_chain','PROT_resid'),summarize,tot_num=sum(av_num))

h3dataE=dcast(h3data,PROT_resid~PROT_chain)
h3dataE=subset(h3dataE,is.na(A))
h3dataE=melt(h3dataE,id.vars=c('PROT_resid'),measure.vars=c('E'),variable.name='PROT_chain',value.name='tot_num')

h3data=arrange(h3data,PROT_chain,PROT_resid,desc(tot_num))
head(h3data)

a<-ggplot(data=h3data,aes(x=PROT_resid,y=(tot_num),fill=PROT_chain,color=PROT_chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(data=subset(h3data,(PROT_chain=='A')),stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_bar(data=h3dataE,stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-15,15),breaks=seq(0,10,by=2.0),labels=seq(0,10,by=2.0),expand=c(0,0))+
scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('A','E'),labels=c('A','E'),values=c('A'='indianred1','E'='green4'),name='Chain')+
scale_color_manual(breaks=c('A','E'),labels=c('A','E'),values=c('A'='indianred1','E'='green4'),name='Chain',guide=FALSE)+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+guides(colour = guide_legend(override.aes = list(shape = NA)))+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Number")+ggtitle("H3, contacts of side chains with bases in minor grooves")+#geom_point(data=data.frame(PROT_resid=100,tot_num=1),color='red',fill='red')+
annotation_custom(h3, ymin=-15, ymax=0, xmin=0.5,xmax=135.5)

ggsave("../analysis_data/pub_int_dna_prot_minor_groove_h3.png",plot=a,height=2.5,width=12)

####H4

h4data=dna_prot[dna_prot$PROT_chain=='B' | dna_prot$PROT_chain=='F',]
h4data=ddply(h4data,c('PROT_chain','PROT_resid'),summarize,tot_num=sum(av_num))

h4dataF=dcast(h4data,PROT_resid~PROT_chain)
h4dataF=subset(h4dataF,is.na(B))
h4dataF=melt(h4dataF,id.vars=c('PROT_resid'),measure.vars=c('F'),variable.name='PROT_chain',value.name='tot_num')

h4data=arrange(h4data,PROT_chain,PROT_resid,desc(tot_num))
head(h4data)

a<-ggplot(data=h4data,aes(x=PROT_resid,y=(tot_num),fill=PROT_chain,color=PROT_chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(data=subset(h4data,(PROT_chain=='B')),stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_bar(data=h4dataF,stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-15,8),breaks=seq(0,10,by=2.0),labels=seq(0,10,by=2.0),expand=c(0,0))+
scale_x_continuous(limits=c(0,103),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('B','F'),labels=c('B','F'),values=c('B'='indianred1','F'='green4'),name='Chain')+
scale_color_manual(breaks=c('B','F'),labels=c('B','F'),values=c('B'='indianred1','F'='green4'),name='Chain',guide=FALSE)+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+guides(colour = guide_legend(override.aes = list(shape = NA)))+#geom_point(data=data.frame(PROT_resid=100,tot_num=1),color='red',fill='red')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Number")+ggtitle("H4, contacts of side chains with bases in minor grooves")+
annotation_custom(h4, ymin=-15, ymax=0, xmin=0.5,xmax=102.5)

ggsave("../analysis_data/pub_int_dna_prot_minor_groove_h4.png",plot=a,height=2.5,width=12)

###H2A


h2adata=dna_prot[dna_prot$PROT_chain=='C' | dna_prot$PROT_chain=='G',]
h2adata=ddply(h2adata,c('PROT_chain','PROT_resid'),summarize,tot_num=sum(av_num))

h2adataF=dcast(h2adata,PROT_resid~PROT_chain)
h2adataF=subset(h2adataF,is.na(C))
h2adataF=melt(h2adataF,id.vars=c('PROT_resid'),measure.vars=c('G'),variable.name='PROT_chain',value.name='tot_num')

h2adata=arrange(h2adata,PROT_chain,PROT_resid,desc(tot_num))
head(h2adata)

a<-ggplot(data=h2adata,aes(x=PROT_resid,y=(tot_num),fill=PROT_chain,color=PROT_chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(data=subset(h2adata,(PROT_chain=='C')),stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_bar(data=h2adataF,stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-15,8),breaks=seq(0,10,by=2.0),labels=seq(0,10,by=2.0),expand=c(0,0))+
scale_x_continuous(limits=c(0,129),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('C','G'),labels=c('C','G'),values=c('C'='indianred1','G'='green4'),name='Chain')+
scale_color_manual(breaks=c('C','G'),labels=c('C','G'),values=c('C'='indianred1','G'='green4'),name='Chain',guide=FALSE)+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+guides(colour = guide_legend(override.aes = list(shape = NA)))+#geom_point(data=data.frame(PROT_resid=100,tot_num=1),color='red',fill='red')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Number")+ggtitle("H2A, contacts of side chains with bases in minor grooves")+
annotation_custom(h2a, ymin=-15, ymax=0, xmin=0.5,xmax=128.5)


ggsave("../analysis_data/pub_int_dna_prot_minor_groove_h2a.png",plot=a,height=2.2,width=12)



###H2B

h2bdata=dna_prot[dna_prot$PROT_chain=='D' | dna_prot$PROT_chain=='H',]
h2bdata=ddply(h2bdata,c('PROT_chain','PROT_resid'),summarize,tot_num=sum(av_num))

h2bdataH=dcast(h2bdata,PROT_resid~PROT_chain)
h2bdataH=subset(h2bdataH,is.na(D))
h2bdataH=melt(h2bdataH,id.vars=c('PROT_resid'),measure.vars=c('H'),variable.name='PROT_chain',value.name='tot_num')

h2bdata=arrange(h2bdata,PROT_chain,PROT_resid,desc(tot_num))
head(h2bdata)

a<-ggplot(data=h2bdata,aes(x=PROT_resid+3,y=(tot_num),fill=PROT_chain,color=PROT_chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(data=subset(h2bdata,(PROT_chain=='D')),stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_bar(data=h2bdataH,stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-15,8),breaks=seq(0,10,by=2.0),labels=seq(0,10,by=2.0),expand=c(0,0))+
scale_x_continuous(limits=c(0,126),labels=c(),breaks=c(),expand=c(0,0))+
scale_fill_manual(breaks=c('D','H'),labels=c('D','H'),values=c('D'='indianred1','H'='green4'),name='Chain')+
scale_color_manual(breaks=c('D','H'),labels=c('D','H'),values=c('D'='indianred1','H'='green4'),name='Chain',guide=FALSE)+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+guides(colour = guide_legend(override.aes = list(shape = NA)))+#geom_point(data=data.frame(PROT_resid=117,tot_num=5),color='red',fill='red')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Number")+ggtitle("H2B, contacts of side chains with bases in minor grooves")+
annotation_custom(h2b, ymin=-15, ymax=0, xmin=0.5,xmax=125.5)


ggsave("../analysis_data/pub_int_dna_prot_minor_groove_h2b.png",plot=a,height=2.2,width=12)




q()



####H2B
h2bdata=df[df$Chain=='D' | df$Chain=='H',]
h2bdata=arrange(h2bdata,Chain,Resid,desc(sasa_fract))

d<-ggplot(data=h2bdata,aes(x=Resid,y=(sasa_fract)*(-100),fill=Chain,color=Chain))+
# geom_tile(aes(fill=RMSD)) + 
geom_bar(data=subset(h2bdata,Chain=='D'),stat='identity',position='dodge',width=0.7)+#scale_y_continuous(limits=c(-101,50),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
geom_point()+scale_y_continuous(limits=c(-101,60),breaks=seq(-100,0,by=20),labels=seq(100,0,by=-20),expand=c(0,0))+
scale_x_continuous(limits=c(0,123),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
# scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
xlab('')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("Relative SASA, %")+#ggtitle("RMSD, C-alpha, nucleosome core. Histone H3, A")+
annotation_custom(h2b, ymin=0, ymax=60, xmin=0.5,xmax=122.5)

ggsave("../analysis_data/pub_prot_sasa_h2b.png",plot=d,height=2,width=12)

q<-arrangeGrob(a,b,c,d,ncol=1)


ggsave("../analysis_data/pub_prot_sasa.png",plot=q,height=8,width=12)



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
