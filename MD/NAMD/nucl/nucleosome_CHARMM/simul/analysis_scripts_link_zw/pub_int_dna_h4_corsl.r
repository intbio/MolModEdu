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
dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')


####H4


dna_prot_h3<-subset(dna_prot,type=='SC' & (PROT_chain%in%c('CHF','CHB')) & (PROT_resid%in%c(8,12,15,16,18,20)))
# dfm=melt(df,id.var=c("Time"))

dna_prot_h3=ddply(dna_prot_h3,c('DNA_chain', 'DNA_resid', 'PROT_chain', 'PROT_resid'),summarize,num=length(param1)/1000)
dna_prot_h3[dna_prot_h3$DNA_chain=='CHJ','DNA_resid']=-dna_prot_h3[dna_prot_h3$DNA_chain=='CHJ','DNA_resid']
head(dna_prot_h3)
dna_prot_h3$PROT_resid=factor(dna_prot_h3$PROT_resid)

theme_set(theme_bw()+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=dna_prot_h3,aes(x=PROT_resid,y=(DNA_resid/10),color=PROT_resid,size=num))+
# geom_tile(aes(fill=RMSD)) + 
geom_point()+
xlab('Time')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("DNA SHL")+ggtitle("H4, contacts of selected residues with DNA")+facet_grid(PROT_chain~.,scales='free')

ggsave("../analysis_data/pub_int_dna_h4_crosl.png",plot=a,height=5,width=12)


q()

dna_prot_h3<-subset(dna_prot,type=='SC' & (PROT_chain%in%c('CHG','CHC')) & (PROT_resid%in%c(3,9,11,13)) & groove=='minor')
# dfm=melt(df,id.var=c("Time"))

dna_prot_h3=ddply(dna_prot_h3,c('DNA_chain', 'DNA_resid', 'PROT_chain', 'PROT_resid', 'Time'),summarize,num=length(param1))
dna_prot_h3[dna_prot_h3$DNA_chain=='CHJ','DNA_resid']=-dna_prot_h3[dna_prot_h3$DNA_chain=='CHJ','DNA_resid']
head(dna_prot_h3)
dna_prot_h3$PROT_resid=factor(dna_prot_h3$PROT_resid)

theme_set(theme_bw()+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=dna_prot_h3,aes(x=Time,y=(DNA_resid/10),color=PROT_resid))+
# geom_tile(aes(fill=RMSD)) + 
geom_point()+
xlab('Time')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("DNA SHL")+ggtitle("H2A, contacts of anchor residues with bases in minor grooves")+facet_grid(PROT_chain~.,scales='free')

ggsave("../analysis_data/pub_int_dna_anch_dyn_h2a.png",plot=a,height=5,width=12)




dna_prot_h3<-subset(dna_prot,type=='SC' & (PROT_chain%in%c('CHH','CHD')) & (PROT_resid%in%c(5-3,29-3,30-3,33-3)) & groove=='minor')
# dfm=melt(df,id.var=c("Time"))
dna_prot_h3$PROT_resid=dna_prot_h3$PROT_resid+3
dna_prot_h3=ddply(dna_prot_h3,c('DNA_chain', 'DNA_resid', 'PROT_chain', 'PROT_resid', 'Time'),summarize,num=length(param1))
dna_prot_h3[dna_prot_h3$DNA_chain=='CHJ','DNA_resid']=-dna_prot_h3[dna_prot_h3$DNA_chain=='CHJ','DNA_resid']
head(dna_prot_h3)
dna_prot_h3$PROT_resid=factor(dna_prot_h3$PROT_resid)

theme_set(theme_bw()+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=dna_prot_h3,aes(x=Time,y=(DNA_resid/10),color=PROT_resid))+
# geom_tile(aes(fill=RMSD)) + 
geom_point()+
xlab('Time')+
# theme(plot.margin = unit(c(1,1,1,1), "cm"))+
ylab("DNA SHL")+ggtitle("H2B, contacts of anchor residues with bases in minor grooves")+facet_grid(PROT_chain~.,scales='free')

ggsave("../analysis_data/pub_int_dna_anch_dyn_h2b.png",plot=a,height=5,width=12)


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
