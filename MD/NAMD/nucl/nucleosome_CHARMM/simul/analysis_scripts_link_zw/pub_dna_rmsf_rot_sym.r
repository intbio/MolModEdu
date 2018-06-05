#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)



# df_cryst<-read.table("1kx5_R.pdb",skip=0,header=TRUE,check.name=FALSE)
# df_cryst<-df_cryst[(df_cryst$Aname %in% c('CA','P')),c('Chname','Resid','Bfactor')]
# chname_conv=data.frame(Chname=c('A','B','C','D','E','F','G','H','I','J'),Chain=c('H3_1','H4_1','H2A_1','H2B_1','H3_2','H4_2','H2A_2','H2B_2','DNA_I','DNA_J'))
# dfc=merge(df_cryst,chname_conv)[,c('Chain','Resid','Bfactor')]
# dfc$RMSF=sqrt(dfc$Bfactor/8/3.1415/3.1415)
# dfcm=dfc[,c('Chain','Resid','RMSF')]

# head(dfcm)
#dfcryst$X<-as.factor(dfcryst$X)
df<-read.table("../analysis_data/rmsf_chains.dat",skip=4,header=TRUE,check.name=FALSE)

df<-df[,c('Resid','DNA_I','DNA_J')]
dfm=melt(df,id.var=c("Resid"))
dfm=subset(dfm, value>0)
dfm=rename(dfm,c("value"="RMSF","variable"="Chain"))
head(dfm)

# dfcm$data='X-ray'
# dfm$data='MD'

# df=rbind(dfm,dfcm)
df=dfm
df=df[(df$Resid>-74)&(df$Resid<75),]
head(df)

# ggplot(data=df,aes(x=as.numeric(as.character(variable)),y=Time))+

df[(df$Chain=='DNA_J'),'Resid']=(-1)*df[(df$Chain=='DNA_J'),'Resid']
dfp=subset(df,Resid>=0)
dfm=subset(df,Resid<=0)
dfm$Resid=(-1)*dfm$Resid
dfp$turn='positive'
dfm$turn='negative'
df=rbind(dfp,dfm)


theme_set(theme_grey(base_size = 18))


a<-ggplot(data=df,aes(x=Resid/10,y=RMSF,color=Chain,alpha=turn))+
geom_line(size=1)+geom_point(size=2)+ylim(0.2,4.0)+
#scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
#scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
scale_color_manual(values=c("red","blue"),name='',breaks=c('DNA_I','DNA_J'),labels=c('Chain I','Chain J'))+
ylab("RMSF, Å")+#ggtitle("RMSF, P-atoms, core nucleosome DNA, FN-model.")+
# geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
# annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)
scale_x_continuous(limits=c(-0.05,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+
scale_alpha_discrete(range = c(0.3, 1.0),name='SHL')


ggsave("../analysis_data/pub_dna_rmsf_sym.png",plot=a,width=15,height=2.5)





dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))

dnaseq2<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
dnaseq2<-gsub("([AG])", "R", dnaseq2)
dnaseq2<-gsub("([CT])", "Y", dnaseq2)

seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))
seqdf2=data.frame(sequence=substring(dnaseq2, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))


#dfcryst$X<-as.factor(dfcryst$X)
df_avr<-read.csv("../analysis_data/dna_rot_df_avr.csv",header=TRUE,check.name=FALSE)
df_cryst<-read.csv("../analysis_data/dna_rot_df_cryst.csv",header=TRUE,check.name=FALSE)

df_avr$DATA='MD'
# head(df_avr)
df_cryst$DATA='X-ray'
# head(df_cryst)
df=rbind(df_avr,df_cryst)
# df=df_avr
head(df)

dfp=subset(df,Basepair>=0)
dfm=subset(df,Basepair<=0)
dfm$Basepair=(-1)*dfm$Basepair
dfp$turn='positive'
dfm$turn='negative'
df=rbind(dfp,dfm)

df[df$Angle>0,'Angle']=df[df$Angle>0,'Angle']-180

a<-ggplot(data=df,aes(x=Basepair/10,y=-Angle,alpha=turn,color=DATA,size=DATA))+
geom_line()+geom_point(size=2)+#ylim(0.2,4.0)+
scale_size_discrete(range=c(1,0.5),guide=FALSE)+
#scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
#scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
scale_color_manual(values=c('X-ray'="red",'MD'="black"),name='',breaks=c('X-ray','MD'),labels=c('X-ray','MD'))+
xlab("Base pair number")+ylab("Angle, deg")+ggtitle("Nucleosomal DNA periodicity plot")+
# geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
# annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)
geom_text(data=seqdf,aes(x=X/10,y=-1.5,label=sequence),alpha=1,color='black',size=3)+
geom_text(data=seqdf2,aes(x=X/10,y=185,label=sequence),alpha=1,color='black',size=3)+
scale_x_continuous(limits=c(0,7.5),breaks = round(seq(0,9.0, by = 0.5),2),labels=c('0','±0.5','±1.0','±1.5','±2.0','±2.5','±3.0','±3.5','±4.0','±4.5','±5.0','±5.5','±6.0','±6.5','±7.0','±7.5','±8.0','±8.5','±9.0'),expand=c(0,0),name='Super helix location (SHL)')+scale_alpha_discrete(range = c(0.3, 1.0),name='SHL')



ggsave("../analysis_data/pub_dna_rot_sym.png",plot=a,height=2.5,width=15)




#plot.new()
# z <- locator(1)
quit()
q<-arrangeGrob(t,s,ncol=1)

img <- readPNG(paste("dna_par_labels/",i,".png",sep=''))
g <- rasterGrob(img, interpolate=TRUE)
+ annotation_custom(g, ymin=0,7xmax=-st2e0+meand)
