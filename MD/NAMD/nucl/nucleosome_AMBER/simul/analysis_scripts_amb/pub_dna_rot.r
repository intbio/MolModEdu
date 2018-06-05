#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)

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

df_avr$DATA='Average'
# head(df_avr)
df_cryst$DATA='X-ray'
# head(df_cryst)
df=rbind(df_avr,df_cryst)
head(df)
df[df$Angle>0,'Angle']=df[df$Angle>0,'Angle']-180
a<-ggplot(data=df,aes(x=Basepair,y=-Angle,color=DATA))+
geom_line(size=1)+geom_point(size=2)+#ylim(0.2,4.0)+
#scale_fill_gradient2(limits=c(0,6),low="blue",mid="green", high="red",midpoint=3.0,guide_legend(title="RMSD, A"))+
#scale_y_continuous(limits=c(0,1210),breaks = round(seq(0,1000, by = 100),1))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
scale_color_manual(values=c('X-ray'="red",'Average'="blue"),name='',breaks=c('X-ray','Average'),labels=c('X-ray','MD average'))+
xlab("Base pair number")+ylab("Angle with superhelical axis, deg")+ggtitle("Nucleosomal DNA periodicity plot")+
# geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
# annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)
geom_text(data=seqdf,aes(x=X,y=-1.1,label=sequence),color='black',size=3)+
geom_text(data=seqdf2,aes(x=X,y=185,label=sequence),color='black',size=3)


ggsave("../analysis_data/pub_dna_rot.png",plot=a,height=4,width=14)
#plot.new()
# z <- locator(1)


