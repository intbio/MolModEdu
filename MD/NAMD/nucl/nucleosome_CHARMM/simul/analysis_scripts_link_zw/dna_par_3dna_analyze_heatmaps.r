#R-script to analyze DNA-param data
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape)
##

#####
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
#seqnumbering starts from 0 here
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
#Get crystal and simulation data
#dfcryst<-read.csv("../analysis_data/dna_param_cryst.csv")
#dfcryst$X<-as.factor(dfcryst$X)
df<-read.csv("../analysis_data/dna_param_big_df.csv")
#df$X<-as.factor(df$X)
#Let's add a sequence column
####################
#Now we will make plot at overall parameter distributions



seqdf_t<-seqdf
seqdf_t$X<-seqdf_t$X-73.0
seqdf_t$Y<-0



ggplot(data=df[df$X>9 & df$X<137,],aes(x=X-73.5,y=Time/10))+
geom_tile(aes(fill=Roll))+scale_fill_gradient2(low="blue",mid="grey", high="red",midpoint=0.0,guide_legend(title="Roll, deg"))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
xlab("Base pair")+ylab("Time, ns")+ggtitle("Roll evolution")+
geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=1)

ggsave("../analysis_data/dna_heatmap_roll.png",height=8,width=8)

ggplot(data=df[df$X>9 & df$X<137,],aes(x=X-73.5,y=Time/10))+
geom_tile(aes(fill=Slide))+scale_fill_gradient2(low="blue",mid="grey", high="red",midpoint=0.0,guide_legend(title="Slide, A"))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
xlab("Base pair")+ylab("Time, ns")+ggtitle("Slide evolution")+
geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=1)

ggsave("../analysis_data/dna_heatmap_slide.png",height=8,width=8)

#plot.new()
# z <- locator(1)


