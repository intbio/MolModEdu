#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)

#dfcryst$X<-as.factor(dfcryst$X)
df<-read.table("../analysis_data/dna_rmsd.dat",skip=4,header=TRUE,check.name=FALSE)

dfm=melt(df,id.var=c("Time"))

ggplot(data=dfm,aes(x=as.numeric(as.character(variable)),y=Time))+
geom_tile(aes(fill=value)) + 
scale_fill_gradient2(low="blue",mid="green", high="red",midpoint=7.0,guide_legend(title="RMSD, A"))+
scale_x_continuous(breaks = round(seq(-70,70, by = 10),1))+
xlab("Base pair")+ylab("Time, ns")+ggtitle("Base pair RMSD relative to crystal structure")

ggsave("../analysis_data/dna_rmsd.png",height=6,width=6)
#plot.new()
# z <- locator(1)


