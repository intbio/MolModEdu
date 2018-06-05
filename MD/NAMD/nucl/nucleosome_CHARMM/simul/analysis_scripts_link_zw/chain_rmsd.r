#R-script to diplay DNA hbonds
library(matrixStats)
png(file="..//analysis_data//chain_rmsd.png",width=10,height=5,units="in",res=600)
df<-read.table("../analysis_data//chain_rmsd.dat",skip=4,header=TRUE)
rmsd<-as.matrix(df[,-1])
min<-colMins(rmsd)
max<-colMaxs(rmsd)
mean<-colMeans(rmsd)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#ratio<-mean/max
data<-rbind(min,mean-min,max-mean)
barx <- barplot(data, beside=FALSE, legend.text=c("Min","Mean","Max"), names.arg=c("H3_1","H3_2","H4_1","H4_2","H2A_1","H2A_2","H2B_1","H2B_2","DNA_I","DNA_J"), col=c("green","blue","red"), axis.lty=1, xlab="Chains", ylab="RMSD, A")
#error.bar(barx,y.means, 1.96*y.sd/10)
#barx
#boxplot.matrix(bonds)

dev.off()