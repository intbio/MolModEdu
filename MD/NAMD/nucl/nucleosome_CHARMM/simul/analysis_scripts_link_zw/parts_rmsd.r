#R-script to diplay DNA hbonds
library(matrixStats)
png(file="..//analysis_data//parts_rmsd.png",width=20,height=8,units="in",res=300)
df<-read.table("../analysis_data//parts_rmsd.dat",skip=4,header=TRUE)
rmsd<-as.matrix(df[500:1000,-1])
min<-colMins(rmsd)
max<-colMaxs(rmsd)
mean<-colMeans(rmsd)
print(min)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
par(omi=c(1,1,1,1))
par(las=2)
data<-rbind(min,mean-min,max-mean)
#barx <- barplot(data, beside=FALSE, legend.text=c("Min","Mean","Max"), names.arg=c("Ac_ca","H3c_ca","H4c_ca","H2Ac_ca","H2Bc_ca","H3_1","H3_2","H4_1","H4_2","H2A_1","H2A_2","H2B_1","H2B_2","DNA_I","DNA_J"), col=c("green","blue","red"), axis.lty=1, xlab="Chains", ylab="RMSD, A")
barx <- barplot(data, beside=FALSE, legend.text=c("Min","Mean","Max"), col=c("green","blue","red"), axis.lty=1, xlab="", ylab="RMSD, A")



dev.off()