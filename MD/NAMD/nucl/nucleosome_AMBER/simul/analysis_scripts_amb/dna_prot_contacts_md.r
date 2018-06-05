#R-script to diplay DNA hbonds
library(matrixStats)
png(file="..//analysis_data//dna_prot_contacts_hist_md.png",width=11,height=5,units="in",res=1200)
df<-read.table("../analysis_data//dna_prot_contacts_md.dat",skip=4,header=TRUE)
bonds<-as.matrix(df[,-1])
min<-colMins(bonds)
max<-colMaxs(bonds)
mean<-colMeans(bonds)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

ratio<-mean/max
data<-rbind(ratio,mean-ratio,max-mean)
barx <- barplot(data, beside=FALSE, legend.text=c("Ratio mean/max","Mean","Max"), names.arg=-73:73, col=c("green","blue","red"), axis.lty=1, xlab="DNA bases", ylab="Number of contacts with protein")
#error.bar(barx,y.means, 1.96*y.sd/10)
#barx
#boxplot.matrix(bonds)

dev.off()