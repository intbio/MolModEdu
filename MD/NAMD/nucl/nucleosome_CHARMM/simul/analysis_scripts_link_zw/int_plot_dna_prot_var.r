#R-script analysis of dna protein interactions with variation during MD
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(gridExtra)
##############
###############
##Loading data frames
#DNA-protein interactions
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
# dna_prot_cryst<-read.csv('../analysis_data/dna_prot_raw_df_cryst.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')

################
#Here we need to filter out HB and IP from contacts, and recalcualte WM interactions.
############
#data frame SQL-like magic

a=split(dna_prot,dna_prot$type)
#get rid of duplicated in contacts via tricky merge and join
m=merge(a$C,a$IP,by.x=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','y'))
m2=merge(m,a$HB,by.x=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),by.y=c("Time","DNA_atom","DNA_chain","DNA_resid","PROT_atom","PROT_chain","PROT_resid"),all.x=TRUE,suffix=c('','z'))
a$C=m2[is.na(m2$typez) & is.na(m2$typey),-(12:20)]
#Let's simlify now the water mediated interactions
a$WM=a$WM[!duplicated(a$WM),]

dna_prot=do.call(rbind,a)

#-----------


#########
nf=length(table(dna_prot$Time))
# nf=750

theme_set(theme_gray(base_size = 18))

################Total number of interactions
q=seq(73,-73,-1)

##MD----------
d1=ddply(subset(dna_prot,type=='C' & DNA_chain=='CHI'),c('Time',"DNA_resid"),summarize, num=length(param1))

d1zt=data.frame(DNA_resid=seq(73,-73,-1),num=c(0))
t=data.frame(Time=seq(0,749,1))
d1z=merge(d1zt,t)
d1<-rbind(d1,d1z)
tot_md=ddply(d1,c('Time',"DNA_resid"),summarize,number=sum(num))

cont_var=ddply(tot_md,c('DNA_resid'),summarize, mean=mean(number),sd=sd(number))


cv<-ggplot(data=cont_var) + ggtitle("DNA-protein interactions in nucleosome, MD, chain I, number of contacts with stdev") + 
xlab("Base pair")+geom_line(aes(x=DNA_resid,y=mean),size=1)+geom_errorbar(aes(x=DNA_resid,ymin=mean-sd/2,ymax=mean+sd/2,color='blue'),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c("blue", "dark green"))#+#+xlim(-70,70)
# geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)

ggsave(filename="../analysis_data/int_dna_chi_cont_var.png",plot=cv,width=15,height=5)

##############################

# dna_prot_cryst$DNA_part
# dna_prot=dna_prot[dna_prot$type %in% c('HB','IM','WM','IP'),]
# theme_set(theme_grey(base_size=30))

# dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
# seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(0,nchar(dnaseq)-1))
# seqdf_t<-seqdf
# seqdf_t$X<-seqdf_t$X-73.0
# seqdf_t$Y<-rep(0,length(seqdf_t$X))
# #Let;s average interactions for every resid.
# q=seq(73,-73,-1)
# d1=ddply(dna_prot_cryst,c("DNA_chain","DNA_resid","type"),function(df) c(num=nrow(df)))
# d1i=d1[d1$DNA_chain=='CHI',]
# d1j=d1[d1$DNA_chain=='CHJ',]
# d1j$DNA_resid<-q[d1j$DNA_resid+74]
# d1<-rbind(d1i,d1j)
# d1=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num))

# # ggplot(data=d1,aes(x=DNA_resid,y=number,color=type))+geom_line()

# pl<-ggplot(data=d1) + ggtitle("DNA-protein interactions in nucleosome: crystal") + 
# xlab("Base pair")+geom_line(aes(x=DNA_resid,y=number,color=type),size=2)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))#+#+xlim(-70,70)
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)

# #collect graphs to one page and print them
# png(file="../analysis_data/int_dna_prot_cryst.png",width=2000,height=500)
# print(pl)

# nf=length(table(dna_prot$Time))
# d1=ddply(dna_prot,c("DNA_chain","DNA_resid","type"),function(df) c(num=nrow(df)))
# d1i=d1[d1$DNA_chain=='CHI',]
# d1j=d1[d1$DNA_chain=='CHJ',]
# d1j$DNA_resid<-q[d1j$DNA_resid+74]
# d1<-rbind(d1i,d1j)
# d1=ddply(d1,c("DNA_resid","type"),summarize,number=sum(num)/nf)

# pl<-ggplot(data=d1) + ggtitle("DNA-protein interactions in nucleosome: MD simulation") + 
# xlab("Base pair ")+geom_line(aes(x=DNA_resid,y=number,color=type),size=2)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))#+
# # geom_text(data=seqdf_t,aes_string(x='X',y='Y',label='sequence'),size=10)#+ylim(15,50)


# png(file="../analysis_data/int_dna_prot.png",width=2000,height=500)
# print(pl)
