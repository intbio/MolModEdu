#R-script analysis of ion intein interactions
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
#INT-intein interactions
# ion_int<-read.csv('../analysis_data/ion_int_raw_df.csv')
ion_int_cryst<-read.csv('../analysis_data/ion_int_avr_df_cryst.csv')
ion_int<-read.csv('../analysis_data/ion_int_avr_df.csv')

dna_prot<-subset(dna_prot,type!='Z')
dna_prot_cryst<-subset(dna_prot_cryst,type!='Z')

#Let's define levels order for better graphics
int_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS','ADE','THY','GUA','CYT')
type_lev=c('IP','HB','VdW','IM','WM')

#Let's assign level order
ion_int_cryst$type<-factor(ion_int_cryst$type,levels=type_lev)
ion_int$type<-factor(ion_int$type,levels=type_lev)

ion_int$INT_part<-factor(ion_int$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))
ion_int_cryst$INT_part<-factor(ion_int_cryst$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))


ion_int_cryst$INT_resname<-factor(ion_int_cryst$INT_resname,levels=int_rn_lev)
ion_int$INT_resname<-factor(ion_int$INT_resname,levels=int_rn_lev)



###Symmetrization
sym<-function(chain,ri){
if(chain=='CHI') {return(ri)}
else {return(ri*(-1))}
}

ion_int_cryst=transform(ion_int_cryst,INT_resid=mapply(sym,INT_chain,INT_resid))
ion_int=transform(ion_int,INT_resid=mapply(sym,INT_chain,INT_resid))
#protein is messed up!!!


####General statistics section:
theme_set(theme_gray(base_size = 15))



##############################
#Now let's make profiles along DNA sequence
dnaseq<-"ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
seqdf=data.frame(sequence=substring(dnaseq, seq(1,nchar(dnaseq)),seq(1,nchar(dnaseq))),X=seq(-73,73,1))


################Total number of interactions by strand
q=seq(73,-73,-1)

###Crystal
d1=ddply(subset(ion_int_cryst,(INT_chain %in% c('CHI','CHJ'))&type %in% c('VdW','HB','IP')),c("INT_chain","INT_resid"),function(df) c(num=nrow(df)))


d1z=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
d2z=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHJ'),num=c(0))


d1<-rbind(d1,d1z,d2z)

tot_cryst=ddply(d1,c("INT_resid",'INT_chain'),summarize,number=sum(num))

c<-ggplot(data=tot_cryst) + ggtitle("Number of iteractions between DNA and ions, chain I and J(-1), total (VdW+IP+HB) by strand: crystal") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+xlim(0,80)
c<-c+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=5)

# ggsave(filename="../analysis_data/int_ion_int_dna_prof_tot_cryst.png",plot=c,width=15,height=5)




##MD----------
d1=ddply(subset(ion_int,(INT_chain %in% c('CHI','CHJ'))&type %in% c('VdW','HB','IP')),c("INT_chain","INT_resid"),function(df) c(num=sum(df$av_num)))


d1z=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
d2z=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHJ'),num=c(0))


d1<-rbind(d1,d1z,d2z)

tot_cryst=ddply(d1,c("INT_resid",'INT_chain'),summarize,number=sum(num))

m<-ggplot(data=tot_cryst) + ggtitle("Number of iteractions between DNA and ions, chain I and J(-1), total (VdW+IP+HB) by strand: MD") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_chain),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+xlim(0,80)
m<-m+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=5)

q<-arrangeGrob(c,m)

ggsave(filename="../analysis_data/int_ion_int_dna_prof_tot_sym.png",plot=q,width=15,height=10)


################Profile of interactions with bases, sugars, phosphates from MD

###Crystal
d1=ddply(subset(ion_int_cryst,(INT_chain %in% c('CHI','CHJ'))&type %in% c('VdW','HB','IP')),c("INT_chain","INT_resid",'INT_part'),function(df) c(num=nrow(df)))


d1zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
d2zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHJ'),num=c(0))

t=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHI'))
t2=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHJ'))

d1z=merge(d1zt,t)
d2z=merge(d2zt,t2)

d1<-rbind(d1,d1z,d2z)

tot_cryst=ddply(d1,c("INT_resid",'INT_chain','INT_part'),summarize,number=sum(num))

i<-ggplot(data=subset(tot_cryst,INT_chain=='CHI')) + ggtitle("Number of iteractions between DNA and ions, total (VdW+IP+HB), chain I: crystal") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+xlim(0,80)
i<-i+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=5)


j<-ggplot(data=subset(tot_cryst,INT_chain=='CHJ')) + ggtitle("Number of iteractions between DNA and ions, total (VdW+IP+HB), chain J(-1): crystal") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+xlim(0,80)
j<-j+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=5)

q<-arrangeGrob(i,j)


ggsave(filename="../analysis_data/int_ion_int_dna_prof_part_cryst_sym.png",plot=q,width=15,height=10)


#MD-----

d1=ddply(subset(ion_int,(INT_chain %in% c('CHI','CHJ'))&type %in% c('VdW','HB','IP')),c("INT_chain","INT_resid",'INT_part'),function(df) c(num=sum(df$av_num)))


d1zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
d2zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHJ'),num=c(0))

t=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHI'))
t2=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHJ'))

d1z=merge(d1zt,t)
d2z=merge(d2zt,t2)

d1<-rbind(d1,d1z,d2z)

tot_md=ddply(d1,c("INT_resid",'INT_chain','INT_part'),summarize,number=sum(num))

i<-ggplot(data=subset(tot_md,INT_chain=='CHI')) + ggtitle("Number of iteractions between DNA and ions, total (VdW+IP+HB), chain I: MD") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+xlim(0,80)
i<-i+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=5)



j<-ggplot(data=subset(tot_md,INT_chain=='CHJ')) + ggtitle("Number of iteractions between DNA and ions, total (VdW+IP+HB), chain J(-1): MD") + 
xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+xlim(0,80)
j<-j+geom_text(data=seqdf,aes(x=X,y=-0.2,label=sequence),size=5)


q<-arrangeGrob(i,j)

ggsave(filename="../analysis_data/int_ion_int_dna_prof_part_sym.png",plot=q,width=15,height=10)


