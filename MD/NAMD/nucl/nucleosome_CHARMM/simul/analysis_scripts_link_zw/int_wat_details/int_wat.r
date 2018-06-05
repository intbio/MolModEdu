#R-script analysis of wat intein interactions
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
# wat_int<-read.csv('../analysis_data/wat_int_raw_df.csv')
wat_int_cryst<-read.csv('../../analysis_data/int_wat_details/wat_int_avr_df_3p9_cryst.csv')
wat_int<-read.csv('../../analysis_data/int_wat_details/wat_int_avr_df_3p9.csv')

#Let's define levels order for better graphics
int_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS','ADE','THY','GUA','CYT')

#Let's assign level order

wat_int$INT_part<-factor(wat_int$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))
wat_int_cryst$INT_part<-factor(wat_int_cryst$INT_part,levels=c('phosphate','sugar','base','backbone','side chain'))


wat_int_cryst$INT_resname<-factor(wat_int_cryst$INT_resname,levels=int_rn_lev)
wat_int$INT_resname<-factor(wat_int$INT_resname,levels=int_rn_lev)


####General statistics section:
theme_set(theme_gray(base_size = 15))

##########Histograms with INT_part classification
data_c=subset(wat_int_cryst,INT_resname %in% c('GUA','THY','CYT','ADE') & INT_part=='base')
data_c$data='Crystal'
data_m=subset(wat_int,INT_resname %in% c('GUA','THY','CYT','ADE') & INT_part=='base')
data_m$data='MD'
data=rbind(data_c,data_m)
data$data=as.factor(data$data)
#---Cryst
a1<-ggplot(data=data,aes(x=INT_resname,fill=data,weight=num_HB))+
geom_bar(aes(color=data,y=..count..),position='dodge')
c<-a1+ggtitle('Average number of contacts of bases with waters by nucleotide \n contact is defined as distance between heavy atoms less than 3.9 A')
# ggsave(filename="../../analysis_data/int_wat_details/int_wat_int_hist_part_cryst.png",plot=a1,width=13,height=5)

##------MD

# a1<-ggplot(data=subset(wat_int,INT_resname %in% c('GUA','THY','CYT','ADE')& INT_part=='base'),aes(x=INT_resname,fill=INT_part,weight=num_HB))+
# geom_bar(aes(color=INT_part,y=..count..),position='stack')
# m<-a1+ggtitle('Number of hydrogen bonds with water by residue or nucleotide: MD')


# q<-arrangeGrob(c,m)

ggsave(filename="../../analysis_data/int_wat_details/dna_wat_int_hist_3p9.png",plot=c,width=13,height=10)

##############################
#Now let's make profiles along DNA sequence


# ################Total number of interactions by strand
# q=seq(73,-73,-1)

# ###Crystal
# d1=ddply(subset(wat_int_cryst,INT_chain %in% c('CHI','CHJ')),c("INT_chain","INT_resid"),function(df) c(num=sum(df$num_HB)))

# # d1=ddply(subset(wat_int_cryst,INT_chain %in% c('CHI','CHJ')),c("INT_chain","INT_resid",'INT_part'),function(df) c(num=sum(df$num_HB)))
# # d1i=d1[d1$INT_chain=='CHI',]
# # d1j=d1[d1$INT_chain=='CHJ',]
# # d1j$INT_resid<-q[d1j$INT_resid+74]
# # #Let's add zeros also

# # d1zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
# # t=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHI'))
# # d1z=merge(d1zt,t)

# # d1<-rbind(d1i,d1j,d1z)

# # tot_cryst=ddply(d1,c("INT_resid",'INT_part'),summarize,number=sum(num))
# #Uncomment to get I+J*(-1) combined
# tot_cryst=d1

# c<-ggplot(data=tot_cryst) + ggtitle("Number of hydrogen bonds between DNA and water, total by strand: crystal") + 
# xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=num,color=INT_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))

# # ggsave(filename="../analysis_data/int_wat_int_dna_prof_tot_cryst.png",plot=c,width=15,height=5)


# ##MD----------

# d1=ddply(subset(wat_int,INT_chain %in% c('CHI','CHJ')),c("INT_chain","INT_resid"),function(df) c(num=sum(df$num_HB)))

# # d1=ddply(subset(wat_int_cryst,INT_chain %in% c('CHI','CHJ')),c("INT_chain","INT_resid",'INT_part'),function(df) c(num=sum(df$num_HB)))
# # d1i=d1[d1$INT_chain=='CHI',]
# # d1j=d1[d1$INT_chain=='CHJ',]
# # d1j$INT_resid<-q[d1j$INT_resid+74]
# # #Let's add zeros also

# # d1zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
# # t=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHI'))
# # d1z=merge(d1zt,t)

# # d1<-rbind(d1i,d1j,d1z)

# # tot_cryst=ddply(d1,c("INT_resid",'INT_part'),summarize,number=sum(num))
# #Uncomment to get I+J*(-1) combined
# tot=d1

# m<-ggplot(data=tot) + ggtitle("Number of hydrogen bonds between DNA and water, total by strand: MD") + 
# xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=num,color=INT_chain),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))
# q<-arrangeGrob(c,m)

# ggsave(filename="../../analysis_data/int_wat_details/int_wat_int_dna_prof_tot.png",plot=q,width=15,height=10)



# ################Profile of interactions with bases, sugars, phosphates
# q=seq(73,-73,-1)

# ###Crystal
# d1=ddply(subset(wat_int_cryst,INT_chain %in% c('CHI','CHJ')),c("INT_chain","INT_resid","INT_part"),function(df) c(num=sum(df$num_HB)))

# # d1j=d1[d1$INT_chain=='CHJ',]
# # d1j$INT_resid<-q[d1j$INT_resid+74]
# # #Let's add zeros also

# d1zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
# d2zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHJ'),num=c(0))

# t=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHI'))
# t2=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHJ'))

# d1z=merge(d1zt,t)
# d2z=merge(d2zt,t2)

# d1<-rbind(d1,d1z,d2z)

# tot_cryst=ddply(d1,c("INT_resid",'INT_part','INT_chain'),summarize,number=sum(num))
# #Uncomment to get I+J*(-1) combined
# # tot_cryst=d1

# i<-ggplot(data=subset(tot_cryst,INT_chain=='CHI')) + ggtitle("Number of hydrogen bonds between DNA and water, chain I: crystal") + 
# xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))
# j<-ggplot(data=subset(tot_cryst,INT_chain=='CHJ')) + ggtitle("Number of hydrogen bonds between DNA and water, chain J: crystal") + 
# xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))
# q<-arrangeGrob(i,j)

# ggsave(filename="../../analysis_data/int_wat_details/int_wat_int_dna_prof_part_cryst.png",plot=q,width=15,height=10)


# ##MD----------

# d1=ddply(subset(wat_int,INT_chain %in% c('CHI','CHJ')),c("INT_chain","INT_resid","INT_part"),function(df) c(num=sum(df$num_HB)))

# # d1j=d1[d1$INT_chain=='CHJ',]
# # d1j$INT_resid<-q[d1j$INT_resid+74]
# # #Let's add zeros also

# d1zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHI'),num=c(0))
# d2zt=data.frame(INT_resid=seq(73,-73,-1),INT_chain=c('CHJ'),num=c(0))

# t=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHI'))
# t2=data.frame(INT_part=c('phosphate','sugar','base'),INT_chain=c('CHJ'))

# d1z=merge(d1zt,t)
# d2z=merge(d2zt,t2)

# d1<-rbind(d1,d1z,d2z)

# tot_md=ddply(d1,c("INT_resid",'INT_part','INT_chain'),summarize,number=sum(num))
# #Uncomment to get I+J*(-1) combined
# # tot_cryst=d1

# i<-ggplot(data=subset(tot_md,INT_chain=='CHI')) + ggtitle("Number of hydrogen bonds between DNA and water, chain I: MD") + 
# xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))
# j<-ggplot(data=subset(tot_md,INT_chain=='CHJ')) + ggtitle("Number of hydrogen bonds between DNA and water, chain J: MD") + 
# xlab("Nucleotide")+geom_line(aes(x=INT_resid,y=number,color=INT_part),size=1)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))
# q<-arrangeGrob(i,j)

# ggsave(filename="../../analysis_data/int_wat_details/int_wat_int_dna_prof_part.png",plot=q,width=15,height=10)

