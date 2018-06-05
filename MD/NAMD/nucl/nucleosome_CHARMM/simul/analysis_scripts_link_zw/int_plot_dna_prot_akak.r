#R-script analysis of prot protein interactions
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(gridExtra)
library(raster)
##############
###############
##Loading data frames
#PROT1-protein interactions
# prot_prot<-read.csv('../analysis_data/prot_prot_raw_df.csv')
dna_prot<-read.csv('../analysis_data/dna_prot_avr_df.csv')
dna_prot<-subset(dna_prot,type=='SC')

# dna_prot=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SC','SB','HB','vdW','IM','WM')
dna_prot$type<-factor(dna_prot$type,levels=type_lev)
dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)

# prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
dna_prot$PROT_part<-factor(dna_prot$PROT_part,levels=c('backbone','side chain'))



img <- readPNG(paste("seq_img/",'H2Afullr',".png",sep=''))
h2ar <- rasterGrob(img, interpolate=TRUE,height=1)
h2ar2 <- rasterGrob(img, interpolate=TRUE,height=1)
h2ar3 <- rasterGrob(img, interpolate=TRUE,height=1)



####General statistics section:
theme_set(theme_gray(base_size = 10))



# Idea 2. Map along all sequences.

#Should correspond to cross-correlatoin maps
#Cryst ----

t=ddply(dna_prot,c('DNA_chain','PROT_chain','DNA_resid','PROT_resid','DNA_part'),summarize,num=sum(av_num))
#symmetrization
#combine DNA
q=seq(93,-93,-1)
d1i=t[t$DNA_chain=='CHI',]
d1j=t[t$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
c<-rbind(d1i,d1j)
c=subset(c,num>1)

# a<-ggplot(data=c,aes(x=PROT2_resid,y=PROT1_resid))+geom_point(aes(color=num))+
# facet_grid(PROT2_chain~.)+scale_colour_gradient(low="blue", high="red")
c_a=subset(c,(PROT_chain=='CHA'))
c_e=subset(c,(PROT_chain=='CHE'))

c_b=subset(c,(PROT_chain=='CHB'))
c_f=subset(c,(PROT_chain=='CHF'))

c_c=subset(c,(PROT_chain=='CHC'))
c_g=subset(c,(PROT_chain=='CHG'))

c_d=subset(c,(PROT_chain=='CHD'))
c_h=subset(c,(PROT_chain=='CHH'))


d=5.21




d=7.6


if(length(c_c$num)){
c<-ggplot(data=c_c,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point(size=3)+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain C')+xlab('DNA')+
ggtitle("Contacts of H2A chain C with DNA,MD")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=-60,xmax=-50,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar2, xmin=-10-d,xmax=-10,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar3, xmin=30-d,xmax=30,ymin=0.5, ymax=128.5)+


scale_y_continuous(limits=c(10.5,17.5),breaks=seq(-90,90,10),minor_breaks = seq(0, 128, 2))+
scale_x_continuous(limits=c(-60.0,-35),breaks=seq(-90,90,5),minor_breaks = seq(-90, 90, 2))+
facet_grid(DNA_part~.)


# ggsave(filename="../analysis_data/int_dna_prot_full_map_C.png",plot=c,width=14,height=7)

}




if(length(c_g$num)){
g<-ggplot(data=c_g,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point(size=3)+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain G')+xlab('DNA')+
ggtitle("Contacts of H2A chain G with DNA, MD")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=30-1,xmax=30,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar2, xmin=-10-1,xmax=-10,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar3, xmin=30-1,xmax=30,ymin=0.5, ymax=128.5)+


scale_y_continuous(limits=c(10.5,17.5),breaks=seq(0,128,5),minor_breaks = seq(0, 128, 2))+
scale_x_continuous(limits=c(20.0,50),breaks=seq(-90,90,5),minor_breaks = seq(-90, 90, 2))+
facet_grid(DNA_part~.)


# ggsave(filename="../analysis_data/int_dna_prot_full_map_G.png",plot=g,width=14,height=7)

}


q<-arrangeGrob(c,g,ncol=1)
ggsave(filename="../analysis_data/int_dna_prot_akak_MD.png",plot=q,width=4,height=8)

#Now for crystal

dna_prot<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')
dna_prot<-subset(dna_prot,type=='SC')

# dna_prot=subset(dna_prot, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SC','SB','HB','vdW','IM','WM')
dna_prot$type<-factor(dna_prot$type,levels=type_lev)
dna_prot$DNA_part<-factor(dna_prot$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot$PROT_resname<-factor(dna_prot$PROT_resname,levels=prot_rn_lev)

# prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
dna_prot$PROT_part<-factor(dna_prot$PROT_part,levels=c('backbone','side chain'))



img <- readPNG(paste("seq_img/",'H2Afullr',".png",sep=''))
h2ar <- rasterGrob(img, interpolate=TRUE,height=1)
h2ar2 <- rasterGrob(img, interpolate=TRUE,height=1)
h2ar3 <- rasterGrob(img, interpolate=TRUE,height=1)



####General statistics section:
theme_set(theme_gray(base_size = 10))



# Idea 2. Map along all sequences.

#Should correspond to cross-correlatoin maps
#Cryst ----

t=ddply(dna_prot,c('DNA_chain','PROT_chain','DNA_resid','PROT_resid','DNA_part'),summarize,num=length(param1))
#symmetrization
#combine DNA
q=seq(93,-93,-1)
d1i=t[t$DNA_chain=='CHI',]
d1j=t[t$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
c<-rbind(d1i,d1j)
c=subset(c,num>1)

# a<-ggplot(data=c,aes(x=PROT2_resid,y=PROT1_resid))+geom_point(aes(color=num))+
# facet_grid(PROT2_chain~.)+scale_colour_gradient(low="blue", high="red")
c_a=subset(c,(PROT_chain=='CHA'))
c_e=subset(c,(PROT_chain=='CHE'))

c_b=subset(c,(PROT_chain=='CHB'))
c_f=subset(c,(PROT_chain=='CHF'))

c_c=subset(c,(PROT_chain=='CHC'))
c_g=subset(c,(PROT_chain=='CHG'))

c_d=subset(c,(PROT_chain=='CHD'))
c_h=subset(c,(PROT_chain=='CHH'))


d=5.21




d=7.6


if(length(c_c$num)){
c<-ggplot(data=c_c,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point(size=3)+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain C')+xlab('DNA')+
ggtitle("Contacts of H2A chain C with DNA,X-ray")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=-65,xmax=-58,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar2, xmin=-10-d,xmax=-10,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar3, xmin=30-d,xmax=30,ymin=0.5, ymax=128.5)+


scale_y_continuous(limits=c(10.5,17.5),breaks=seq(-90,90,10),minor_breaks = seq(0, 128, 2))+
scale_x_continuous(limits=c(-60.0,-35),breaks=seq(-90,90,5),minor_breaks = seq(-90, 90, 2))+
facet_grid(DNA_part~.)


# ggsave(filename="../analysis_data/int_dna_prot_full_map_C.png",plot=c,width=14,height=7)

}




if(length(c_g$num)){
g<-ggplot(data=c_g,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point(size=3)+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain G')+xlab('DNA')+
ggtitle("Contacts of H2A chain G with DNA, X-ray")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=30-1,xmax=30,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar2, xmin=-10-1,xmax=-10,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar3, xmin=30-1,xmax=30,ymin=0.5, ymax=128.5)+


scale_y_continuous(limits=c(10.5,17.5),breaks=seq(0,128,5),minor_breaks = seq(0, 128, 2))+
scale_x_continuous(limits=c(20.0,50),breaks=seq(-90,90,5),minor_breaks = seq(-90, 90, 2))+
facet_grid(DNA_part~.)


# ggsave(filename="../analysis_data/int_dna_prot_full_map_G.png",plot=g,width=14,height=7)

}


q<-arrangeGrob(c,g,ncol=1)
ggsave(filename="../analysis_data/int_dna_prot_akak_cryst.png",plot=q,width=4,height=8)



