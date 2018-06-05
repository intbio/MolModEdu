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
dna_prot_cryst<-read.csv('../analysis_data/dna_prot_avr_df_cryst.csv')
dna_prot_cryst<-subset(dna_prot_cryst,type=='SC')

# dna_prot_cryst=subset(dna_prot_cryst, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SC','SB','HB','vdW','IM','WM')
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)

# prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
dna_prot_cryst$PROT_part<-factor(dna_prot_cryst$PROT_part,levels=c('backbone','side chain'))



img <- readPNG(paste("seq_img/",'H3fullr',".png",sep=''))
h3r <- rasterGrob(img, interpolate=TRUE,height=1)
h3r2 <- rasterGrob(img, interpolate=TRUE,height=1)


img <- readPNG(paste("seq_img/",'H4fullr',".png",sep=''))
h4r <- rasterGrob(img, interpolate=TRUE,height=1)


img <- readPNG(paste("seq_img/",'H2Afullr',".png",sep=''))
h2ar <- rasterGrob(img, interpolate=TRUE,height=1)
h2ar2 <- rasterGrob(img, interpolate=TRUE,height=1)
h2ar3 <- rasterGrob(img, interpolate=TRUE,height=1)


img <- readPNG(paste("seq_img/",'H2Bfullr',".png",sep=''))
h2br <- rasterGrob(img, interpolate=TRUE,height=1)
h2br2 <- rasterGrob(img, interpolate=TRUE,height=1)


####General statistics section:
theme_set(theme_gray(base_size = 15))



# Idea 2. Map along all sequences.

#Should correspond to cross-correlatoin maps
#Cryst ----

t=ddply(dna_prot_cryst,c('DNA_chain','PROT_chain','DNA_resid','PROT_resid'),summarize,num=length(param1))
#symmetrization
#combine DNA
q=seq(93,-93,-1)
d1i=t[t$DNA_chain=='CHI',]
d1j=t[t$DNA_chain=='CHJ',]
d1j$DNA_resid<-q[d1j$DNA_resid+94]
c<-rbind(d1i,d1j)
c=subset(c,num>0)

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


if(length(c_a$num)){
a<-ggplot(data=c_a,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H3 chain A')+xlab('DNA')+
ggtitle("Interactions of H3 chain A with DNA, X-ray")+
geom_hline(yintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h3r, xmin=-30-d,xmax=-30,ymin=0.5, ymax=135.5)+
annotation_custom(h3r2, xmin=60-d,xmax=60,ymin=0.5, ymax=135.5)+
# annotation_custom(h3r3, xmin=-10-d3,xmax=-10,ymin=0.5, ymax=135.5)+
# annotation_custom(h3r4, xmin=0-d4,xmax=0,ymin=0.5, ymax=135.5)+
# annotation_custom(h3r5, xmin=10-d5,xmax=10,ymin=0.5, ymax=135.5)+


scale_x_continuous(limits=c(-40,80),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))+
scale_y_continuous(limits=c(0.5,135.5),breaks=seq(0,140,10),minor_breaks = seq(0, 136, 2))

ggsave(filename="../analysis_data/int_dna_prot_full_map_A_cryst.png",plot=a,width=14,height=7)

}

if(length(c_b$num)){
b<-ggplot(data=c_b,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H4 chain B')+xlab('DNA')+
ggtitle("Interactions of H4 chain B with DNA, X-ray")+
geom_hline(yintercept = c(24,29,30,41,49,76,82,93), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h4r, xmin=-40.53, xmax=-34, ymin=0.5,ymax=102.5)+
# annotation_custom(h4r2, xmin=-42.0, xmax=-34, ymin=0.5,ymax=102.5)+

scale_y_continuous(limits=c(0.5,102.5),breaks=seq(0,110,10),minor_breaks = seq(0, 102, 2))+
scale_x_continuous(limits=c(-40.0,0),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))

ggsave(filename="../analysis_data/int_dna_prot_full_map_B_cryst.png",plot=b,width=7,height=7)

}


d=7.6


if(length(c_c$num)){
c<-ggplot(data=c_c,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain C')+xlab('DNA')+
ggtitle("Interactions of H2A chain C with DNA, X-ray")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=-80-d,xmax=-80,ymin=0.5, ymax=128.5)+
annotation_custom(h2ar2, xmin=-10-d,xmax=-10,ymin=0.5, ymax=128.5)+
annotation_custom(h2ar3, xmin=30-d,xmax=30,ymin=0.5, ymax=128.5)+


scale_y_continuous(limits=c(0.5,128.5),minor_breaks = seq(16, 127, 2))+
scale_x_continuous(limits=c(-90.0,40),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))


ggsave(filename="../analysis_data/int_dna_prot_full_map_C_cryst.png",plot=c,width=14,height=7)

}



if(length(c_d$num)){
d<-ggplot(data=c_d,aes(x=DNA_resid,y=PROT_resid+3,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2B chain D')+xlab('DNA')+
ggtitle("Interactions of H2B chain D with DNA, X-ray")+
geom_hline(yintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2br, xmin=-60-1, xmax=-60, ymin=2.5,ymax=128.5)+
annotation_custom(h2ar2, xmin=20-1,xmax=20,ymin=2.5, ymax=128.5)+

scale_y_continuous(limits=c(2.5,128.5),breaks=seq(0,128,10),minor_breaks = seq(0, 128, 2))+
scale_x_continuous(limits=c(-70.0,40.0),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))

ggsave(filename="../analysis_data/int_dna_prot_full_map_D_cryst.png",plot=d,width=14,height=7)

}




if(length(c_e$num)){
e<-ggplot(data=c_e,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H3 chain E')+xlab('DNA')+
ggtitle("Interactions of H3 chain E with DNA, X-ray")+
geom_hline(yintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h3r, xmin=-95-1,xmax=-95,ymin=0.5, ymax=135.5)+
annotation_custom(h3r2, xmin=-10-1,xmax=-10,ymin=0.5, ymax=135.5)+


scale_x_continuous(limits=c(-100,30),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))+
scale_y_continuous(limits=c(0.5,135.5),breaks=seq(0,140,10),minor_breaks = seq(0, 136, 2))

ggsave(filename="../analysis_data/int_dna_prot_full_map_E_cryst.png",plot=e,width=14,height=7)

}

if(length(c_f$num)){
f<-ggplot(data=c_f,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H4 chain F')+xlab('DNA')+
ggtitle("Interactions of H4 chain F with DNA, X-ray")+
geom_hline(yintercept = c(24,29,30,41,49,76,82,93), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h4r, xmin=-10, xmax=0, ymin=0.5,ymax=102.5)+
# annotation_custom(h4r2, xmin=-42.0, xmax=-34, ymin=0.5,ymax=102.5)+

scale_y_continuous(limits=c(0.5,102.5),breaks=seq(0,110,10),minor_breaks = seq(0, 102, 2))+
scale_x_continuous(limits=c(-10.0,40),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))

ggsave(filename="../analysis_data/int_dna_prot_full_map_F_cryst.png",plot=f,width=7,height=7)

}





if(length(c_g$num)){
g<-ggplot(data=c_g,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain G')+xlab('DNA')+
ggtitle("Interactions of H2A chain G with DNA, X-ray")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=30-1,xmax=30,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar2, xmin=-10-1,xmax=-10,ymin=0.5, ymax=128.5)+
# annotation_custom(h2ar3, xmin=30-1,xmax=30,ymin=0.5, ymax=128.5)+


scale_y_continuous(limits=c(0.5,128.5),minor_breaks = seq(16, 127, 2))+
scale_x_continuous(limits=c(20.0,80),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))


ggsave(filename="../analysis_data/int_dna_prot_full_map_G_cryst.png",plot=g,width=14,height=7)

}



if(length(c_h$num)){
h<-ggplot(data=c_h,aes(x=DNA_resid,y=PROT_resid+3,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2B chain H')+xlab('DNA')+
ggtitle("Interactions of H2B chain H with DNA, X-ray")+
geom_hline(yintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2br, xmin=-40-1, xmax=-40, ymin=2.5,ymax=128.5)+
annotation_custom(h2ar2, xmin=30-1,xmax=30,ymin=2.5, ymax=128.5)+

scale_y_continuous(limits=c(2.5,128.5),breaks=seq(0,128,10),minor_breaks = seq(0, 128, 2))+
scale_x_continuous(limits=c(-50.0,60.0),breaks=seq(-90,90,10),minor_breaks = seq(-90, 90, 2))

ggsave(filename="../analysis_data/int_dna_prot_full_map_H_cryst.png",plot=h,width=14,height=7)

}
