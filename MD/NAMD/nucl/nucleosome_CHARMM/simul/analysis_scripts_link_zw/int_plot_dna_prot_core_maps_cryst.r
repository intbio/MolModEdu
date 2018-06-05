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

dna_prot_cryst=subset(dna_prot_cryst, (((PROT_chain %in% c('CHA','CHE'))&(PROT_resid > 43)&(PROT_resid < 132))|((PROT_chain %in% c('CHB','CHF'))&(PROT_resid > 23)&(PROT_resid < 99))|((PROT_chain %in% c('CHC','CHG'))&(PROT_resid > 15)&(PROT_resid < 118))|((PROT_chain %in% c('CHD','CHH'))&(PROT_resid > 29)&(PROT_resid < 124))))
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')
type_lev=c('SC','SB','HB','vdW','IM','WM')
dna_prot_cryst$type<-factor(dna_prot_cryst$type,levels=type_lev)
dna_prot_cryst$DNA_part<-factor(dna_prot_cryst$DNA_part,levels=c('phosphate','sugar','base'))
dna_prot_cryst$PROT_resname<-factor(dna_prot_cryst$PROT_resname,levels=prot_rn_lev)

# prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
dna_prot_cryst$PROT_part<-factor(dna_prot_cryst$PROT_part,levels=c('backbone','side chain'))



img <- readPNG(paste("seq_img/",'H3r',".png",sep=''))
h3r <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4r',".png",sep=''))
h4r <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Ar',".png",sep=''))
h2ar <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Br',".png",sep=''))
h2br <- rasterGrob(img, interpolate=TRUE,width=1)


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

c_b=subset(c,(PROT_chain=='CHB')&(PROT_resid >=24) & (PROT_resid <=98))
c_f=subset(c,(PROT_chain=='CHF')&(PROT_resid >=24) & (PROT_resid <=98))

c_c=subset(c,(PROT_chain=='CHC')&(PROT_resid >=16) & (PROT_resid <=117))
c_g=subset(c,(PROT_chain=='CHG')&(PROT_resid >=16) & (PROT_resid <=117))

c_d=subset(c,(PROT_chain=='CHD')&(PROT_resid >=30) & (PROT_resid <=120))
c_h=subset(c,(PROT_chain=='CHH')&(PROT_resid >=30) & (PROT_resid <=120))



if(length(c_a$num)){
a<-ggplot(data=c_a,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H3 chain A')+xlab('DNA')+
ggtitle("Interactions of H3 chain A with DNA")+
geom_hline(yintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h3r, xmin=-49.4,xmax=-33,ymin=43.5, ymax=131.5)+

scale_x_continuous(limits=c(-49.0,73),breaks=seq(-70,70,10),minor_breaks = seq(-70, 70, 5))+
scale_y_continuous(limits=c(43.5,131.5),minor_breaks = seq(44, 131, 2))

ggsave(filename="../analysis_data/int_dna_prot_core_map_A_cryst.png",plot=a,width=7,height=7)

}

if(length(c_b$num)){
b<-ggplot(data=c_b,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H4 chain B')+xlab('DNA')+
ggtitle("Interactions of H3 chain A with DNA")+
geom_hline(yintercept = c(24,29,30,41,49,75,82,93), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h4r, xmin=-37.52, xmax=-30, ymin=23.5,ymax=98.5)+

scale_y_continuous(limits=c(20,98.5),minor_breaks = seq(24, 98, 2))+
scale_x_continuous(limits=c(-40.0,10),breaks=seq(-70,70,10),minor_breaks = seq(-70, 70, 2))

ggsave(filename="../analysis_data/int_dna_prot_core_map_B_cryst.png",plot=b,width=7,height=7)

}


if(length(c_c$num)){
c<-ggplot(data=c_c,aes(x=DNA_resid,y=PROT_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+

ylab('H2A chain C')+xlab('DNA')+
ggtitle("Interactions of H2A chain C with DNA")+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2ar, xmin=-70.6,xmax=-65,ymin=15.5, ymax=117.5)+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 2))+
scale_x_continuous(limits=c(-75.0,-20),breaks=seq(-70,70,10),minor_breaks = seq(-70, 70, 2))


ggsave(filename="../analysis_data/int_dna_prot_core_map_C_cryst.png",plot=c,width=7,height=7)

}


if(length(c_d$num)){
d<-ggplot(data=c_d,aes(x=DNA_resid,y=PROT_resid+3,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2B chain D')+xlab('DNA')+
ggtitle("Interactions of H2B chain D with DNA")+
geom_hline(yintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
# geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2br, xmin=-66.36, xmax=-60, ymin=32.5,ymax=123.5)+
# annotation_custom(h2ar, xmin=22.07,xmax=32.5,ymin=15.5, ymax=117.5)+

scale_y_continuous(limits=c(30,123.5),minor_breaks = seq(33, 123, 2))+
scale_x_continuous(limits=c(-70.0,-20.0),breaks=seq(-70,70,10),minor_breaks = seq(-70, 70, 2))

ggsave(filename="../analysis_data/int_dna_prot_core_map_D_cryst.png",plot=d,width=7,height=7)

}

quit()


if(length(c_e$num)){
e<-ggplot(data=c_e,aes(x=PROT2_resid,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2A chain C')+xlab('H3 chain E')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h3, ymin=117.5, ymax=127, xmin=43.5,xmax=131.5)+
annotation_custom(h2ar, xmin=33.5,xmax=43.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(33,131.5),minor_breaks = seq(43, 131, 5))+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 5))

ggsave(filename="../analysis_data/int_prot_prot_full_cryst_H2AZmap_E.png",plot=e,width=7,height=7)

}



if(length(c_f$num)){
f<-ggplot(data=c_f,aes(x=PROT2_resid,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2A chain C')+xlab('H4 chain F')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(24,29,30,41,49,75,82,93), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h4, ymin=117.5, ymax=127, xmin=23.5,xmax=98.5)+
annotation_custom(h2ar, xmin=15.0,xmax=23.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(13.5,98.5),minor_breaks = seq(24, 98, 5))+
scale_y_continuous(limits=c(14,127),minor_breaks = seq(16, 127, 5))

ggsave(filename="../analysis_data/int_prot_prot_full_cryst_H2AZmap_F.png",plot=f,width=7,height=7)

}




if(length(c_g$num)){
g<-ggplot(data=c_g,aes(x=PROT2_resid,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+

ylim(15.5,127)+xlim(2.5,117.5)+ylab('H2A chain C')+xlab('H2A chain G')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2a, ymin=117.5, ymax=127, xmin=15.5,xmax=117.5)+
annotation_custom(h2ar, xmin=3.83,xmax=15.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(2.5,117.5),minor_breaks = seq(0, 117.5, 5))+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 5))


ggsave(filename="../analysis_data/int_prot_prot_full_cryst_H2AZmap_G.png",plot=g,width=7,height=7)

}



if(length(c_h$num)){
h<-ggplot(data=c_h,aes(x=PROT2_resid+3,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylim(15.5,127)+xlim(20.5,123.5)+ylab('H2A chain C')+xlab('H2B chain H')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2b, ymin=117, ymax=127, xmin=32.5,xmax=123.5)+
annotation_custom(h2ar, xmin=22.07,xmax=32.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(20.5,123.5),minor_breaks = seq(33, 123, 5))+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 5))

ggsave(filename="../analysis_data/int_prot_prot_full_cryst_H2AZmap_H.png",plot=h,width=7,height=7)

}

# print(c_cb)
# a<-ggplot(data=c_ca,aes(x=PROT2_resid,y=PROT1_resid))+geom_point(aes(color=num))+
# scale_colour_gradient(low="blue", high="red")+
# ylim(16,117)+ylab('H2A residue number')

# a<-a+ggtitle('Interactions between histones in nucleosome, num = SB+HB+vdW: crystal')
# annotation_custom(h2ar, ymin=15.5, ymax=117.5, xmin=-10,xmax=10)+
# annotation_custom(h3, ymin=1.75, ymax=2.0, xmin=43.5,xmax=131.5)
# q<-arrangeGrob(a,b,e,f,c,d,g,h,ncol=2)
	q<-arrangeGrob(e,f,c,d,g,h,ncol=3)
ggsave(filename="../analysis_data/int_prot_prot_full_cryst_H2AZmap.png",plot=q,width=21,height=14)

quit()

#MD-----

nn=data.frame(CHAIN_name=c('CHA','CHB','CHC','CHD','CHE','CHF','CHG','CHH'),HIST_name=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

t=ddply(prot_prot[prot_prot$type %in% c('vdW','HB','SB'),],c('PROT1_chain','PROT2_chain','PROT1_resid','PROT2_resid'),summarize,num=sum(av_num))
t2=merge(t,nn,by.x=c('PROT1_chain'),by.y=c('CHAIN_name'))
t3=merge(t2,nn,by.x=c('PROT2_chain'),by.y=c('CHAIN_name'),suffixes=c('_1','_2'))
t4=subset(t3,select = !(colnames(t3) %in% c('PROT1_chain','PROT2_chain')))
# t5=melt(t4,id=c('HIST_name_1','HIST_name_2'))
# t6=transform(t5,X=1)
# c_md=rbind(t4,data.frame(HIST_name_1='H2B H',HIST_name_2='H3 A',PROT1_resid=60,PROT2_resid=60,num=NA))
c_md=t4

c_md$HIST_name_1<-factor(c_md$HIST_name_1,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))
c_md$HIST_name_2<-factor(c_md$HIST_name_2,levels=c('H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'))

c_mdr=c_md
names(c_mdr)<-c("PROT1_resid" ,"PROT2_resid", "num","HIST_name_2","HIST_name_1")
c_md=rbind(c_md,c_mdr)

a<-ggplot(data=c_md,aes(x=PROT1_resid,y=PROT2_resid))+geom_point(aes(color=num))+
facet_grid(HIST_name_1~HIST_name_2)+scale_colour_gradient(low="blue", high="red")

a<-a+ggtitle('Interactions between histones in nucleosome, num = SB+HB+vdW: MD')
ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap.png",plot=a,width=7,height=7)




