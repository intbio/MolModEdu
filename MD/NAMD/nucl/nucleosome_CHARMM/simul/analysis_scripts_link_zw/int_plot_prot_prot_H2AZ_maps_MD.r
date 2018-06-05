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
# prot_prot_cryst<-read.csv('../analysis_data/prot_prot_avr_df_cryst.csv')
prot_prot<-read.csv('../analysis_data/prot_prot_avr_df.csv')

#Let's define levels order for better graphics
prot_rn_lev=c('ARG','LYS','THR','SER','TYR','HSE','GLN','ASN','VAL','ILE','ALA','GLY','PRO','PHE','LEU','GLU','ASP','MET','CYS')

# prot_prot_cryst=subset(prot_prot_cryst,type!='SC')
prot_prot=subset(prot_prot,type!='SC')

#exclude intrachain interactions

# prot_prot_cryst=subset(prot_prot_cryst,PROT1_chain!=PROT2_chain)
# prot_prot=subset(prot_prot,PROT1_chain!=PROT2_chain)


type_lev=c('SB','HB','vdW','IM','WM')

#residues important for H2A.Z
	# vdw.changeSelection("segname CHC and noh and resid 30 33 35 37 38 39 40 41 48 62 71 73 75 76 79 89 94 95 97 98 100 104 110 117")
#Let's concentrate on them
h2azimp=c(30,33,35,37,38,39,40,41,48,62,71,73,75,76,79,89,94,95,97,98,100,104,110,117)

prot_prot=subset(prot_prot,((PROT1_chain=='CHC')&(PROT1_resid %in% h2azimp))|((PROT2_chain=='CHC')&(PROT2_resid %in% h2azimp)))



#Let's assign level order
# prot_prot_cryst$type<-factor(prot_prot_cryst$type,levels=type_lev)
prot_prot$type<-factor(prot_prot$type,levels=type_lev)

prot_prot$PROT1_part<-factor(prot_prot$PROT1_part,levels=c('backbone','side chain'))
# prot_prot_cryst$PROT1_part<-factor(prot_prot_cryst$PROT1_part,levels=c('backbone','side chain'))


# prot_prot_cryst$PROT1_resname<-factor(prot_prot_cryst$PROT1_resname,levels=prot_rn_lev)
prot_prot$PROT1_resname<-factor(prot_prot$PROT1_resname,levels=prot_rn_lev)

# prot_prot_cryst$PROT2_resname<-factor(prot_prot_cryst$PROT2_resname,levels=prot_rn_lev)
prot_prot$PROT2_resname<-factor(prot_prot$PROT2_resname,levels=prot_rn_lev)



img <- readPNG(paste("seq_img/",'H3',".png",sep=''))
h3 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H4',".png",sep=''))
h4 <- rasterGrob(img, interpolate=TRUE,width=1)


img <- readPNG(paste("seq_img/",'H2Ar',".png",sep=''))
h2ar <- rasterGrob(img, interpolate=TRUE,width=1)

img <- readPNG(paste("seq_img/",'H2A',".png",sep=''))
h2a <- rasterGrob(img, interpolate=TRUE,width=1)

img <- readPNG(paste("seq_img/",'H2B',".png",sep=''))
h2b <- rasterGrob(img, interpolate=TRUE,width=1)


####General statistics section:
theme_set(theme_gray(base_size = 15))



# Idea 2. Map along all sequences.

#Should correspond to cross-correlatoin maps
#Cryst ----

# t=ddply(prot_prot_cryst,c('PROT1_chain','PROT2_chain','PROT1_resid','PROT2_resid'),summarize,num=table(type)['SB']+table(type)['HB']+table(type)['vdW'])

t=ddply(prot_prot[prot_prot$type %in% c('vdW','HB','SB'),],c('PROT1_chain','PROT2_chain','PROT1_resid','PROT2_resid','PROT1_part','PROT2_part'),summarize,num=sum(av_num))

t=subset(t,num>0)
#symmetrization
t_r=t
names(t_r)<-c('PROT2_chain','PROT1_chain','PROT2_resid','PROT1_resid','PROT2_part','PROT1_part','num')
c=rbind(t,t_r)

#unsymmetry and only side chains for H2A
c=subset(c,(PROT1_chain=='CHC')&(PROT1_resid >=16) & (PROT1_resid <=117)&(PROT1_part=='side chain'))

# c=subset(c,(PROT1_chain=='CHC')&(PROT1_resid >=16) & (PROT1_resid <=117))

# a<-ggplot(data=c,aes(x=PROT2_resid,y=PROT1_resid))+geom_point(aes(color=num))+
# facet_grid(PROT2_chain~.)+scale_colour_gradient(low="blue", high="red")
c_a=subset(c,(PROT2_chain=='CHA')&(PROT2_resid >=44) & (PROT2_resid <=131))
c_e=subset(c,(PROT2_chain=='CHE')&(PROT2_resid >=44) & (PROT2_resid <=131))

c_b=subset(c,(PROT2_chain=='CHB')&(PROT2_resid >=24) & (PROT2_resid <=98))
c_f=subset(c,(PROT2_chain=='CHF')&(PROT2_resid >=24) & (PROT2_resid <=98))

c_c=subset(c,(PROT2_chain=='CHC')&(PROT2_resid >=16) & (PROT2_resid <=117))
c_g=subset(c,(PROT2_chain=='CHG')&(PROT2_resid >=16) & (PROT2_resid <=117))

c_d=subset(c,(PROT2_chain=='CHD')&(PROT2_resid >=30) & (PROT2_resid <=120))
c_h=subset(c,(PROT2_chain=='CHH')&(PROT2_resid >=30) & (PROT2_resid <=120))



if(length(c_a$num)){
a<-ggplot(data=c_a,aes(x=PROT2_resid,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2A chain C')+xlab('H3 chain A')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(44,57,63,77,85,114,120,131), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h3, ymin=117.5, ymax=127, xmin=43.5,xmax=131.5)+
annotation_custom(h2ar, xmin=33.5,xmax=43.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(33,131.5),minor_breaks = seq(43, 131, 5))+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 5))

ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_A.png",plot=a,width=7,height=7)

}

if(length(c_b$num)){
b<-ggplot(data=c_b,aes(x=PROT2_resid,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylab('H2A chain C')+xlab('H4 chain B')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(24,29,30,41,49,75,82,93), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h4, ymin=117.5, ymax=127, xmin=23.5,xmax=98.5)+
annotation_custom(h2ar, xmin=15.0,xmax=23.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(13.5,98.5),minor_breaks = seq(24, 98, 5))+
scale_y_continuous(limits=c(14,127),minor_breaks = seq(16, 127, 5))

ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_B.png",plot=b,width=7,height=7)

}


if(length(c_c$num)){

#unsymmetrize
c_c=subset(c_c,PROT1_resid %in% h2azimp)

c<-ggplot(data=c_c,aes(x=PROT2_resid,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+

ylim(15.5,127)+xlim(2.5,117.5)+ylab('H2A chain C')+xlab('H2A chain C')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2a, ymin=117.5, ymax=127, xmin=15.5,xmax=117.5)+
annotation_custom(h2ar, xmin=3.83,xmax=15.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(2.5,117.5),minor_breaks = seq(0, 117.5, 5))+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 5))


ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_C.png",plot=c,width=7,height=7)

}


if(length(c_d$num)){
d<-ggplot(data=c_d,aes(x=PROT2_resid+3,y=PROT1_resid,color=num))+
geom_point()+scale_colour_gradient(low="blue", high="red",name="Number")+
ylim(15.5,127)+xlim(20.5,123.5)+ylab('H2A chain C')+xlab('H2B chain D')+
ggtitle("Interactions of H2A/H2A.Z differential residues")+
geom_vline(xintercept = c(37,49,55,84,90,102,103,123), colour="green", linetype = "longdash",size=0.5)+
geom_hline(yintercept = c(16,22,26,37,46,73,79,89,90,97), colour="green", linetype = "longdash",size=0.5)+

annotation_custom(h2b, ymin=117, ymax=127, xmin=32.5,xmax=123.5)+
annotation_custom(h2ar, xmin=22.07,xmax=32.5,ymin=15.5, ymax=117.5)+

scale_x_continuous(limits=c(20.5,123.5),minor_breaks = seq(33, 123, 5))+
scale_y_continuous(limits=c(15.5,127),minor_breaks = seq(16, 127, 5))

ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_D.png",plot=d,width=7,height=7)

}



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

ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_E.png",plot=e,width=7,height=7)

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

ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_F.png",plot=f,width=7,height=7)

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


ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_G.png",plot=g,width=7,height=7)

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

ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap_H.png",plot=h,width=7,height=7)

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
ggsave(filename="../analysis_data/int_prot_prot_full_H2AZmap.png",plot=q,width=21,height=14)

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




