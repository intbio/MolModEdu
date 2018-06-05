#R-script to analyze average dna propterties
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(ggplot2)
library(reshape2)
library(xtable)
library(plyr)
library(gridExtra)
###################
##Loading data frames
#DNA-protein interactions
# dna_prot<-read.csv('../analysis_data/dna_prot_raw_df.csv')
saxs_avr<-read.csv('../analysis_data/saxs_df_md_avr.csv')

# d_c=melt(dna_cryst,id.vars=c('BPnum'),measure.vars=c('x','y','z','Roll','Twist','Slide','chi_1','chi_2'))
d=melt(saxs_avr,id.vars=c('Vector'),measure.vars=c('Int_in_sol_av','Int_in_vac_av'))
# d_c$data='X-ray'
# d_md$data='MD average'
# d=rbind(d_c,d_md)

# print(d)
theme_set(theme_bw(base_size = 18))
###Roll and z-profiles
r<-ggplot(data=d,aes(x=Vector,y=value,color=variable)) + ggtitle("Theoretical SAXS intensity in solution") + 
xlab("q, 1/Å")+geom_line(size=1)+xlim(0,0.1)
# +geom_point(aes(color=data),size=3)+
# scale_x_continuous(breaks = round(seq(-80,80, by = 10),1))+scale_color_manual(values=c('blue','red','green'))+
# theme(panel.grid.major = element_line(colour = "grey90",size=2),panel.grid.minor = element_blank())+
# geom_text(data=seqdf,aes(x=X,y=min(d[d$variable=='Roll'&d$BPnum<147,]$value,na.rm=TRUE)*1.1,label=sequence),size=3)
ggsave(filename="../analysis_data/saxs.png",plot=r,width=15,height=10)


quit()



