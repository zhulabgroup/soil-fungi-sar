#effect_all_clima_soil

library(ggplot2)

a=ggplot()+
  geom_segment(data=effect_all_clima_soil[3:15,],color="tan",size=0.8,aes(y=1:13,yend=1:13,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=effect_all_clima_soil[3:15,],aes(x=Estimate,y=1:13),size=4)+
  scale_y_continuous(breaks=1:13,label=c("SoilOC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","MAP","Pre.seas.","SPEI","Fun.rich","MAT"))+
  xlab("Effect size and 95% CI")+
  ggtitle("NEON+DoB (N=483)")+
  ylab("")+
  theme(legend.position = "bottom", text = element_text(size=20), axis.text.x = element_text(hjust = 1),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))
# for the NEON site
#effect_neon_clima_soil
b=ggplot()+
  geom_segment(data=effect_neon_clima_soil[c(2,4:16),],color="tan",size=0.8,aes(y=1:14,yend=1:14,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=effect_neon_clima_soil[c(2,4:16),],aes(x=Estimate,y=1:14),size=4)+

  scale_y_continuous(breaks=1:14,label=c("SoilOC","pH","SoilN","Pla.rich.","Sand","Bio2","Bio8","Bio18","Tem.seas.","MAP","Pre.seas.","SPEI","Fun.rich","MAT"))+
  xlab("Effect size and 95% CI")+
  ggtitle("NEON (N=396)")+
  ylab("")+
  theme(legend.position = "bottom", text = element_text(size=20), axis.text.x = element_text(hjust = 1),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))

# climate, soil plant and root model

c=ggplot()+
  geom_segment(data=effect_neon_clima_soil_root[c(2,4:19),],color="tan",size=0.8,aes(y=1:17,yend=1:17,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=effect_neon_clima_soil_root[c(2,4:19),],aes(x=Estimate,y=1:17),size=4)+
  
  scale_y_continuous(breaks=1:17,label=c("SoilOC","pH","SoilN","Pla.rich.","Sand","Bio2","Bio8","Bio18","MAP","Pre.seas.","SPEI","Fun.rich","MAT",expression(Root[mass]),expression(Root[d13C]),expression(Root[C]),expression(Root[C:N])))+
  xlab("Effect size and 95% CI")+
  ggtitle("NEON (N=104)")+
  ylab("")+
  theme(legend.position = "bottom", text = element_text(size=20), axis.text.x = element_text(hjust = 1),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))

           