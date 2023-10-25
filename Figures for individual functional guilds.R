# creat plot
## plots for ACM
a=ggplot()+
  geom_segment(data=acm_effect[3:21,],color="tan",size=0.8,aes(y=1:19,yend=1:19,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=acm_effect[3:21,],pch=21,aes(x=Estimate,y=1:19,fill=sig),size=4)+
  
  scale_y_continuous(breaks=1:19,label=c("SoilOC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","MAP","Pre.seas.","SPEI","Pla.rich", "Fun.rich","MAT",expression(Root[mass]),expression(Root[d15N]),expression(Root[d13C]),expression(Root[C]),expression(Root[CN])))+
  xlab("")+
  ggtitle("ACM (N=87)")+
  scale_fill_manual("",breaks=c("sig","no"),values=c("forestgreen","gray"))+
  ylab("")+
  theme(legend.position = "", text = element_text(size=20), plot.title =element_text(hjust=0.5),axis.text.x = element_text(hjust = 0),axis.text.y = element_text(hjust = 0),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))

b=ggplot()+
  geom_segment(data=ecm_effect[3:20,],color="tan",size=0.8,aes(y=1:18,yend=1:18,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=ecm_effect[3:20,],pch=21,aes(x=Estimate,y=1:18,fill=sig),size=4)+
  
  scale_y_continuous(breaks=1:18,label=c("SoilOC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","Pre.seas.","SPEI", "Pla.rich", "Fun.rich","MAT",expression(Root[mass]),expression(Root[d15N]),expression(Root[d13C]),expression(Root[C]),expression(Root[CN])))+
  xlab("")+
  ggtitle("ECM (N=104)")+
  scale_fill_manual("",breaks=c("sig","no"),values=c("forestgreen","gray"))+
  ylab("")+
  theme(legend.position = "", plot.title =element_text(hjust=0.5) , text = element_text(size=20), axis.text.x = element_text(hjust = 0),axis.text.y = element_text(hjust = 0),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))


c=ggplot()+
  geom_segment(data=soilsap_effect[3:20,],color="tan",size=0.8,aes(y=1:18,yend=1:18,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=soilsap_effect[3:20,],pch=21,aes(x=Estimate,y=1:18,fill=sig),size=4)+
  
  scale_y_continuous(breaks=1:18,label=c("SoilOC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","Pre.seas.","SPEI", "Pla.rich", "Fun.rich","MAT",expression(Root[mass]),expression(Root[d15N]),expression(Root[d13C]),expression(Root[C]),expression(Root[CN])))+
  xlab("Effect size and 95% CI")+
  ggtitle("Soilsapro. (N=104)")+
  scale_fill_manual("",breaks=c("sig","no"),values=c("forestgreen","gray"))+
  ylab("")+
  theme(legend.position = "", text = element_text(size=20), plot.title =element_text(hjust=0.5),axis.text.x = element_text(hjust = 0),axis.text.y = element_text(hjust = 0),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))

d=ggplot()+
  geom_segment(data=litsap_effect[3:20,],color="tan",size=0.8,aes(y=1:18,yend=1:18,x=Estimate  -1.96*Std..Error,xend=Estimate+1.96*Std..Error))+
  geom_vline(xintercept=0,color="red",linetype="dashed",size=1.2)+
  geom_point(data=litsap_effect[3:20,],pch=21,aes(x=Estimate,y=1:18,fill=sig),size=4)+
  
  scale_y_continuous(breaks=1:18,label=c("SoilOC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","Pre.seas.","SPEI", "Pla.rich", "Fun.rich","MAT",expression(Root[mass]),expression(Root[d15N]),expression(Root[d13C]),expression(Root[C]),expression(Root[CN])))+
  xlab("Effect size and 95% CI")+
  ggtitle("Lit.sapro. (N=104)")+
  scale_fill_manual("",breaks=c("sig","no"),values=c("forestgreen","gray"))+
  ylab("")+
  theme(legend.position = "", text = element_text(size=20), plot.title =element_text(hjust=0.5),axis.text.x = element_text(hjust = 0),axis.text.y = element_text(hjust = 0),axis.title.y = element_text(face= "italic",size=20),axis.title.x = element_text(size=20),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))
