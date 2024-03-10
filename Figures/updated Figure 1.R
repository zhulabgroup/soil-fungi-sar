# load the required data
load("/model_data.RData")
load("/core_rich.RData")# core-level fungal richness for each plot

1.# climate and soil model

data_corich=merge(model_data[,1:20],core_rich,by="plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10)
data_corich[, 5:21] <- apply(data_corich[, 5:21], 2, range01) %>% data.frame()

# only bold and cec were excluded because of colinearity

mod <- lmer(z ~ corich+ funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+bio8+bio12+bio15+ bio18 + spei  + (1 | siteIDD / plotID), data = data_corich)
effect_CS=summary(mod)
effect_CS=effect_CS$coefficients
effect_CS=data.frame(effect_CS)[2:dim(effect_CS)[1],]

p01=ggplot()+
  geom_point(data=effect_CS,aes(x= Estimate,y=1:dim(effect_CS)[1]),
             color=rev(c("royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_CS,size=0.8,color=rev(c("royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_CS$Estimate-1.96*effect_CS$Std..Error,y=1:dim(effect_CS)[1],xend=effect_CS$Estimate+1.96*effect_CS$Std..Error,yend=1:dim(effect_CS)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c( "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ", "Tem.seas","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  annotate("text",x=0.28,y=11,label="",size=10)+
  annotate("text",x=0.182048,y=12,label="*",size=10)+
  annotate("text",x=-0.1592048,y=1,label="*",size=10)+
  annotate("text",x=0.425192048,y=9,label="**",size=10)+
  ggtitle("Clima.+Soil\n (N=483)")

2. # when climate, soil and plant richness were included

data_corich=merge(model_data[,1:20],core_rich,by="plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10&richness>0)
data_corich[, 5:21] <- apply(data_corich[, 5:21], 2, range01) %>% data.frame()
# onlyd bold and cec were included

mod <- lmer(z ~ corich+ funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+bio8+bio12+bio15+ bio18 + spei  +richness+ (1 | siteIDD / plotID), data = data_corich)
effect_CSP=summary(mod)
effect_CSP=effect_CSP$coefficients
effect_CSP=data.frame(effect_CSP)[2:dim(effect_CSP)[1],]

p02=ggplot()+
  geom_point(data=effect_CSP,aes(x= Estimate,y=1:dim(effect_CSP)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_CSP,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_CSP$Estimate-1.96*effect_CSP$Std..Error,y=1:dim(effect_CSP)[1],xend=effect_CSP$Estimate+1.96*effect_CSP$Std..Error,yend=1:dim(effect_CSP)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:15,labels = rev(c( "Pla.rich","Spei","Pre.WQ","Pre.seas.","MAP","MTWQ", "Tem.seas","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  annotate("text",x=0.28,y=11,label="",size=10)+
  annotate("text",x=0.182048,y=12,label="*",size=10)+
  annotate("text",x=-0.1592048,y=1,label="*",size=10)+
  annotate("text",x=0.425192048,y=9,label="**",size=10)+
  annotate("text",x=0.425192048,y=15,label="***",size=10)+
  ggtitle("Clima.+Soil+Pla.\n (N=438)")

3.# when climate, soil, plant richness and root traits were included

data_corich=merge(model_data[,1:27],core_rich,by="plotID")
data_corich=subset(data_corich,siteIDD != "GUAN" & z < 10 & richness > 0 & fine > 0 & rootc > 0)
data_corich[, 5:27] <- apply(data_corich[, 5:27], 2, range01) %>% data.frame()
ggcorrplot(cor(data_corich[, c(5:28)]), hc.order = TRUE, type = "lower", lab = TRUE)

mod <- lmer(z ~ corich+ funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+bio8+bio12+bio15+ bio18 + spei  +richness+ fine+d15N+d13C+rootc+rootcn+(1 | siteIDD / plotID), data = data_corich)

effect_CSPR=summary(mod)
effect_CSPR=effect_CSPR$coefficients
effect_CSPR=data.frame(effect_CSPR)[2:dim(effect_CSPR)[1],]

p03=ggplot()+
  geom_point(data=effect_CSPR,aes(x= Estimate,y=1:dim(effect_CSPR)[1]),
             color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_CSPR,size=0.8,color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_CSPR$Estimate-1.96*effect_CSPR$Std..Error,y=1:dim(effect_CSPR)[1],xend=effect_CSPR$Estimate+1.96*effect_CSPR$Std..Error,yend=1:dim(effect_CSPR)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:20,labels = rev(c(expression("Root"["cn"]),expression("Root"["c"]),"d13C" ,"d15N",expression("Root"["fmass"]),"Pla.rich","Spei","Pre.WQ","Pre.seas.","MAP","MTWQ", "Tem.seas","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  
  annotate("text",x=0.2592048,y=2,label="*",size=10)+
  annotate("text",x=0.325192048,y=15,label="***",size=10)+
  ggtitle("Clima.+Soil+Pla.+\n Root (N=104)")
# create a df to combine all the three figures

effect_ALL=merge(effect_CSPR[,c("va","Estimate","Std..Error")],effect_CS[,c("va","Estimate","Std..Error")],by="va",all.x=TRUE,sort=FALSE)
effect_ALL=merge(effect_ALL,effect_CSP[,c("va","Estimate","Std..Error")],by="va",all.x=TRUE,sort=FALSE)
names(effect_ALL)=c("va","es","sd","es1","sd1","es2","sd2")
effect_ALL[is.na(effect_ALL)]=0

# for the climate and soil model

p001=ggplot()+
  geom_point(data=effect_ALL,aes(x= es1,y=1:dim(effect_ALL)[1]),
             color=rev(c("gray","gray","gray","gray","gray","gray","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_ALL,size=0.8,color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_ALL$es1-1.96*effect_ALL$sd1,y=1:dim(effect_ALL)[1],xend=effect_ALL$es1+1.96*effect_ALL$sd1,yend=1:dim(effect_ALL)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:20,labels = rev(c(expression("Root"["cn"]),expression("Root"["c"]),"d13C" ,"d15N",expression("Root"["fmass"]),"Pla.rich","Spei","Pre.WQ","Pre.seas.","MAP","MTWQ", "Tem.seas","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  annotate("text",x=0.3592048,y=9,label="*",size=10)+
  annotate("text",x=-0.15325192048,y=1,label="*",size=10)+
  annotate("text",x=0.182053025192048,y=12,label="**",size=10)+
  ggtitle("Clima.+Soil \n (N=483)")+
  xlab("")

# for the climate,soil and plant model
p002=ggplot()+
  geom_point(data=effect_ALL,aes(x= es2,y=1:dim(effect_ALL)[1]),
             color=rev(c("gray","gray","gray","gray","gray","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_ALL,size=0.8,color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_ALL$es2-1.96*effect_ALL$sd2,y=1:dim(effect_ALL)[1],xend=effect_ALL$es2+1.96*effect_ALL$sd2,yend=1:dim(effect_ALL)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:20,labels = rev(c(expression("Root"["cn"]),expression("Root"["c"]),"d13C" ,"d15N",expression("Root"["fmass"]),"Pla.rich","Spei","Pre.WQ","Pre.seas.","MAP","MTWQ", "Tem.seas","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  annotate("text",x=0.2592048,y=2,label="",size=10)+
  annotate("text",x=0.3025192048,y=15,label="***",size=10)+
  annotate("text",x=-0.185325192048,y=1,label="**",size=10)+
  annotate("text",x=0.3029851325192048,y=9,label="*",size=10)+
  annotate("text",x=0.15325192048,y=12,label="*",size=10)+
  annotate("text",x=0.275325192048,y=11,label="**",size=10)+
  ggtitle("Clima.+Soil+Pla.\n (N=438)")+
  xlab("Effect size and 95% CI")


# for the climate, soil and plant model

p003=ggplot()+
  geom_point(data=effect_ALL,aes(x= es,y=1:dim(effect_ALL)[1]),
             color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_ALL,size=0.8,color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_ALL$es-1.96*effect_ALL$sd,y=1:dim(effect_ALL)[1],xend=effect_ALL$es+1.96*effect_ALL$sd,yend=1:dim(effect_ALL)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:20,labels = rev(c(expression("Root"["cn"]),expression("Root"["c"]),"d13C" ,"d15N",expression("Root"["fmass"]),"Pla.rich","Spei","Pre.WQ","Pre.seas.","MAP","MTWQ", "Tem.seas","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  
  annotate("text",x=0.2592048,y=2,label="*",size=10)+
  annotate("text",x=0.325192048,y=15,label="***",size=10)+
  ggtitle("Clima.+Soil+Pla.+\n Root (N=104)")+
  xlab("")


