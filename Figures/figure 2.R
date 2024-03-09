# codes that create figure 2. these models don't include root traits.

#________________ when the estimated c was included in the model____________
library(lme4)
library(lmerTest)
library(sjPlot)
library(dplyr)
library(MuMIn)
library(ggcorrplot)
library(cowplot)
library(ggplot2)

1. # for arbuscular mycorrhizal fungi
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio1 + bio2 + bio4 + bio8 + bio12 + bio15 + bio18 + spei + richness + funrich + (1 | siteIDD / plotID), data = acm_model_rich)
step(mod)
mod_acm_rich <- lmer(z ~ c + organicCPercent + bio15 + funrich + (1 | siteIDD / plotID), data = acm_model_rich)

p1 <- plot_model(mod_acm_rich, axis.labels = c("Fun.rich", "Pre.seas.", "SoilC"), colors = c("royalblue1", "red"), rm.terms = "c", title = "ACM (N=319)", axis.lim = c(-1, 1))

2. # for ectomycorrhizal fungi
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = ecm_model_rich)
step(mod)
mod_ecm_rich <- lmer(z ~ c + ph + funrich + (1 | siteIDD / plotID), data = ecm_model_rich)

p2 <- plot_model(mod_ecm_rich, axis.labels = c("Fun.rich", "pH"), color = c("royalblue1", "red"), rm.terms = "c", title = "ECM (N=438)")

3. # for soil saprotroph

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = soilsap_model_rich)
step(mod)
mod_soilsap_rich <- lmer(z ~ c + nitrogen + sand + bio2 + bio18 + bio12 + spei + richness + funrich + (1 | siteIDD / plotID), data = soilsap_model_rich)

p3 <- plot_model(mod_soilsap_rich, axis.labels = c("Fun.rich", "Pla.rich", "Spei", "MAP", "Pre.WQ", "MDR", "Sand", "SoilN"), color = c("royalblue1", "red"), rm.terms = "c", title = "Soilsap. (N=438)")

4. # for plant pathogens

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = plapat_model_rich)
step(mod)
mod_plapat_rich <- lmer(z ~ c + organicCPercent + sand + bio15 + spei + richness + funrich + (1 | siteIDD / plotID), data = plapat_model_rich)
p4 <- plot_model(mod_plapat_rich, axis.labels = c("Fun.rich", "Pla.rich", "Spei", "Pre.seas", "Sand", "SoilC"), colors = c("royalblue1", "red"), rm.terms = "c", title = "Pla.patho. (N=427)")

5. # for litter saprotroph
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = litsap_model_rich)
step(mod)
mod_litsap_rich <- mod <- lmer(z ~ c + sand + bio12 + bio15 + spei + richness + funrich + (1 | siteIDD / plotID), data = litsap_model_rich)
p5 <- plot_model(mod_litsap_rich, axis.labels = c("Fun.rich", "Pla.rich", "Spei", "Pre.seas.", "MAP", "Sand"), colors = c("royalblue1", "red"), rm.terms = "c", title = "Lit.sap. (N=438)", axis.lim = c(-1, 1))

6. # for wood saprotroph
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = woosap_model_rich)
step(mod)
mod_woosap_rich <- lmer(z ~ c + sand + bio2 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), woosap_model_rich)
p6 <- plot_model(mod_woosap_rich, axis.labels = c("MAT", "Fun.rich", "Pla.rich", "Spei", "MAP", "MDR", "Sand"), colors = c("royalblue1", "red"), rm.terms = "c", title = "Woo.sap. (N=427)")

#________________ when the estimated c was replaced by the core-level richness____________-
# all variables were presented without model selection

range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

# for arbuscular mycorrhizal fungi
data=subset(rich_guild_com,aguild=="arbuscular_mycorrhizal")
data=merge(acm_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_acm=summary(mod)

effect_acm=effect_acm$coefficients
effect_acm=data.frame(effect_acm)[2:dim(effect_acm)[1],]

p1=ggplot()+
  geom_point(data=effect_acm,aes(x= Estimate,y=1:dim(effect_acm)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_acm,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_acm$Estimate-1.96*effect_acm$Std..Error,y=1:dim(effect_acm)[1],xend=effect_acm$Estimate+1.96*effect_acm$Std..Error,yend=1:dim(effect_acm)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.28,y=11,label="*",size=10)+
  ggtitle("ACM (N=319)")


# for ecm fungi

data=subset(rich_guild_com,aguild=="ectomycorrhizal")
data=merge(ecm_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_ecm=summary(mod)
effect_ecm=effect_ecm$coefficients
effect_ecm=data.frame(effect_ecm)[2:dim(effect_ecm)[1],]

p2=ggplot()+
  geom_point(data=effect_ecm,aes(x= Estimate,y=1:dim(effect_ecm)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_ecm,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_ecm$Estimate-1.96*effect_ecm$Std..Error,y=1:dim(effect_ecm)[1],xend=effect_ecm$Estimate+1.96*effect_ecm$Std..Error,yend=1:dim(effect_ecm)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.38,y=14,label="***",size=8)+
  annotate("text",x=0.28,y=11,label="**",size=8)+
  annotate("text",x=-0.28,y=1,label="***",size=8)+
  ggtitle("ECM (N=438)")

# for soil saprotrophic fungi

data=subset(rich_guild_com,aguild=="soil_saprotroph" )
data=merge(soilsap_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_soilsap=summary(mod)
effect_soilsap=effect_soilsap$coefficients
effect_soilsap=data.frame(effect_soilsap)[2:dim(effect_soilsap)[1],]

p3=ggplot()+
  geom_point(data=effect_soilsap,aes(x= Estimate,y=1:dim(effect_soilsap)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_soilsap,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_soilsap$Estimate-1.96*effect_soilsap$Std..Error,y=1:dim(effect_soilsap)[1],xend=effect_soilsap$Estimate+1.96*effect_soilsap$Std..Error,yend=1:dim(effect_soilsap)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.3,y=14,label="***",size=8)+
  annotate("text",x=0.208,y=11,label="*",size=8)+
  annotate("text",x=-0.32,y=1,label="***",size=8)+
  ggtitle("Soilsap. (N=438)")

# for plant pathogens

data=subset(rich_guild_com,aguild=="plant_pathogen"  )
data=merge(plapat_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
effect_plapat=summary(mod)
effect_plapat=effect_plapat$coefficients
effect_plapat=data.frame(effect_plapat)[2:dim(effect_plapat)[1],]

p4=ggplot()+
  geom_point(data=effect_plapat,aes(x= Estimate,y=1:dim(effect_plapat)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_plapat,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_plapat$Estimate-1.96*effect_plapat$Std..Error,y=1:dim(effect_plapat)[1],xend=effect_plapat$Estimate+1.96*effect_plapat$Std..Error,yend=1:dim(effect_plapat)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.3,y=14,label="",size=8)+
  annotate("text",x=0.28,y=11,label="*",size=8)+
  annotate("text",x=-0.5,y=1,label="***",size=8)+
  ggtitle("Pla.patho. (N=427)")

## for litter saprotroph

data=subset(rich_guild_com,aguild=="litter_saprotroph"   )
data=merge(litsap_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
effect_litsap=summary(mod)
effect_litsap=effect_litsap$coefficients
effect_litsap=data.frame(effect_litsap)[2:dim(effect_litsap)[1],]

p5=ggplot()+
  geom_point(data=effect_litsap,aes(x= Estimate,y=1:dim(effect_litsap)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_litsap,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_litsap$Estimate-1.96*effect_litsap$Std..Error,y=1:dim(effect_litsap)[1],xend=effect_litsap$Estimate+1.96*effect_litsap$Std..Error,yend=1:dim(effect_litsap)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.3,y=14,label="***",size=8)+
  annotate("text",x=0.28,y=11,label="*",size=8)+
  annotate("text",x=0.35,y=10,label="*",size=8)+
  annotate("text",x=-0.45,y=1,label="***",size=8)+
  ggtitle("Litt.sap. (N=438)")

#
data=subset(rich_guild_com,aguild=="wood_saprotroph"    )
data=merge(woosap_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
effect_woosap=summary(mod)
effect_woosap=effect_woosap$coefficients
effect_woosap=data.frame(effect_woosap)[2:dim(effect_woosap)[1],]

p6=ggplot()+
  geom_point(data=effect_woosap,aes(x= Estimate,y=1:dim(effect_woosap)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_woosap,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_woosap$Estimate-1.96*effect_woosap$Std..Error,y=1:dim(effect_woosap)[1],xend=effect_woosap$Estimate+1.96*effect_woosap$Std..Error,yend=1:dim(effect_woosap)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.3,y=14,label="",size=8)+
  annotate("text",x=0.28,y=11,label="",size=8)+
  annotate("text",x=0.35,y=10,label="",size=8)+
  annotate("text",x=-0.25,y=1,label="*",size=8)+
  ggtitle("Wood.sap. (N=427)")

plot_grid(p1,p2,p3,p4,p5,p6,ncol=2,labels=c("(a)","(b)","(c)","(d)","(e)","(f)"),label_x = 0.25)

# for epiphyte

data=subset(rich_guild_com,aguild=="epiphyte" )

data=merge(epiphy_model_rich,data,by="plotID")
data=unique(data)
data[,6:29]=apply(data[,6:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
effect_epiphy=summary(mod)
effect_epiphy=effect_epiphy$coefficients
effect_epiphy=data.frame(effect_epiphy)[2:dim(effect_epiphy)[1],]

p7=ggplot()+
  geom_point(data=effect_epiphy,aes(x= Estimate,y=1:dim(effect_epiphy)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_epiphy,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_epiphy$Estimate-1.96*effect_epiphy$Std..Error,y=1:dim(effect_epiphy)[1],xend=effect_epiphy$Estimate+1.96*effect_epiphy$Std..Error,yend=1:dim(effect_epiphy)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.3,y=14,label="",size=8)+
  annotate("text",x=0.28,y=11,label="",size=8)+
  annotate("text",x=0.35,y=10,label="",size=8)+
  annotate("text",x=-0.25,y=2,label="*",size=8)+
  ggtitle("Epiphyte (N=392)")

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = para_model_rich)
effect_para=summary(mod)
effect_para=effect_para$coefficients
effect_para=data.frame(effect_para)[2:dim(effect_para)[1],]

p8=ggplot()+
  geom_point(data=effect_para,aes(x= Estimate,y=1:dim(effect_para)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_para,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_para$Estimate-1.96*effect_para$Std..Error,y=1:dim(effect_para)[1],xend=effect_para$Estimate+1.96*effect_para$Std..Error,yend=1:dim(effect_para)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
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
  annotate("text",x=0.35,y=14,label="**",size=8)+
  annotate("text",x=-0.5,y=13,label="**",size=8)+
  annotate("text",x=0.35,y=11,label="***",size=8)+
  annotate("text",x=0.35,y=10,label="",size=8)+
  annotate("text",x=-0.30,y=1,label="*",size=8)+
  ggtitle("Parasitic (N=409)")


##
data=subset(rich_guild_com,aguild=="soil_saprotroph")
data=merge(soilsap_model_rich,data,by="plotID")

data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

mod <- lmer(z ~ Observed + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = data)

step(mod)
mod_soilsap_rich <- lmer(z ~ c + nitrogen + sand + bio2 + bio18 + bio12 + spei + richness + funrich + (1 | siteIDD / plotID), data = soilsap_model_rich)

p3 <- plot_model(mod_soilsap_rich, axis.labels = c("Fun.rich", "Pla.rich", "Spei", "MAP", "Pre.WQ", "MDR", "Sand", "SoilN"), color = c("royalblue1", "red"), rm.terms = "c", title = "Soilsap. (N=438)")

