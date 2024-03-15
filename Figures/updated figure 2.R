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
# all variables were presented without model selection, and root traits were excluded

range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

# for arbuscular mycorrhizal fungi
data=subset(rich_guild_com,aguild=="arbuscular_mycorrhizal")
data=merge(acm_model_rich,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(5:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-soilC[for an earlier version the cec and bio4 were not included in the model]
mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
effect_acm=summary(mod)
effect_acm=effect_acm$coefficients
effect_acm=data.frame(effect_acm)[2:dim(effect_acm)[1],]
# best-fit model
step(mod)
mod_fit=lmer(z ~ organicCPercent + ph + nitrogen + bio15 + (1 | siteIDD/plotID),data=data)
effect_best_acm=summary(mod_fit)
effect_best_acm=effect_best_acm$coefficients
effect_best_acm=data.frame(effect_best_acm)[2:dim(effect_best_acm)[1],]

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
ggcorrplot(cor(data[, c(5:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-soilC
mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_ecm=summary(mod)
effect_ecm=effect_ecm$coefficients
effect_ecm=data.frame(effect_ecm)[2:dim(effect_ecm)[1],]
#best model
mod_best=lmer(z ~ Observed + bio4 + bio12 + bio15 + spei + richness + (1 | siteIDD/plotID),data=data)
effect_best_ecm=summary(mod_best)
effect_best_ecm=effect_best_ecm$coefficients
effect_best_ecm=data.frame(effect_best_ecm)[2:dim(effect_best_ecm)[1],]


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
ggcorrplot(cor(data[, c(5:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
#bold-soilC
mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_soilsap=summary(mod)
effect_soilsap=effect_soilsap$coefficients
effect_soilsap=data.frame(effect_soilsap)[2:dim(effect_soilsap)[1],]
step(mod)

mod_best=lmer(z ~ Observed + nitrogen + bio4 + bio12 + bio15 + richness + (1 | siteIDD/plotID),data=data)

effect_best_soilsap=summary(mod_best)#unconverge
effect_best_soilsap=effect_best_soilsap$coefficients
effect_best_soilsap=data.frame(effect_best_soilsap)[2:dim(effect_best_soilsap)[1],]

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
ggcorrplot(cor(data[, c(5:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-soilC
mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_plapat=summary(mod)
effect_plapat=effect_plapat$coefficients
effect_plapat=data.frame(effect_plapat)[2:dim(effect_plapat)[1],]
#best model

mod_best=lmer(z ~ Observed + (1 | siteIDD/plotID),data=data)
effect_best_plapat=summary(mod_best)
effect_best_plapat=effect_best_plapat$coefficients
effect_best_plapat=data.frame(effect_best_plapat)[2:dim(effect_best_plapat)[1],]

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
ggcorrplot(cor(data[, c(5:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_litsap=summary(mod)
effect_litsap=effect_litsap$coefficients
effect_litsap=data.frame(effect_litsap)[2:dim(effect_litsap)[1],]
# mod_best
mod_best=lmer(z ~ Observed + cec + nitrogen + sand + bio2 + bio4 + bio12 + richness + (1 | siteIDD/plotID),data=data)
effect_best_litsap=summary(mod_best)
effect_best_litsap=effect_best_litsap$coefficients
effect_best_litsap=data.frame(effect_best_litsap)[2:dim(effect_best_litsap)[1],]


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
ggcorrplot(cor(data[, c(5:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-soilC
mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)


effect_woosap=summary(mod)
effect_woosap=effect_woosap$coefficients
effect_woosap=data.frame(effect_woosap)[2:dim(effect_woosap)[1],]
#best model
mod_best=lmer(z ~ Observed + (1 | siteIDD/plotID),data=data)
effect_best_woosap=summary(mod_best)
effect_best_woosap=effect_best_woosap$coefficients
effect_best_woosap=data.frame(effect_best_woosap)[2:dim(effect_best_woosap)[1],]

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
ggcorrplot(cor(data[, c(6:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-soilC
mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)
mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_epiphy=summary(mod)
effect_epiphy=effect_epiphy$coefficients
effect_epiphy=data.frame(effect_epiphy)[2:dim(effect_epiphy)[1],]
# best model selection
mod_best=lmer(z ~ funrich + bio2 + (1 | siteIDD/plotID),data=data)
effect_best_epiphy=summary(mod_best)
effect_best_epiphy=effect_best_epiphy$coefficients
effect_best_epiphy=data.frame(effect_best_epiphy)[2:dim(effect_best_epiphy)[1],]

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
# for saprotroph fungi


data=merge(para_model,rich_papra_mean,by="plotID")
data=subset(data,siteIDD!="GUAN"&z<10&richness>0)
data=unique(data)
ggcorrplot(cor(data[, c(6:20,29)]), hc.order = TRUE, type = "lower", lab = TRUE)
data[,c(6:21,29)]=apply(data[,c(6:21,29)],2,range01)%>%data.frame()
#bold-soilC

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = para_model_rich)
mod <- lmer(z ~ Observed + funrich+organicCPercent +cec+ ph + nitrogen + sand +bio1+ bio2 +bio4+ bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = data)

effect_para=summary(mod)
effect_para=effect_para$coefficients
effect_para=data.frame(effect_para)[2:dim(effect_para)[1],]

mod_best=lmer(z ~ Observed + organicCPercent + bio4 + bio12 + bio15 + spei + richness + (1 | siteIDD/plotID),data=data)
effect_best_para=summary(mod_best)
effect_best_para=effect_best_para$coefficients
effect_best_para=data.frame(effect_best_para)[2:dim(effect_best_para)[1],]

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

# combine all the effect sizes
effect_acm=cbind(va=rownames(effect_acm),effect_acm)
effect_ecm=cbind(va=rownames(effect_ecm),effect_ecm)
effect_soilsap=cbind(va=rownames(effect_soilsap),effect_soilsap)
effect_plapat=cbind(va=rownames(effect_plapat),effect_plapat)
effect_litsap=cbind(va=rownames(effect_litsap),effect_litsap)
effect_woosap=cbind(va=rownames(effect_woosap),effect_woosap)
effect_epiphy=cbind(va=rownames(effect_epiphy),effect_epiphy)
effect_para=cbind(va=rownames(effect_para),effect_para)

effect_noroot=left_join(effect_acm[,c(1,2,6)],effect_ecm[,c(1,2,6)], by="va")
effect_noroot=left_join(effect_noroot,effect_soilsap[,c(1,2,6)],by="va")
effect_noroot=left_join(effect_noroot,effect_plapat[,c(1,2,6)],by="va")
effect_noroot=left_join(effect_noroot,effect_litsap[,c(1,2,6)],by="va")
effect_noroot=left_join(effect_noroot,effect_woosap[,c(1,2,6)],by="va")
effect_noroot=left_join(effect_noroot,effect_epiphy[,c(1,2,6)],by="va")
effect_noroot=left_join(effect_noroot,effect_para[,c(1,2,6)],by="va")
names(effect_noroot)=c("va","acm","pacm", "ecm","pecm","soilsap","psoilsap","plapat","pplapat","litsap", "plitsap","woosap", "pwoosap","epiphy","pepiphy","para","ppara")


effect_noroot_value=effect_noroot[,c(2,4,6,8,10,12,14,16)]
effect_noroot_pva=effect_noroot[,c(3,5,7,9,11,13,15,17)]
rownames(effect_noroot_value)=rownames(effect_acm)
rownames(effect_noroot_pva)=rownames(effect_acm)
effect_noroot_pva=as.matrix(effect_noroot_pva)
effect_noroot_value=as.matrix(effect_noroot_value)# need to be convered as a matrix

effect_noroot_value=melt(effect_noroot_value)
effect_noroot_pva=melt(effect_noroot_pva)
effect_noroot_pva$value[effect_noroot_pva$value<0.01]="**"
effect_noroot_pva$value[effect_noroot_pva$value>0.01&effect_noroot_pva$value<0.05]="*"
effect_noroot_pva$value[effect_noroot_pva$value>0.05]=""

effect_noroot_value=cbind(effect_noroot_value,effect_noroot_pva$value)

P1=ggplot(effect_noroot_value, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  geom_text(aes(x = X1, y = X2, label = effect_noroot_pva$value))+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
  scale_y_discrete(breaks=as.character(unique(effect_noroot_value$X2)),labels=c("ACM(N=319)","ECM(N=438)","Soil saprotroph(N=438)","Plant pathogen(N=427)","Litter saprotroph(N=438)","Wood saprotroph(N=427)","Epiphyte(N=392)","Parasitic(N=409)"))+
  scale_x_discrete(breaks=as.character(unique(effect_noroot_value$X1)),labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.","MDR","MAT","Sand","SoilN","pH","Cec","SoilC","Plot.rich","Core.rich")))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b=-0.5, unit="cm"))+
  ggtitle("Full model")+
  xlab("")+
  ylab("")
  


######
effect_best_acm=cbind(va=rownames(effect_best_acm),effect_best_acm)
effect_best_ecm=cbind(va=rownames(effect_best_ecm),effect_best_ecm)
effect_best_soilsap=cbind(va=rownames(effect_best_soilsap),effect_best_soilsap)
effect_best_plapat=cbind(va=rownames(effect_best_plapat),effect_best_plapat)
effect_best_litsap=cbind(va=rownames(effect_best_litsap),effect_best_litsap)
effect_best_woosap=cbind(va=rownames(effect_best_woosap),effect_best_woosap)
effect_best_epiphy=cbind(va=rownames(effect_best_epiphy),effect_best_epiphy)
effect_best_para=cbind(va=rownames(effect_best_para),effect_best_para)


va=cbind(va= c("Observed"  , "funrich" , "organicCPercent", "cec"   , "ph"   , "nitrogen"  , "sand", "bio1"  ,"bio2"  ,  "bio4" , "bio8" ,"bio12" , "bio15" ,"bio18" , "spei" , "richness"  ),code=c(1:16))%>%data.frame()

effect_best_noroot=left_join(va,effect_best_acm[,c(1,2,6)], by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_ecm[,c(1,2,6)], by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_soilsap[,c(1,2,6)],by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_plapat[,c(1,2,6)],by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_litsap[,c(1,2,6)],by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_woosap[,c(1,2,6)],by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_epiphy[,c(1,2,6)],by="va")
effect_best_noroot=left_join(effect_best_noroot,effect_best_para[,c(1,2,6)],by="va")
# delete the code
effect_best_noroot=effect_best_noroot[,-2]

names(effect_best_noroot)=c("va","acm","pacm", "ecm","pecm","soilsap","psoilsap","plapat","pplapat","litsap", "plitsap","woosap", "pwoosap","epiphy","pepiphy","para","ppara")
rownames(effect_best_noroot)=effect_best_noroot$va



effect_best_noroot_value=effect_best_noroot[,c(2,4,6,8,10,12,14,16)]
effect_best_noroot_pva=effect_best_noroot[,c(3,5,7,9,11,13,15,17)]
effect_best_noroot_value= melt(as.matrix(effect_best_noroot_value))
effect_best_noroot_pva= melt(as.matrix(effect_best_noroot_pva))
effect_best_noroot_value=cbind(effect_best_noroot_value,effect_best_noroot_pva$value)
names(effect_best_noroot_value)[4]="pv"
effect_best_noroot_value$pv[effect_best_noroot_value$pv<0.001]="***"
effect_best_noroot_value$pv[effect_best_noroot_value$pv<0.01&effect_best_noroot_value$pv>0.001]="**"
effect_best_noroot_value$pv[effect_best_noroot_value$pv<0.05&effect_best_noroot_value$pv>0.01]="*"
effect_best_noroot_value$pv[effect_best_noroot_value$pv>0.05]=""
effect_best_noroot_value$value[effect_best_noroot_value$value=="NA"]="0"
#

P2=ggplot(effect_best_noroot_value, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  geom_text(aes(x = X1, y = X2, label = pv))+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000",na.value = "gray50")+
  
  scale_y_discrete(breaks=as.character(unique(effect_best_noroot_value$X2)),labels=c("ACM(N=319)","ECM(N=438)","Soil saprotroph(N=438)","Plant pathogen(N=427)","Litter saprotroph(N=438)","Wood saprotroph(N=427)","Epiphyte(N=392)","Parasitic(N=409)"))+
  scale_x_discrete(breaks=as.character(unique(effect_best_noroot_value$X1)),labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.","MDR","MAT","Sand","SoilN","pH","Cec","SoilC","Plot.rich","Core.rich")))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.margin = margin(t=-0.1, unit="cm"),# reduce the space between individual plots
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15,angle=270,hjust=0,
        color=rev(c("seagreen1", "royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","peru","purple","purple"))))+
  xlab("")+
  ylab("")+
  ggtitle("Best-fit model")

# 
P1=ggplotGrob(P1)
P2=ggplotGrob(P2)


plot_grid(P1,P2,ncol=1,labels=c("(a)","(b)",label_x = 0.8),rel_heights = c(1,1.3),label_size = 14)

##

mod <- lmer(z ~ Observed + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = data)

step(mod)
mod_soilsap_rich <- lmer(z ~ c + nitrogen + sand + bio2 + bio18 + bio12 + spei + richness + funrich + (1 | siteIDD / plotID), data = soilsap_model_rich)

p3 <- plot_model(mod_soilsap_rich, axis.labels = c("Fun.rich", "Pla.rich", "Spei", "MAP", "Pre.WQ", "MDR", "Sand", "SoilN"), color = c("royalblue1", "red"), rm.terms = "c", title = "Soilsap. (N=438)")

