library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(sjPlot)
library(ggcorrplot)


range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

# when root traits were included for eight individual guilds
load("rich_guild_com.RData")# core-level data for different guilds
load("rich_papra_mean.RData")# core-level data for parasitic fungi only (included several guidls and was determined speratedly)
load("acm_model.RData")
load("ecm_model.RData")
load("soilsap_model.RData")
load("plapat_model.RData")
load("woosap_model.RData")
load("litsap_model.RData")
load("epiphy_model.RData")
load("para_model.RData")

1.# for arbuscular mycorrhizal fungi
data=subset(rich_guild_com,aguild=="arbuscular_mycorrhizal")
data=merge(acm_model,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()

ggcorrplot(cor(data[, c(5:29)]), hc.order = TRUE, type = "lower", lab = TRUE)

#bold-soilC show colinearity patterns
#rootn-rootcn
mod <- lmer(z ~ Observed+funrich+organicCPercent +cec+ ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei  + coarse+fine + d13C + d15N+rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_acm=summary(mod)
effect_acm=effect_acm$coefficients
effect_acm=data.frame(effect_acm)[2:dim(effect_acm)[1],]

mod_best=lmer(z ~ funrich + organicCPercent + cec + sand + bio1 + bio2 + bio8 + bio12 + d13C + (1 | plotID:siteIDD),data=data)

effect_bestroot_acm=summary(mod_best)
effect_bestroot_acm=effect_bestroot_acm$coefficients
effect_bestroot_acm=data.frame(effect_bestroot_acm)[2:dim(effect_bestroot_acm)[1],]

ggplot()+
  geom_point(data=effect_acm,aes(x= Estimate,y=1:dim(effect_acm)[1]),
             color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1", "royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_acm,size=0.8,color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1", "seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","peru","purple","purple")),aes(x=effect_acm$Estimate-1.96*effect_acm$Std..Error,y=1:dim(effect_acm)[1],xend=effect_acm$Estimate+1.96*effect_acm$Std..Error,yend=1:dim(effect_acm)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:22,labels = rev(c("Pla.rich", expression("Root"["cn"]), expression("Root"["c"]),  "d15N","d13C", expression("Root"["fmass"]), expression("Root"["cmass"]),  "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.","MDR","MAT","Sand","SoilN","pH","Cec","SoilC","Plot.rich","Core.rich")))+
  theme( panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-1,1)+
  ylab("")+
  annotate("text",x=0.28,y=11,label="*",size=10)+
  ggtitle("ACM (N=319)")



2.# for ectomycorrhizal fungi

data=subset(rich_guild_com,aguild=="ectomycorrhizal")
data=merge(ecm_model,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(5:29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-cec
#coarse-fine root mass
#rootcn-rootn
mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C + d15N+rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_ecm=summary(mod)
effect_ecm=effect_ecm$coefficients
effect_ecm=data.frame(effect_ecm)[2:dim(effect_ecm)[1],]

step(mod)

 mod_best=lmer(z ~ Observed + cec + nitrogen + sand + bio2 + richness + (1 | siteIDD/plotID),data=data)

 effect_bestroot_ecm=summary(mod_best)
 effect_bestroot_ecm=effect_bestroot_ecm$coefficients
 effect_bestroot_ecm=data.frame(effect_bestroot_ecm)[2:dim(effect_bestroot_ecm)[1],]
 
3. # for soil saprotrophs
data=subset(rich_guild_com,aguild=="soil_saprotroph" )
data=merge(soilsap_model,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(5:29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#bold-soilC
#coarse-fine
#rootcn-rootn

mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C +d15N+ rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_soilsap=summary(mod)
effect_soilsap=effect_soilsap$coefficients
effect_soilsap=data.frame(effect_soilsap)[2:dim(effect_soilsap)[1],]

mod_best=lmer(z ~ funrich + cec + nitrogen + sand + bio4 + bio8 + bio12 + spei + d15N + rootc + richness + (1 | siteIDD/plotID),data=data)
effect_bestroot_soilsap=summary(mod_best)
effect_bestroot_soilsap=effect_bestroot_soilsap$coefficients
effect_bestroot_soilsap=data.frame(effect_bestroot_soilsap)[2:dim(effect_bestroot_soilsap)[1],]


4.# for plant pathogens

data=subset(rich_guild_com,aguild=="plant_pathogen")
data=merge(plapat_model,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(5:29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#coarse-fine
#bold-soilC
#rootn-rootcn

mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C +d15N+ rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)
# the low-p value variables were selected
effect_plapat=summary(mod)
effect_plapat=effect_plapat$coefficients
effect_plapat=data.frame(effect_plapat)[2:dim(effect_plapat)[1],]

mod_best=lmer(z ~ Observed + cec + nitrogen + sand + bio1 + bio4 + bio12 + rootc + richness + (1 | plotID:siteIDD),data=data)
effect_bestroot_plapat=summary(mod_best)
effect_bestroot_plapat=effect_bestroot_plapat$coefficients
effect_bestroot_plapat=data.frame(effect_bestroot_plapat)[2:dim(effect_bestroot_plapat)[1],]


5. # for litter_saprotroph
data=subset(rich_guild_com,aguild=="litter_saprotroph")
data=merge(litsap_model,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(5:29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#coarse-fine
#bold-soilc
#rootn-rootcn

mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C+d15N + rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_litsap=summary(mod)
effect_litsap=effect_litsap$coefficients
effect_litsap=data.frame(effect_litsap)[2:dim(effect_litsap)[1],]

mod_best=lmer(z ~ Observed + cec + nitrogen + sand + bio1 + bio4 + bio12 + fine + rootc + richness + (1 | siteIDD/plotID),data=data)

effect_bestroot_litsap=summary(mod_best)
effect_bestroot_litsap=effect_bestroot_litsap$coefficients
effect_bestroot_litsap=data.frame(effect_bestroot_litsap)[2:dim(effect_bestroot_litsap)[1],]


# the low-p value variables were selected
6.# for wood saptrophic
data=subset(rich_guild_com,aguild=="wood_saprotroph"    )
data=merge(woosap_model,data,by="plotID")
data=unique(data)
data[,28:29]=apply(data[,28:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(5:29)]), hc.order = TRUE, type = "lower", lab = TRUE)

#coarse-fine
#bold-soilc
#rootn-rootcn

mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C +d15N+ rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_woosap=summary(mod)
effect_woosap=effect_woosap$coefficients
effect_woosap=data.frame(effect_woosap)[2:dim(effect_woosap)[1],]

mod_best=lmer(z ~ cec + nitrogen + sand + bio4 + bio8 + bio12 + bio15 + rootc + richness + (1 | plotID:siteIDD),data=data)
effect_bestroot_woosap=summary(mod_best)
effect_bestroot_woosap=effect_bestroot_woosap$coefficients
effect_bestroot_woosap=data.frame(effect_bestroot_woosap)[2:dim(effect_bestroot_woosap)[1],]

7.# for epiphytes
data=subset(rich_guild_com,aguild=="epiphyte" )
data=merge(epiphy_model,data,by="plotID")
data=unique(data)
data=subset(data,siteIDD!="GUAN"& z < 10 & richness > 0 & fine > 0 & rootc > 0)
data[,6:29]=apply(data[,6:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(6:29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#coarse-fine
#bold-soilc
#rootn-rootcn

mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C + d15N+rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_epiphy=summary(mod)
effect_epiphy=effect_epiphy$coefficients
effect_epiphy=data.frame(effect_epiphy)[2:dim(effect_epiphy)[1],]
mod_best=lmer(z ~ cec + nitrogen + sand + bio1 + bio4 + bio8 + bio12 + bio15 + rootc + richness + (1 | plotID:siteIDD),data=data)

effect_bestroot_epiphy=summary(mod_best)
effect_bestroot_epiphy=effect_bestroot_epiphy$coefficients
effect_bestroot_epiphy=data.frame(effect_bestroot_epiphy)[2:dim(effect_bestroot_epiphy)[1],]

8.# for parasitic fungi

data=merge(para_model,rich_papra_mean,by="plotID")
data=unique(data)
data=subset(data,siteIDD!="GUAN"& z < 10 & richness > 0 & fine > 0 & rootc > 0)

data[,6:29]=apply(data[,6:29],2,range01)%>%data.frame()
ggcorrplot(cor(data[, c(6:29)]), hc.order = TRUE, type = "lower", lab = TRUE)
#coarse-fine
#bold-soilC
#rootn-rootcn

mod <- lmer(z ~ Observed+funrich+organicCPercent +cec + ph + nitrogen + sand + bio1+bio2 +bio4+ bio8 +bio12 + bio15 +bio18 + spei +fine + d13C +d15N+ rootc + rootcn +richness + (1 | siteIDD / plotID), data = data)

effect_para=summary(mod)
effect_para=effect_para$coefficients
effect_para=data.frame(effect_para)[2:dim(effect_para)[1],]
step(mod)

mod_best=lmer(z ~ Observed + nitrogen + sand + d15N + richness + (1 | plotID:siteIDD),data=data)

effect_bestroot_para=summary(mod_best)
effect_bestroot_para=effect_bestroot_para$coefficients
effect_bestroot_para=data.frame(effect_bestroot_para)[2:dim(effect_bestroot_para)[1],]

# combined the effects of the eight guilds
effect_acm=cbind(va=rownames(effect_acm),effect_acm)
effect_ecm=cbind(va=rownames(effect_ecm),effect_ecm)
effect_soilsap=cbind(va=rownames(effect_soilsap),effect_soilsap)
effect_plapat=cbind(va=rownames(effect_plapat),effect_plapat)
effect_litsap=cbind(va=rownames(effect_litsap),effect_litsap)
effect_woosap=cbind(va=rownames(effect_woosap),effect_woosap)
effect_epiphy=cbind(va=rownames(effect_epiphy),effect_epiphy)
effect_para=cbind(va=rownames(effect_para),effect_para)

dm1=left_join(effect_acm[,c(1,2,6)],effect_ecm[,c(1,2,6)], by="va")
dm1=left_join(dm1,effect_soilsap[,c(1,2,6)],by="va")
dm1=left_join(dm1,effect_plapat[,c(1,2,6)],by="va")
dm1=left_join(dm1,effect_litsap[,c(1,2,6)],by="va")
dm1=left_join(dm1,effect_woosap[,c(1,2,6)],by="va")
dm1=left_join(dm1,effect_epiphy[,c(1,2,6)],by="va")
dm1=left_join(dm1,effect_para[,c(1,2,6)],by="va")

names(dm1)=c("va","acm","pacm", "ecm","pecm","soilsap","psoilsap","plapat","pplapat","litsap", "plitsap","woosap", "pwoosap","epiphy","pepiphy","para","ppara")
dm2=dm1[,c(2,4,6,8,10,12,14,16)]
dmp=dm1[,c(3,5,7,9,11,13,15,17)]

rownames(dm2)=dm1$va
rownames(dmp)=dm1$va
dm2[is.na(dm2)]=0
dmp[is.na(dmp)]="-"
va=rep(rownames(dm2),times=8)

dm3=cbind(va,melt(dm2))
dmpp=cbind(va,melt(dmp))# the data must be a

# oder the display of different variables
dm3$va=factor(dm3$va,levels=unique(dm3$va))

dmpp$value[dmpp$value<0.01]="**"
dmpp$value[dmpp$value>=0.01&dmpp$value<0.05]="*"
dmpp$value[dmpp$value>0.05]=""

dm3=cbind(dm3,dmpp$value)

ggplot(dm3, aes(x = va, y = variable, fill = value)) +
geom_tile(color = "white", lwd = 1,linetype = 1)+
  geom_text(aes(x = va, y = variable, label = dm3$sig))+
scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
scale_y_discrete(breaks=as.character(unique(dm3$variable)),labels=c("ACM(N=87)","ECM(N=104)","Soil saprotroph(N=104)","Plant pathogen(N=103)","Litter saprotroph(N=104)","Wood saprotroph(N=103)","Epiphyte(N=103)","Parasitic(N=98)"))+
scale_x_discrete(breaks=as.character(unique(dm3$va)),labels = rev(c("Pla.rich", expression("Root"["cn"]), expression("Root"["c"]),  "d15N","d13C", expression("Root"["fmass"]), expression("Root"["cmass"]),  "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.","MDR","MAT","Sand","SoilN","pH","Cec","SoilC","Plot.rich","Core.rich")))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.margin = margin(b=-0.8, unit="cm"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.ticks =element_blank())+
        xlab("")+
       ylab("")+
  ggtitle("Full model")

  ## after modeling selection

effect_bestroot_acm=cbind(va=rownames(effect_bestroot_acm),effect_bestroot_acm)
effect_bestroot_ecm=cbind(va=rownames(effect_bestroot_ecm),effect_bestroot_ecm)
effect_bestroot_soilsap=cbind(va=rownames(effect_bestroot_soilsap),effect_bestroot_soilsap)
effect_bestroot_plapat=cbind(va=rownames(effect_bestroot_plapat),effect_bestroot_plapat)
effect_bestroot_litsap=cbind(va=rownames(effect_bestroot_litsap),effect_bestroot_litsap)
effect_bestroot_woosap=cbind(va=rownames(effect_bestroot_woosap),effect_bestroot_woosap)
effect_bestroot_epiphy=cbind(va=rownames(effect_bestroot_epiphy),effect_bestroot_epiphy)
effect_bestroot_para=cbind(va=rownames(effect_bestroot_para),effect_bestroot_para)


# all the valuables
va=cbind(va=as.character(unique(dm3$va)),code=1:22)%>%data.frame()

effect_best_root=left_join(va,effect_bestroot_acm[,c(1,2,6)], by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_acm[,c(1,2,6)], by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_soilsap[,c(1,2,6)],by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_plapat[,c(1,2,6)],by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_litsap[,c(1,2,6)],by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_woosap[,c(1,2,6)],by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_epiphy[,c(1,2,6)],by="va")
effect_best_root=left_join(effect_best_root,effect_bestroot_para[,c(1,2,6)],by="va")
# remove the code column
effect_best_root=effect_best_root[,-2]

names(effect_best_root)=c("va","acm","pacm", "ecm","pecm","soilsap","psoilsap","plapat","pplapat","litsap", "plitsap","woosap", "pwoosap","epiphy","pepiphy","para","ppara")

rownames(effect_best_root)=effect_best_root$va
effect_best_root_value=effect_best_root[,c(2,4,6,8,10,12,14,16)]
effect_best_root_pva=effect_best_root[,c(3,5,7,9,11,13,15,17)]

effect_best_root_value=melt(as.matrix(effect_best_root_value))
effect_best_root_pva=melt(as.matrix(effect_best_root_pva))

effect_best_root_value=cbind(effect_best_root_value,effect_best_root_pva$value)
names(effect_best_root_value)[4]="pva"

effect_best_root_value$pva[effect_best_root_value$pva<0.001]="***"
effect_best_root_value$pva[effect_best_root_value$pva<0.01&effect_best_root_value$pva>0.001]="**"
effect_best_root_value$pva[effect_best_root_value$pva<0.05&effect_best_root_value$pva>0.01]="*"
effect_best_root_value$pva[effect_best_root_value$pva>0.05]=""

# plots

ggplot(effect_best_root_value, aes(x = X1, y = X2, fill = value)) +
  geom_tile(color = "white", lwd = 1,linetype = 1)+
  geom_text(aes(x =X1, y = X2, label =pva))+
  scale_fill_gradient2("Eeffect size",low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")+
  scale_y_discrete(breaks=as.character(unique(effect_best_root_value$X2)),labels=c("ACM(N=87)","ECM(N=104)","Soil saprotroph(N=104)","Plant pathogen(N=103)","Litter saprotroph(N=104)","Wood saprotroph(N=103)","Epiphyte(N=103)","Parasitic(N=98)"))+
  scale_x_discrete(breaks=as.character(unique(effect_best_root_value$X1)),labels = rev(c("Pla.rich", expression("Root"["cn"]), expression("Root"["c"]),  "d15N","d13C", expression("Root"["fmass"]), expression("Root"["cmass"]),  "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","Tem.seas.","MDR","MAT","Sand","SoilN","pH","Cec","SoilC","Plot.rich","Core.rich")))+
  theme(panel.border = element_rect(fill=NA,size=1,color="black"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15,angle=270,hjust=0,
        color=rev(c("seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1","seagreen1", "royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","peru","purple","purple"))))+
  xlab("")+
  ylab("")+
  ggtitle("Best-fit model")

# combine the plots
P1=ggplotGrob(P1)
P2=ggplotGrob(P2)
plot_grid(P1,P2,ncol=1,labels=c("(a)","(b)",label_x = 0.8),rel_heights = c(1,1.3),label_size = 14)

