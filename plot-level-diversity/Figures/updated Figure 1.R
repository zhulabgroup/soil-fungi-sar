# load the required datasets
library(lme4)
library(lmerTest)
library(glmm.hp)
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(patchwork)
library(ggcorrplot)
library(glmm.hp)

load("/model_data.RData")
load("/core_rich.RData") # core-level fungal richness mean for each plot

1. # climate and soil model

data_corich <- merge(model_data[, 1:20], core_rich, by = "plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10)
data_corich[, 5:21] <- apply(data_corich[, 5:21], 2, range01) %>% data.frame()

ggcorrplot(cor(data_corich[, c(5:18,20:21)]), hc.order = TRUE, type = "lower", lab = TRUE)

#cec-soilC
#bold-soilC
# bold and cec were excluded from the model

mod <- lmer(z ~ corich + funrich + organicCPercent + ph + nitrogen + sand + bio1 + bio2 + bio4 + bio8 + bio12 + bio15 + bio18 + spei + (1 | siteIDD / plotID), data = data_corich)

effect_CS <- summary(mod)
effect_CS <- effect_CS$coefficients
effect_CS <- data.frame(effect_CS)[2:dim(effect_CS)[1], ]

# best model
step(mod)

mod_best=lmer(z ~ corich + bio4 + bio12 + bio15 + (1 | siteIDD/plotID),data=data_corich)

effect_best_CS <- summary(mod_best)
effect_best_CS <- effect_best_CS$coefficients
effect_best_CS <- data.frame(effect_best_CS)[2:dim(effect_best_CS)[1], ]


# the variables were categorized into different types

# construct the model


data_corich$climate=data_corich$bio4*data_corich$bio12*data_corich$bio15

data_corich$myco=data_corich$corich


mod=lmer(z~climate+myco+ (1 | siteIDD/plotID),data=data_corich)

varpcs=glmm.hp(mod,commonality=TRUE)

varpcs=data.frame(va=rownames(varpcs$commonality.analysis),varpcs$commonality.analysis)

varpcs$va=gsub(" ","",varpcs$va)

p1=ggplot(data=subset(varpcs,va!="Total"&Fractions>0),aes(x = 1, y = Fractions, fill = va),alpha=0.5)+
  geom_bar(stat = "identity", pch=21,color="black",position = "fill",width=0.6)+
  theme(legend.position = c(0.5,0.8), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  ylab("Contribution to the explained variance")+
  xlab("")+
  scale_fill_manual("Components",breaks=c("Uniquetoclimate","Uniquetomyco"),
                    labels=c("Clima.","Fung."),values=c("mediumpurple","tomato"))+
  ylim(0,1.4)



# plot the variance partitioning results
plot(p1,color=1:2,labels=c(1:2))
dk=data.frame(a=c(1:4),af=c(5:8))

pv1=ggplot(dk)+
  geom_point(x=0.5,y=0.75,size=80,pch=21,alpha=0.2,color="black",fill="gray")+
  geom_point(x=1,y=0.75,size=80,pch=21,alpha=0.2,color="black",fill="pink")+
  xlim(0,1.5)+
  ylim(0,1.5)+
  theme( panel.background = element_blank(), 
         panel.border = element_rect(fill=NA,size=1,color="black"),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank(),
         plot.title=element_text(hjust=0.5,face="bold",size=18)
        )+
  annotate("text",x=0.25,y=0.75,label="0.03",size=6)+
annotate("text",x=0.75,y=0.75,label="0.01",size=6)+
annotate("text",x=1.1,y=0.75,label="0.00",size=6)+
  annotate("text",x=1.1,y=1.4,label="Community",size=6)+
  annotate("text",x=0.25,y=1.4,label="Climate",size=6)+
  xlab("")+
  ylab("")+
  ggtitle("Clima.+Soil\n (N=483)")
  
  
  

###
 ggplot() +
  geom_point(
    data = effect_CS, aes(x = Estimate, y = 1:dim(effect_CS)[1]),
    color = rev(c("royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), size = 3
  ) +
  geom_segment(data = effect_CS, size = 0.8, color = rev(c("royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_CS$Estimate - 1.96 * effect_CS$Std..Error, y = 1:dim(effect_CS)[1], xend = effect_CS$Estimate + 1.96 * effect_CS$Std..Error, yend = 1:dim(effect_CS)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:14, labels = rev(c("Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich"))) +
  theme(axis.text.y = element_text(colour = rev(c("royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")))) +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = 0.28, y = 11, label = "", size = 10) +
  annotate("text", x = 0.182048, y = 12, label = "*", size = 10) +
  annotate("text", x = -0.1592048, y = 1, label = "*", size = 10) +
  annotate("text", x = 0.425192048, y = 9, label = "**", size = 10) +
  ggtitle("Clima.+Soil\n (N=483)")

2. # when climate, soil and plant richness were included

data_corich <- merge(model_data[, 1:20], core_rich, by = "plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10 & richness > 0)
data_corich[, 5:21] <- apply(data_corich[, 5:21], 2, range01) %>% data.frame()
ggcorrplot(cor(data_corich[, c(5:21)]), hc.order = TRUE, type = "lower", lab = TRUE)

# only bold and soilC were excluded

mod <- lmer(z ~ corich + funrich + organicCPercent + ph + nitrogen + sand + bio1 + bio2 + bio4 + bio8 + bio12 + bio15 + bio18 + spei + richness + (1 | siteIDD / plotID), data = data_corich)
effect_CSP <- summary(mod)
effect_CSP <- effect_CSP$coefficients
effect_CSP <- data.frame(effect_CSP)[2:dim(effect_CSP)[1], ]


#best model
step(mod)

mod_best=lmer(z ~ corich + bio4 + bio12 + bio15 + richness + (1 | siteIDD/plotID),data=data_corich)

effect_best_CSP <- summary(mod_best)
effect_best_CSP <- effect_best_CSP$coefficients
effect_best_CSP <- data.frame(effect_best_CSP)[2:dim(effect_best_CSP)[1], ]

#


data_corich$climate=data_corich$bio4*data_corich$bio12*data_corich$bio15
data_corich$plant=data_corich$richness
data_corich$myco=data_corich$corich
# fitting model

mod=lmer(z~climate+myco+plant+ (1 | siteIDD/plotID),data=data_corich)

varpcsp=glmm.hp(mod,commonality=TRUE)
varpcsp=data.frame(va=rownames(varpcsp$commonality.analysis),varpcsp$commonality.analysis)
varpcsp$va=gsub(" ","",varpcsp$va)

p2=ggplot(data=subset(varpcsp,va!="Total"&Fractions>0),aes(x = 1, y = Fractions, fill = va),alpha=0.5)+
  geom_bar(stat = "identity", pch=21,color="black",position = "fill",width=0.16)+
  scale_fill_manual("Components",breaks=c("Uniquetoclimate","Uniquetomyco","Uniquetoplant","Commontoclimate,andplant",
                                  "Commontoclimate,myco,andplant"),
  labels=c("Clima.","Fung.","Plant","Clima.+Plant","Clima.+Fung.\n+Plant"),values=c("mediumpurple","tomato","royalblue1","tan","purple"))+
  theme(legend.position = c(0.5,0.8), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 5, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  ylab("Contribution to the explained variance")+
  xlab("")+
  ylim(0,1.5)





# create the graph of different components


pv2=ggplot(dk)+
  geom_point(x=0.5,y=0.75,size=80,pch=21,alpha=0.2,color="black",fill="lightblue1")+
  geom_point(x=1,y=0.75,size=80,pch=21,alpha=0.2,color="black",fill="sandybrown")+
  geom_point(x=0.75,y=0.3,size=80,pch=21,alpha=0.2,color="black",fill="seagreen")+
  
  theme( panel.background = element_blank(), 
         panel.border = element_rect(fill=NA,size=1,color="black"),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20),
         #axis.text.x = element_blank(),
         #axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank(),
         plot.title=element_text(hjust=0.5,face="bold",size=18)
  )+
  annotate("text",x=0.25,y=0.75,label="0.03",size=6)+
  annotate("text",x=0.75,y=0.75,label="-0.01",size=6)+
  annotate("text",x=1.1,y=0.75,label="0.02",size=6)+
  annotate("text",x=1.1,y=1.24,label="Community",size=6)+
  annotate("text",x=0.25,y=1.24,label="Climate",size=6)+
  annotate("text",x=0.75,y=0.25,label="0.02",size=6)+
  annotate("text",x=0.75,y=-0.15,label="Plant",size=6)+
  annotate("text",x=0.55,y=0.5,label="0.01",size=6)+
  annotate("text",x=0.85,y=0.5,label="-0.01",size=6)+
  xlab("")+
  ylab("")+
  ggtitle("Clima.+Soil+Pla.\n (N=438)")+
  xlim(0,1.5)+
  ylim(-0.2,1.5)



ggplot() +
  geom_point(
    data = effect_CSP, aes(x = Estimate, y = 1:dim(effect_CSP)[1]),
    color = rev(c("seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), size = 3
  ) +
  geom_segment(data = effect_CSP, size = 0.8, color = rev(c("seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_CSP$Estimate - 1.96 * effect_CSP$Std..Error, y = 1:dim(effect_CSP)[1], xend = effect_CSP$Estimate + 1.96 * effect_CSP$Std..Error, yend = 1:dim(effect_CSP)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:15, labels = rev(c("Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich"))) +
  theme(axis.text.y = element_text(colour = rev(c("seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")))) +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = 0.28, y = 11, label = "", size = 10) +
  annotate("text", x = 0.182048, y = 12, label = "*", size = 10) +
  annotate("text", x = -0.1592048, y = 1, label = "*", size = 10) +
  annotate("text", x = 0.425192048, y = 9, label = "**", size = 10) +
  annotate("text", x = 0.425192048, y = 15, label = "***", size = 10) +
  ggtitle("Clima.+Soil+Pla.\n (N=438)")

3. # when climate, soil, plant richness and root traits were included

data_corich <- merge(model_data[, 1:27], core_rich, by = "plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10 & richness > 0 & fine > 0 & rootc > 0)
data_corich[, 5:27] <- apply(data_corich[, 5:27], 2, range01) %>% data.frame()

ggcorrplot(cor(data_corich[, c(5:28)]), hc.order = TRUE, type = "lower", lab = TRUE)
#coarse-fine
#cec and soilC
# rootn-rootcn
#bio1-bio4
#bold-soilC


mod <- lmer(z ~ corich + funrich + organicCPercent + ph + nitrogen + sand + bio1 + bio2 + bio4 + bio8 + bio12 + bio15 + bio18 + spei + richness + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = data_corich)

mod <- lmer(z ~ corich + funrich + organicCPercent + ph + nitrogen + sand + bio1 + bio2 + bio8 + bio12 + bio15 + bio18 + spei + richness + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = data_corich)

effect_CSPR <- summary(mod)
effect_CSPR <- effect_CSPR$coefficients
effect_CSPR <- data.frame(effect_CSPR)[2:dim(effect_CSPR)[1], ]



mod_best=lmer(z ~ corich + funrich + bio1 + richness + rootcn + (1 | siteIDD/plotID),data=data_corich)

effect_best_CSPR <- summary(mod_best)
effect_best_CSPR <- effect_best_CSPR$coefficients
effect_best_CSPR <- data.frame(effect_best_CSPR)[2:dim(effect_best_CSPR)[1], ]

# three different types of predictors

data_corich$climate=data_corich$bio1
data_corich$plant=data_corich$richness*data_corich$rootn
data_corich$myco=data_corich$corich*data_corich$funrich


mod=lmer(z ~climate + myco + plant + (1 | siteIDD/plotID),data=data_corich)

varpcsr=glmm.hp(mod,commonality=TRUE)

varpcsr=data.frame(va=rownames(varpcsr$commonality.analysis),varpcsr$commonality.analysis)

varpcsr$va=gsub(" ","",varpcsr$va)

p3=ggplot(data=subset(varpcsr,va!="Total"&Fractions>0),aes(x = 1, y = Fractions, fill = va),alpha=0.5)+
  geom_bar(stat = "identity", pch=21,color="black",position = "fill",width=0.6)+
  theme(legend.position = c(0.5,0.8), 
        legend.key.size = unit(0.15, "inches"),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  ylab("Contribution to the explained variance")+
  xlab("")+
  ylim(0,1.5)+
  scale_fill_manual("Components",breaks=c("Commontoclimate,andplant","Uniquetoclimate","Uniquetomyco","Uniquetoplant"),
                    labels=c("Clima.+Plant","Clima.","Fung.","Plant"),values=c("tan","mediumpurple","tomato","royalblue1"))







# partitioning plots


pv3=ggplot(dk)+
  geom_point(x=0.5,y=0.75,size=80,pch=21,alpha=0.2,color="black",fill="lightblue1")+
  geom_point(x=1,y=0.75,size=80,pch=21,alpha=0.2,color="black",fill="sandybrown")+
  geom_point(x=0.75,y=0.3,size=80,pch=21,alpha=0.2,color="black",fill="seagreen")+
  
  theme( panel.background = element_blank(), 
         panel.border = element_rect(fill=NA,size=1,color="black"),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank(),
         plot.title=element_text(hjust=0.5,face="bold",size=18)
  )+
  annotate("text",x=0.25,y=0.75,label="0.002",size=6)+
  annotate("text",x=0.75,y=0.75,label="0.0001",size=6)+
  annotate("text",x=1.1,y=0.75,label="-0.0002",size=6)+
  annotate("text",x=1.1,y=1.24,label="Community",size=6)+
  annotate("text",x=0.25,y=1.24,label="Climate",size=6)+
  annotate("text",x=0.75,y=0.25,label="0.02",size=6)+
  annotate("text",x=0.75,y=-0.15,label="Plant",size=6)+
  annotate("text",x=0.55,y=0.5,label="-0.001",size=6)+
  annotate("text",x=0.85,y=0.5,label="0.002",size=6)+
  xlab("")+
  ylab("")+
  ggtitle("Clima.+Soil+Pla.\n (N=438)")+
  xlim(0,1.5)+
  ylim(-0.2,1.5)


p03 <- ggplot() +
  geom_point(
    data = effect_CSPR, aes(x = Estimate, y = 1:dim(effect_CSPR)[1]),
    color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), size = 3
  ) +
  geom_segment(data = effect_CSPR, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_CSPR$Estimate - 1.96 * effect_CSPR$Std..Error, y = 1:dim(effect_CSPR)[1], xend = effect_CSPR$Estimate + 1.96 * effect_CSPR$Std..Error, yend = 1:dim(effect_CSPR)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich"))) +
  theme(axis.text.y = element_text(colour = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")))) +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = 0.2592048, y = 2, label = "*", size = 10) +
  annotate("text", x = 0.325192048, y = 15, label = "***", size = 10) +
  ggtitle("Clima.+Soil+Pla.+\n Root (N=104)")
# create a df to combine all the three figures

effect_ALL <- merge(effect_CSPR[, c("va", "Estimate", "Std..Error")], effect_CS[, c("va", "Estimate", "Std..Error")], by = "va", all.x = TRUE, sort = FALSE)
effect_ALL <- merge(effect_ALL, effect_CSP[, c("va", "Estimate", "Std..Error")], by = "va", all.x = TRUE, sort = FALSE)
names(effect_ALL) <- c("va", "es", "sd", "es1", "sd1", "es2", "sd2")
effect_ALL[is.na(effect_ALL)] <- 0

# when all the variables were included
# for the climate and soil model

 p01=ggplot() +
  geom_point(
    data = effect_ALL, aes(x = es1, y = 1:dim(effect_ALL)[1]),
    color = rev(c("gray", "gray", "gray", "gray", "gray", "gray", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), size = 3
  ) +
  geom_segment(data = effect_ALL, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL$es1 - 1.96 * effect_ALL$sd1, y = 1:dim(effect_ALL)[1], xend = effect_ALL$es1 + 1.96 * effect_ALL$sd1, yend = 1:dim(effect_ALL)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich"))) +
  theme(axis.text.y = element_text(colour = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")))) +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = 0.3592048, y = 9, label = "*", size = 10) +
  annotate("text", x = -0.15325192048, y = 1, label = "*", size = 10) +
  annotate("text", x = 0.182053025192048, y = 12, label = "**", size = 10) +
  ggtitle("Clima.+Soil \n (N=483)") +
  xlab("")

# for the climate,soil and plant model
p002 <- ggplot() +
  geom_point(
    data = effect_ALL, aes(x = es2, y = 1:dim(effect_ALL)[1]),
    color = rev(c("gray", "gray", "gray", "gray", "gray", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), size = 3
  ) +
  geom_segment(data = effect_ALL, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL$es2 - 1.96 * effect_ALL$sd2, y = 1:dim(effect_ALL)[1], xend = effect_ALL$es2 + 1.96 * effect_ALL$sd2, yend = 1:dim(effect_ALL)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich"))) +
  theme(axis.text.y = element_text(colour = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")))) +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = 0.2592048, y = 2, label = "", size = 10) +
  annotate("text", x = 0.3025192048, y = 15, label = "***", size = 10) +
  annotate("text", x = -0.185325192048, y = 1, label = "**", size = 10) +
  annotate("text", x = 0.3029851325192048, y = 9, label = "*", size = 10) +
  annotate("text", x = 0.15325192048, y = 12, label = "*", size = 10) +
  annotate("text", x = 0.275325192048, y = 11, label = "**", size = 10) +
  ggtitle("Clima.+Soil+Pla.\n (N=438)") +
  xlab("Effect size and 95% CI")


# for the climate, soil and plant model

p003 <- ggplot() +
  geom_point(
    data = effect_ALL, aes(x = es, y = 1:dim(effect_ALL)[1]),
    color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), size = 3
  ) +
  geom_segment(data = effect_ALL, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL$es - 1.96 * effect_ALL$sd, y = 1:dim(effect_ALL)[1], xend = effect_ALL$es + 1.96 * effect_ALL$sd, yend = 1:dim(effect_ALL)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich"))) +
  theme(axis.text.y = element_text(colour = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")))) +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = 0.2592048, y = 2, label = "*", size = 10) +
  annotate("text", x = 0.325192048, y = 15, label = "***", size = 10) +
  ggtitle("Clima.+Soil+Pla.+\n Root (N=104)") +
  xlab("")
# for the best fitted model

effect_best_CS=cbind(va=rownames(effect_best_CS),effect_best_CS)
effect_best_CSP=cbind(va=rownames(effect_best_CSP),effect_best_CSP)
effect_best_CSPR=cbind(va=rownames(effect_best_CSPR),effect_best_CSP)

effect_ALL_best= left_join(effect_ALL[,1:2],effect_best_CS[,c(1,2,3)],by="va")
effect_ALL_best= left_join(effect_ALL_best,effect_best_CSP[,c(1,2,3)],by="va")
effect_ALL_best= left_join(effect_ALL_best,effect_best_CSPR[,c(1,2,3)],by="va")

names(effect_ALL_best)=c("va","es","esCS","sdCS","esCSP","sdCSP","esCSPR","SDCSPR")
effect_ALL_best[is.na(effect_ALL_best)]=0

#create the plots
# for the climate and soil model

ggplot() +
  geom_point(
    data = effect_ALL_best, aes(x = esCS, y = 1:dim(effect_ALL_best)[1]),
    color = rev(c("gray", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "royalblue1", "royalblue1", "gray", "royalblue1", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "purple")), size = 4
  ) +
  geom_segment(data = effect_ALL_best, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL_best$esCS - 1.96 * effect_ALL_best$sdCS, y = 1:dim(effect_ALL_best)[1], xend = effect_ALL_best$esCS + 1.96 * effect_ALL_best$sdCS, yend = 1:dim(effect_ALL_best)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich")))  +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.text.y = element_text(size = 15, hjust=0,color = "black",margin = margin(l = 20, r = -90)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = -0.20592048, y = 1, label = "**", size = 10) +
  annotate("text", x = 0.325192048, y = 15, label = "", size = 10) +
  annotate("text", x = 0.2825192048, y = 9, label = "*", size = 10) +
  annotate("text", x = 0.325192048, y = 11, label = "***", size = 10) +
  annotate("text", x = 0.225192048, y = 12, label = "*", size = 10) +
  
  ggtitle("Clima.+Soil (N=483)") +
  xlab("")
# for the climate, soil and plant model

ggplot() +
  geom_point(
    data = effect_ALL_best, aes(x = esCSP, y = 1:dim(effect_ALL_best)[1]),
    color = rev(c("gray", "gray", "gray", "gray", "gray", "seagreen1", "gray", "gray", "royalblue1", "royalblue1", "gray", "royalblue1", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "purple")), size = 4
  ) +
  geom_segment(data = effect_ALL_best, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL_best$esCSP - 1.96 * effect_ALL_best$sdCSP, y = 1:dim(effect_ALL_best)[1], xend = effect_ALL_best$esCSP + 1.96 * effect_ALL_best$sdCSP, yend = 1:dim(effect_ALL_best)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich")))  +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15,hjust=0, color = "black",margin = margin(l = 20, r = -90)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = -0.240592048, y = 1, label = "***", size = 10) +
  annotate("text", x = 0.325192048, y = 15, label = "***", size = 10) +
  annotate("text", x = 0.2825192048, y = 9, label = "**", size = 10) +
  annotate("text", x = 0.325192048, y = 11, label = "***", size = 10) +
  annotate("text", x = 0.225192048, y = 12, label = "*", size = 10) +
  annotate("text", x = 0.225192048, y = 12, label = "*", size = 10) +
  ggtitle("Clima.+Soil+Pla.(N=438)") +
  xlab("")
# for the climate soil plant and root model

ggplot() +
  geom_point(
    data = effect_ALL_best, aes(x = esCSPR, y = 1:dim(effect_ALL_best)[1]),
    color = rev(c("seagreen1", "gray", "gray", "gray", "gray", "seagreen1", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "royalblue1", "gray", "gray", "gray", "gray", "purple", "purple")), size = 4
  ) +
  geom_segment(data = effect_ALL_best, size = 0.8, color = rev(c("seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL_best$esCSPR - 1.96 * effect_ALL_best$SDCSPR, y = 1:dim(effect_ALL_best)[1], xend = effect_ALL_best$esCSPR + 1.96 * effect_ALL_best$SDCSPR, yend = 1:dim(effect_ALL_best)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = 1:20, labels = rev(c(expression("Root"["cn"]), expression("Root"["c"]), "d13C", "d15N", expression("Root"["fmass"]), "Pla.rich", "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich", "Core.rich")))  +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    axis.ticks.length=unit(-0.25, "cm"),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15, hjust=0,color = "black",margin = margin(l = 20, r = -90)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = -0.230592048, y = 1, label = "***", size = 10) +
  annotate("text", x = 0.1325192048, y = 15, label = "*", size = 10) +
  annotate("text", x = 0.2825192048, y = 20, label = "***", size = 10) +
  annotate("text", x = 0.325192048, y = 11, label = "", size = 10) +
  annotate("text", x = 0.225192048, y = 2, label = "**", size = 10) +
  annotate("text", x = 0.2825192048, y = 7, label = "***", size = 10) +
  ggtitle("Clima.+Soil+Pla.+Root(N=104)") +
  xlab("")

cowplot::plot_grid(pp1,pp2,pp3,ncol=3,labels=c("(a)","(b)","(c)"),label_x = 0.2)


#
p1=ggplotGrob(p1)
pp1=ggplotGrob(pp1)
p1$heights=pp1$heights
p2=ggplotGrob(p2)
pp2=ggplotGrob(pp2)
p3=ggplotGrob(p3)
pp3=ggplotGrob(pp3)
p2$heights=pp2$heights
p3$heights=pp3$heights

# when root traits were not considered
# need to exclude the root and plant variables

effect_ALL_best1=effect_ALL_best[!(effect_ALL_best$va%in%c("richness",   "fine" , "d15N" , "d13C", "rootc"  , "rootcn" )),]

p01=ggplot() +
  geom_point(
    data = effect_ALL_best1, aes(x = esCS, y = 1:dim(effect_ALL_best1)[1]),
    color = rev(c("gray", "gray", "royalblue1", "royalblue1", "gray", "royalblue1", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "purple")), size = 5
  ) +
  geom_segment(data = effect_ALL_best1, size = 0.8, color = rev(c( "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL_best1$esCS - 1.96 * effect_ALL_best1$sdCS, y = 1:dim(effect_ALL_best1)[1], xend = effect_ALL_best1$esCS + 1.96 * effect_ALL_best1$sdCS, yend = 1:dim(effect_ALL_best1)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  
  scale_y_continuous(breaks = 1:14, labels = rev(c( "Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas.", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich.", "Core.rich.")))  +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.text.y = element_text(size = 15, hjust=0,color = "black",margin = margin(l = 20, r = -90)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = -0.20592048, y = 1, label = "**", size = 10) +
  annotate("text", x = 0.325192048, y = 15, label = "", size = 10) +
  annotate("text", x = 0.2825192048, y = 9, label = "*", size = 10) +
  annotate("text", x = 0.325192048, y = 11, label = "***", size = 10) +
  annotate("text", x = 0.152025192048, y = 12, label = "*", size = 10) +
  
  ggtitle("(a)Clima.+Soil (N=483)") +
  xlab("")+
  xlab("Effect size ± 95% CI")
# 
effect_ALL_best2=effect_ALL_best[!(effect_ALL_best$va%in%c(  "fine" , "d15N" , "d13C", "rootc"  , "rootcn" )),]

p02=ggplot() +
  geom_point(
    data = effect_ALL_best2, aes(x = esCSP, y = 1:dim(effect_ALL_best2)[1]),
    color = rev(c("seagreen1","gray", "gray", "royalblue1", "royalblue1", "gray", "royalblue1", "gray", "gray", "gray", "gray", "gray", "gray", "gray", "purple")), size = 5
  ) +
  geom_segment(data = effect_ALL_best2, size = 0.8, color = rev(c("seagreen1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "royalblue1", "peru", "peru", "peru", "peru", "purple", "purple")), aes(x = effect_ALL_best2$esCSP - 1.96 * effect_ALL_best2$sdCSP, y = 1:dim(effect_ALL_best2)[1], xend = effect_ALL_best2$esCSP + 1.96 * effect_ALL_best2$sdCSP, yend = 1:dim(effect_ALL_best2)[1])) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  
  scale_y_continuous(breaks = 1:15, labels = rev(c("Pla.rich." ,"Spei", "Pre.WQ", "Pre.seas.", "MAP", "MTWQ", "Tem.seas.", "MDR", "MAT", "Sand", "SoilN", "pH", "SoilC", "Plot.rich.", "Core.rich.")))  +
  theme(
    legend.key.size = unit(0.18, "inches"),
    legend.position = c(0.4, 0.85),
    legend.text.align = 0, panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.ticks.length=unit(-0.25, "cm"),
    axis.text.y = element_text(size = 15, hjust=0,color = "black",margin = margin(l = 20, r = -90)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 15)
  ) +
  xlim(-0.8, 0.5) +
  ylab("") +
  annotate("text", x = -0.20592048, y = 1, label = "**", size = 10) +
  annotate("text", x = 0.2625192048, y = 15, label = "***", size = 10) +
  annotate("text", x = 0.25825192048, y = 9, label = "*", size = 10) +
  annotate("text", x = 0.325192048, y = 11, label = "***", size = 10) +
  annotate("text", x = 0.152025192048, y = 12, label = "*", size = 10) +
  
  ggtitle("(b)Clima.+Soil+Plant (N=438)") +
  xlab("Effect size ± 95% CI")
# the variance 

#
p1=ggplotGrob(p1)
p2=ggplotGrob(p2)
p01=ggplotGrob(p01)
p02=ggplotGrob(p02)

p1$heights=p01$heights
p2$heights=p02$heights

plot_grid(p1,p01,p2,p02,  ncol = 4,rel_widths = c(0.4,1,0.4,1))

