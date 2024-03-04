# model selection and plot creation when root traits were included
library(sjPlot)
library(ggplot2)
library(effects)
library(lme4)
library(lmerTest)
library(cowplot)
library(styler)

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'Arial',   #To change the font type
          axis.title.size = 2.0,  #To change axis title size
          axis.textsize.x = 1,  #To change x axis text size
          axis.textsize.y = 1,
          title.size = 2,
          title.align= "center")  #To change y axis text size

1. # for abscular mycorrhizal fungi

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+ bio8  + bio12 + bio15 ++ bio18 + spei + richness + funrich + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = acm_model)
step(mod)
mod_acm=lmer(z ~ c + ph + richness + funrich + (1 | siteIDD/plotID),data=acm_model)
p1=plot_model(mod_acm,axis.labels = c("Fun.rich","Pla.rich","pH"),colors="red",rm.terms = "c",title="ACM (N=87)",axis.lim=c(-1, 1),dot.size = 3)


2. # for ectomycorrhizal fungi

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = ecm_model)
step(mod)
mod_ecm=lmer(z ~ c + ph + funrich + d13C + rootcn + (1 | siteIDD/plotID),data=ecm_model)
p2=plot_model(mod_ecm,axis.labels = c(expression("Root"["cn"]),"d13C","Fun.rich","pH"),color=c("blue","red"),rm.terms = "c",title="ACM (N=104)",axis.lim=c(-1, 1),dot.size = 3)


3. # for soil saprotroph
mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1 +fine+ d15N +d13C +rootc +rootcn +(1 |siteIDD/plotID),data=soilsap_model)
step(mod)
mod_soilsap=lmer(z ~ c + bio12 + richness + funrich + bio1 + (1 | siteIDD/plotID),data=soilsap_model)
p3=plot_model(mod_soilsap,axis.labels = c("MAT","Fun.rich","Pla.rich","MAP"),color=c("blue","red"),rm.terms = "c",title="Soilsap. (N=104)",dot.size = 3)

4. # for plant pathogens
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = plapat_model)
step(mod)
mod_plapat=lmer(z ~ c + organicCPercent + bio15 + funrich + fine + (1 | siteIDD/plotID),data=plapat_model)
mod_plapat=lmer(z ~ c + organicCPercent + sand + bio15 + spei + funrich + fine + (1 | siteIDD/plotID),data=plapat_model)
p4=plot_model(mod_plapat,axis.labels = c(expression("Root"["fmass"]),"Fun.rich","Spei","Pre.seas","Sand","SoilC"),colors=c("blue","red"),rm.terms = "c",title="Pla.patho. (N=103)",dot.size = 3)


5.# for litter saprotroph
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = litsap_model)
step(mod)
mod_litsap=lmer(z ~ c + funrich + (1 | siteIDD/plotID),data=litsap_model)
p5=plot_model(mod_litsap,axis.labels = c("Fun.rich"),colors=c("red","red"),rm.terms = "c",title="Lit.sap. (N=104)",axis.lim = c(-1,1),dot.size = 3)

6.# for wood saprotroph
mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1 +fine+ d15N +d13C +rootc +rootcn +(1 |siteIDD/plotID),data=woosap_model)
step(mod)
mod_woosap=lmer(z ~ c + bio15 + funrich + bio1 + rootcn + (1 | siteIDD/plotID),woosap_model)
p6=plot_model(mod_woosap,axis.labels = c(expression("Root"["cn"]),"MAT", "Fun.rich","Pre.seas."),colors=c("blue","red"),rm.terms = "c",title="Woo.sap. (N=103)",dot.size = 3)

#
plot_grid(p1,p2,p3,p4,p5,p6,ncol=2,labels="AUTO")
