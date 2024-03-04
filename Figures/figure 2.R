# codes that create figure 2. these models don't include root traits.

1.# for abscular mycorrhizal fungi
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+ bio8  + bio12 + bio15 + bio18 + spei + richness + funrich  + (1 | siteIDD / plotID), data = acm_model_rich)
step(mod)
mod_acm_rich=lmer(z ~ c + organicCPercent + bio15 + funrich + (1 | siteIDD/plotID),data=acm_model_rich)

p1=plot_model(mod_acm_rich,axis.labels = c("Fun.rich","Pre.seas.","SoilC"),colors=c("blue","red"),rm.terms="c",title="ACM (N=319)",axis.lim=c(-1, 1))

2.# for ectomycorrhizal fungi
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = ecm_model_rich)
step(mod)
mod_ecm_rich=lmer(z ~ c + ph + funrich + (1 | siteIDD/plotID),data=ecm_model_rich)

p2=plot_model(mod_ecm_rich,axis.labels = c("Fun.rich","pH"),color=c("blue","red"),rm.terms="c",title="ECM (N=438)")

3.# for soil sapprotroph

mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1  +(1 |siteIDD/plotID),data=soilsap_model_rich)
step(mod)
mod_soilsap_rich=lmer(z ~ c + nitrogen + sand + bio2 + bio18 + bio12 + spei + richness + funrich + (1 | siteIDD/plotID),data=soilsap_model_rich)

p3=plot_model(mod_soilsap_rich,axis.labels = c("Fun.rich","Pla.rich","Spei", "MAP","Pre.WQ","MDR","Sand","SoilN"),color=c("blue","red"),rm.terms = "c",title="Soilsap. (N=438)")

4. # for plant pathogens

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1  + (1 | siteIDD / plotID), data = plapat_model_rich)
step(mod)
mod_plapat_rich <- lmer(z ~ c + organicCPercent + sand + bio15 + spei + richness + funrich + (1 | siteIDD/plotID), data = plapat_model_rich)
p4=plot_model(mod_plapat_rich,axis.labels = c("Fun.rich", "Pla.rich","Spei","Pre.seas","Sand","SoilC"),colors=c("blue","red"),rm.terms = "c",title="Pla.patho. (N=427)")

5. # for litter saprotroph
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 +  (1 | siteIDD / plotID), data = litsap_model_rich)
step(mod)
mod_litsap_rich <- mod <- lmer(z ~   c + sand + bio12 + bio15 + spei + richness + funrich + (1 | siteIDD / plotID), data = litsap_model_rich)
p5=plot_model(mod_litsap_rich,axis.labels = c("Fun.rich","Pla.rich","Spei","Pre.seas.","MAP","Sand"),colors=c("blue","red"),rm.terms = "c",title="Lit.sap. (N=438)",axis.lim = c(-1,1))

6. # for wood saprotroph
mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1  +(1 |siteIDD/plotID),data=woosap_model_rich)
step(mod)
mod_woosap_rich=lmer(z ~ c + sand + bio2 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD/plotID),woosap_model_rich)
p6=plot_model(mod_woosap_rich,axis.labels = c("MAT", "Fun.rich","Pla.rich","Spei","MAP", "MDR","Sand"),colors=c("blue","red"),rm.terms = "c",title="Woo.sap. (N=427)")

