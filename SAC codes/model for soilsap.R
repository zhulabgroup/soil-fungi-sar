# model for soilsap
a=list()
for (i in 1:dim(soilsap_c_ranall)[1])
{
  d=data.frame(t(soilsap_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(soilsap_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}
plotID <- rep(soilsap_c_ranall$a1, each = 30)
soilsap_c_ranall_30 <- cbind(plotID, b)
## add the z 

a <- list()
for (i in 1:dim(soilsap_z_ranall)[1])
{
  d <- data.frame(t(soilsap_z_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(soilsap_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

soilsap_z_ranall_30 <- b



soilsap_z_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values
names(soilsap_z_ranall_30)[2]="z"
soilsap_model=cbind(soilsap_c_ranall_30,soilsap_z_ranall_30["z"])

names(soilsap_model)[2]="logc"
soilsap_model=merge(soilsap_model,model_var,by="plotID")

soilsap_model=subset(soilsap_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&richness>0)# only 104 plots from 33 sites
soilsap_model<- cbind(soilsap_model, c = 2.71828^soilsap_model$logc)

soilsap_model_rich=subset(soilsap_model,siteIDD!="GUAN"&z<10&richness>0)# only 104 plots from 33 sites

soilsap_model_rich<- cbind(soilsap_model_rich, c = 2.71828^soilsap_model_rich$logc)

#standardized data
soilsap_model[,c(5:28)]=apply(soilsap_model[,c(5:28)],2,range01)
soilsap_model_rich[,c(5:28)]=apply(soilsap_model_rich[,c(5:28)],2,range01)
# colinearity
ggcorrplot(cor(soilsap_model[,c(2,3,5:27)]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil OC and was excluded,cec was related to soil OC and was removed
#bio4 was related with bio1 and was removed
#coarse (removed) was related with fine and was removed
#cec(removed) and soil OC
#root cn and rootn(removed)
#bold (removed)and soil OC
#bio1 and bio4(removed)
3
# build a model for the soilsap guild,with 104 plots included, soil pH and fundiversity

mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1 +fine+ d15N +d13C +rootc +rootcn +(1 |siteIDD/plotID),data=soilsap_model)

step(mod)

mod_soilsap=lmer(z ~ c + bio12 + richness + funrich + bio1 + (1 | siteIDD/plotID),data=soilsap_model)

plot_model(mod_soilsap)

p3=plot_model(mod_soilsap,axis.labels = c("MAT","Fun.rich","Pla.rich","MAP"),color=c("blue","red"),rm.terms = "c",title="Soilsap. (N=104)")


## when root traits were excluded

mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1  +(1 |siteIDD/plotID),data=soilsap_model_rich)

step(mod)

mod_soilsap_rich=lmer(z ~ c + nitrogen + sand + bio2 + bio18 + bio12 + spei + richness + funrich + (1 | siteIDD/plotID),data=soilsap_model_rich)

plot_model(mod_soilsap_rich)

p30=plot_model(mod_soilsap_rich,axis.labels = c("Fun.rich","Pla.rich","Spei", "MAP","Pre.WQ","MDR","Sand","SoilN"),color=c("blue","red"),rm.terms = "c",title="Soilsap. (N=438)")


pp$data <- transform(pp$data, term = factor(term, levels = c("bio1","bio12", "funrich", "richness")))

print(pp)

##
mod=lmer(z ~ logc+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ rich+funrich +bio1 +fine+ d15N +d13C +rootc +rootcn +(1 |plotID),data=soilsap_model)
soilsap_effect_plotran <- summary(mod)
soilsap_effect_plotran <- data.frame(soilsap_effect_plotran$coefficients)
sig <- round(soilsap_effect_plotran$Pr...t..,3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"

soilsap_effect_plotran <- cbind(soilsap_effect_plotran, sig)

