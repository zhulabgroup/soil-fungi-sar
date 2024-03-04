# model for plapat

range01 <- function(x) 
  {
 return((x - min(x)) / (max(x) - min(x)))
}
  
a <- list()
for (i in 1:dim(plapat_c_ranall)[1])
{
  d <- data.frame(t(plapat_c_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(plapat_c_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(plapat_c_ranall$a1, each = 30)
plapat_c_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values
#add z
a <- list()
for (i in 1:dim(plapat_z_ranall)[1])
{
  d <- data.frame(t(plapat_z_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(plapat_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}


plapat_z_ranall_30 <- b # for each plot,with 30 estimated z values

names(plapat_z_ranall_30)="z"

plapat_model <- cbind(plapat_c_ranall_30, plapat_z_ranall_30["z"])

names(plapat_model)[2] <- "logc"
plapat_model <- merge(plapat_model, model_var, by = "plotID")
plapat_model <- subset(plapat_model, siteIDD != "GUAN" & z < 10 & fine > 0 & rootc > 0 & richness > 0) # only 104 plots from 33 sites

plapat_model <- cbind(plapat_model, c = 2.71828^plapat_model$logc)
plapat_model_rich <- subset(plapat_model, siteIDD != "GUAN" & z < 10 & richness > 0) # only 104 plots from 33 sites


# standardized data
plapat_model[, c(5:28)] <- apply(plapat_model[, c(5:28)], 2, range01)
plapat_model_rich[, c(5:28)] <- apply(plapat_model_rich[, c(5:28)], 2, range01)
# colinearity
ggcorrplot(cor(plapat_model[, c(2, 3, 5:27)]), hc.order = TRUE, type = "lower", lab = TRUE) #
# bold was related with soil OC and was excluded,cec was related to soil OC and was removed
# bio4 was related with bio1 and was removed
# coarse (removed) was related with fine and was removed
# cec(removed) and soil OC
# root cn and rootn(removed)
# bold (removed)and soil OC
# bio1 and bio4(removed)
3
# build a model for the plapat guild,with 104 plots included, rich and d13

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = plapat_model)

step(mod)

mod_plapat=lmer(z ~ c + organicCPercent + bio15 + funrich + fine + (1 | siteIDD/plotID),data=plapat_model)

mod_plapat=lmer(z ~ c + organicCPercent + sand + bio15 + spei + funrich + fine + (1 | siteIDD/plotID),data=plapat_model)

p4=plot_model(mod_plapat,axis.labels = c(expression("Root"["fmass"]),"Fun.rich","Spei","Pre.seas","Sand","SoilC"),colors=c("blue","red"),rm.terms = "c",title="Pla.patho. (N=103)")

# when root traits were excluded
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1  + (1 | siteIDD / plotID), data = plapat_model_rich)

step(mod)

mod_plapat_rich <- lmer(z ~ c + organicCPercent + sand + bio15 + spei + richness + funrich + (1 | siteIDD/plotID), data = plapat_model_rich)

p4=plot_model(mod_plapat_rich,axis.labels = c("Fun.rich", "Pla.rich","Spei","Pre.seas","Sand","SoilC"),colors=c("blue","red"),rm.terms = "c",title="Pla.patho. (N=103)")




##
# #
mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | plotID), data = plapat_model)
plapat_effect_plotran <- summary(mod)
plapat_effect_plotran <- data.frame(plapat_effect_plotran$coefficients)
sig <- round(plapat_effect_plotran$Pr...t..,3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"

plapat_effect_plotran <- cbind(plapat_effect_plotran, sig)

