# model for ECM
a <- list()
for (i in 1:dim(ecm_c_ranall)[1])
{
  d <- data.frame(t(ecm_c_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(ecm_c_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(ecm_c_ranall$a1, each = 30)

ecm_c_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# add the z value

a <- list()
for (i in 1:dim(ecm_z_ranall)[1])
{
  d <- data.frame(t(ecm_z_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(ecm_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}


ecm_z_ranall_30 <- b # for each plot,with 30 estimated z values
 # for each plot,with 30 estimated z values
names(ecm_z_ranall_30)[1] <- "z"
# cbind the c and z value
# when used the nested model, i loaded the saved data
ecm_model <- cbind(ecm_c_ranall_30, ecm_z_ranall_30["z"])
names(ecm_model)[2] <- "logc"
ecm_model <- merge(ecm_model, model_var, by = "plotID")
ecm_model <- subset(ecm_model, siteIDD != "GUAN" & z < 10 & fine > 0 & rootc > 0 & richness > 0) # only 104 plots from 33 sites
ecm_model_rich <- subset(ecm_model, siteIDD != "GUAN" & z < 10 & richness > 0) # only 104 plots from 33 sites

# head(ecm_model)
# consider the climate and soil data
ecm_model1 <- cbind(ecm_model1, c = 2.71828^ecm_model1$logc)
ecm_model1 <- ecm_model1[, c(1:28)]

ecm_model_rich <- cbind(ecm_model_rich, c = 2.71828^ecm_model_rich$logc)
# standardized data
ecm_model1[, c(5:28)] <- apply(ecm_model1[, c(5:28)], 2, range01)

ecm_model_rich[, c(5:28)] <- apply(ecm_model_rich[, c(5:28)], 2, range01)

# colinearity
ggcorrplot(cor(ecm_model1[, c(5:28)]), hc.order = TRUE, type = "lower", lab = TRUE) #
# bold was related with soil OC and hence was excluded,cec was related to soil C and was removed
# bio4 was related with bio1 and was removed
# coarse (removed) and fine were related
# cec(removed) and soil OC
# root cn and rootn(removed)
# bold (removed)and soil OC
# bio1 and bio4(removed)
3
# build a model for the ecm guild,with 104 plots included, soil pH and fundive

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = ecm_model1)

step(mod)

mod_ecm=lmer(z ~ c + ph + funrich + d13C + rootcn + (1 | siteIDD/plotID),data=ecm_model1)

p2=plot_model(mod_ecm,axis.labels = c(expression("Root"["cn"]),"d13C","Fun.rich","pH"),color=c("blue","red"),rm.terms="logc",title="ECM (N=104)")

plot_model(mod_ecm)

## when root traits were excluded

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD / plotID), data = ecm_model_rich)

step(mod)

mod_ecm_rich=lmer(z ~ c + ph + funrich + (1 | siteIDD/plotID),data=ecm_model_rich)

plot_model(mod_ecm_rich)

p20=plot_model(mod_ecm_rich,axis.labels = c("Fun.rich","pH"),color=c("blue","red"),rm.terms="c",title="ECM (N=438)")



ecm_effect <- summary(mod)
ecm_effect <- data.frame(ecm_effect$coefficients)
sig <- ecm_effect$Pr...t..
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
ecm_effect <- cbind(ecm_effect, sig)

# with plot as the random terms

mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | plotID), data = ecm_model1)
ecm_effect_plotran <- summary(mod)
ecm_effect_plotran <- data.frame(ecm_effect_plotran$coefficients)
sig <- ecm_effect_plotran$Pr...t..
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
ecm_effect_plotran <- cbind(ecm_effect_plotran, sig)

