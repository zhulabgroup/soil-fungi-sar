# model for litsap
a <- list()
for (i in 1:dim(litsap_c_ranall)[1])
{
  d <- data.frame(t(litsap_c_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(litsap_c_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(litsap_c_ranall$a1, each = 30)

litsap_c_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

#add the z values
a <- list()
for (i in 1:dim(litsap_z_ranall)[1])
{
  d <- data.frame(t(litsap_z_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(litsap_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}
litsap_z_ranall_30=b
names(litsap_z_ranall_30)="z"

litsap_model <- cbind(litsap_c_ranall_30, litsap_z_ranall_30["z"])
names(litsap_model)[2] <- "logc"
litsap_model <- merge(litsap_model, model_var, by = "plotID")
litsap_model <- subset(litsap_model, siteIDD != "GUAN" & z < 10 & fine > 0 & rootc > 0 & richness > 0) # only 104 plots from 33 sites

litsap_model  <- cbind(litsap_model, c = 2.71828^litsap_model $logc)

# standardized data
litsap_model[, c( 5:28)] <- apply(litsap_model[, c(5:28)], 2, range01)
# colinearity
ggcorrplot(cor(litsap_model[, c(2, 3, 5:27)]), hc.order = TRUE, type = "lower", lab = TRUE) #
# bold was related with soil OC and was excluded,cec was related to soil OC and was removed
# bio4 was related with bio1 and was removed
# coarse (removed) was related with fine and was removed
# cec(removed) and soil OC
# root cn and rootn(removed)
# bold (removed)and soil OC
# bio1 and bio4(removed)
3
# build a model for the litsap guild,with 104 plots included, rich and d13

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + richness + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = litsap_model)

step(mod)

mod_litsap=lmer(z ~ c + funrich + (1 | siteIDD/plotID),data=litsap_model)

p5=plot_model(mod_litsap,axis.labels = c("Fun.rich"),colors=c("blue","red"),rm.terms = "c",title="Lit.sap. (N=104)",axis.lim = c(-1,1))


##
mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + (1 | plotID), data = litsap_model)
