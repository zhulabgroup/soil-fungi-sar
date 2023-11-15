# for specific functional guilds
# for the ACM functional guild
# get the explanatory variables from the full data based on the plotID
head(model_data)
model_var <- model_data[, c(2, 3, 6:28)]
model_var <- unique(model_var)

a <- list()
for (i in 1:dim(acm_c_ranall)[1])
{
  d <- data.frame(t(acm_c_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(acm_c_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(acm_c_ranall$a1, each = 30)

acm_c_ranall_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values

# add the z value

a <- list()
for (i in 1:dim(acm_z_ranall)[1])
{
  d <- data.frame(t(acm_z_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(acm_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}



acm_z_ranall_30 <- b # for each plot,with 30 estimated z values

acm_model <- cbind(acm_c_ranall_30, acm_z_ranall_30["id"])
names(acm_model) <- c("plotID", "logc", "z")
acm_model <- merge(acm_model, model_var, by = "plotID")
acm_model <- subset(acm_model, siteIDD != "GUAN" & z < 10 & fine > 0 & rootc > 0 & rich > 0) # only 87 plots from 33 sites
# head(acm_model)
# consider the climate and soil data


acm_model1 <- cbind(acm_model, c = 2.71828^acm_model$logc)

range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}


# standardized data
acm_model1[, c(2, 5:27)] <- apply(acm_model1[, c(2, 5:27)], 2, range01)

# testing colinearity
ggcorrplot(cor(acm_model1[, c(2, 5:27)]), hc.order = TRUE, type = "lower", lab = TRUE) #

# bold was related with soil c and hence was excluded
# root n (excluded) and root cn were correlated
# build a model for the acm guild,355 plots
# decide to remove cec to reduce model complexity

mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + bio1 + (1 | siteIDD / plotID), data = acm_model1)
acm_effect <- summary(mod)
acm_effect <- data.frame(acm_effect$coefficients)
sig <- acm_effect$Pr...t..
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
acm_effect <- cbind(acm_effect, sig)

## use plotid as the only random effect

mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + bio1 + (1 | plotID), data = acm_model1)
acm_effect_plotran <- summary(mod)
acm_effect_plotran <- data.frame(acm_effect_plotran$coefficients)
sig <- acm_effect_plotran$Pr...t..
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
acm_effect_plotran <- cbind(acm_effect_plotran, sig)




# for the random slope mode, takes a long time to compute
mod_slop_acm <- lmer(z ~ logc + organicCPercent + ph + nitrogen + cec + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + bio1 + (0 + logc + organicCPercent + ph + nitrogen + cec + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + rich + funrich + bio1 + fine + d15N + d13C + rootc + rootcn + bio1 | siteIDD / plotID), data = acm_model2)
