# for specific functional guilds
# for the ACM functional guild
# get the explanatory variables from the full data based on the plotID
library(ggcorrplot)
head(model_data)
model_var <- model_data[, c(1,2,5:27)]
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
acm_model <- subset(acm_model, siteIDD != "GUAN" & z < 10 & fine > 0 & rootc > 0 & richness > 0) # only 87 plots from 33 sites
# does not include root traits but plant richness
acm_model_rich <- subset(acm_model, siteIDD != "GUAN" & z < 10 &richness>0) # only 87 plots from 33 sites

# head(acm_model)
# consider the climate and soil data


acm_model <- cbind(acm_model, c = 2.71828^acm_model$logc)

acm_model_rich <- cbind(acm_model_rich, c = 2.71828^acm_model_rich$logc)

range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}


# standardize the data

acm_model[, c(5:28)] <- apply(acm_model[, c(5:28)], 2, range01)

acm_model_rich[, c(5:28)] <- apply(acm_model_rich[, c(5:28)], 2, range01)

# testing colinearity
ggcorrplot(cor(acm_model1[, c(5:28)]), hc.order = TRUE, type = "lower", lab = TRUE) #

# bold was related with soil c and hence was excluded
# root n (excluded) and root cn were correlated
# build a model for the acm guild,355 plots
# decide to remove cec to reduce model complexity

# plotid and site id were nested as the random effects

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+ bio8  + bio12 + bio15 ++ bio18 + spei + richness + funrich + fine + d15N + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = acm_model)
step(mod)
mod_acm=lmer(z ~ c + ph + richness + funrich + (1 | siteIDD/plotID),data=acm_model)

# creat the effect size plot
set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'Arial',   #To change the font type
          axis.title.size = 1.5,  #To change axis title size
          axis.textsize.x = 1,  #To change x axis text size
          axis.textsize.y = 1,
          title.size = 1.5,
          title.align= "center")  #To change y axis text size

p1=plot_model(mod_acm,axis.labels = c("Fun.rich","Pla.rich","pH"),colors="red",rm.terms = "c",title="ACM (N=87)",axis.lim=c(-1, 1))

## for the model does not include root traits
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio4+ bio8  + bio12 + bio15 + bio18 + spei + richness + funrich  + (1 | siteIDD / plotID), data = acm_model_rich)

step(mod)
mod_acm_rich=lmer(z ~ c + organicCPercent + bio15 + funrich + (1 | siteIDD/plotID),data=acm_model_rich)

p10=plot_model(mod_acm_rich,axis.labels = c("Fun.rich","Pre.seas.","SoilC"),colors=c("blue","red"),rm.terms="c",title="ACM (N=319)",axis.lim=c(-1, 1))


# for the prediction
effects_ph <- effects::effect(term= "ph", mod= mod_acm)
summary(effects_ph) 
x_ph <- as.data.frame(effects_ph)
# creat a plot showing the prediction
ggplot() + 
  geom_point(data=acm_model1, aes(ph, z,color=siteIDD),alpha=0.5) + 
  geom_point(data=x_ph, aes(x=ph, y=fit), color="black") +
  geom_line(data=x_ph, aes(x=ph, y=fit), color="black") +
  geom_ribbon(data= x_ph, aes(x=ph, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  labs(x="pH", y=expression(italic(z)))+
  guides(color=FALSE)







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
