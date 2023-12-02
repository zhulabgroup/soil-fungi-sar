# step 8 is about model build but before doing this we need to combine different datasets for modeling
# the resulting data that included z, c and other variables will be used for the construction of guild-based model
# for this data, if we select a plot, there are 30 rows corredponding to the 30 simulated values
#load the required packages
library(ggcorrplot)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggeffects)
library(MuMIn)


plot_loca_all_soil_climate_mean # soil variables# based on step 2 and 7#
plot_plant_rich # plant data# based on step 3#
rich_plot # fungal richness# based on step 4#
core_mass_type # root mass data#based on step 5#
root.chemi.mean # root trait data# based on step 6#

# merging these data sets generates an object of 'd5' that will be used at line 56(or around).
# computing the z and the log(c) value
# this is based on the calculation in step 1#
a <- list()
for (i in 1:dim(all_z_ranall)[1]) # here the initial all_z was saved as "all_z.ranll"
{
  a[[i]] <- t(all_z_ranall[i, 3:32])
}

b <- a[[1]]
for (i in 2:dim(all_z_ranall)[1]) # here the initial "all_z" was saved as "all_z.ranll"
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(all_z_ranall$a1, each = 30)

all_z_30 <- cbind(plotID, b) # for each plot,with 30 estimated z values
all_z_30 <- data.frame(all_z_30)
all_z_30$V2 <- as.numeric(all_z_30$V2)
# for the c values

a <- list()
for (i in 1:dim(all_c_ranall)[1])
{
  a[[i]] <- t(all_c_ranall[i, 3:32])
}

b <- a[[1]]
for (i in 2:dim(all_c_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

plotID <- rep(all_c_ranall$a1, each = 30)
all_c_30 <- cbind(plotID, b) # for each plot,with 30 estimated c values, which can be used to estimated species per unit sample(estimated=log(c))

all_c_30 <- data.frame(all_c_30)
all_c_30$V2 <- as.numeric(all_c_30$V2)
# combind the z and c
com_z_c <- cbind(all_z_30, all_c_30[, 2])
com_z_c <- subset(com_z_c, V2 < 10) ## need to remove the plots with <2 cores
names(com_z_c) <- c("plotID", "z", "log(c)")
com_z_c <- cbind(com_z_c, c = 2.71828^com_z_c$`log(c)`) # the estimated c and z are negatively correlated
model_data <- merge(com_z_c[, c(1, 2, 4)], d5, by = "plotID", all.x = TRUE)#

# adding the site
siteIDD <- substr(model_data$plotID, 1, 4)
model_data <- cbind(siteIDD, model_data)
model_data <- subset(model_data, z < 10) # some sites have explanatory variables but none z,

write.csv(model_data, "model_data.csv")
model_data <- model_data[, -5]
model_data$rich <- as.numeric(model_data$rich)

write.csv(model_data, "model_data.csv")# this is a key data set saved for the downstream analyses.



1. # climate and soil model:both dob and neon sites were included, here only the plotID was treated as a random effect
mode.data1 <- model_data[, c(1:18,20)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data1 <- subset(mode.data1, siteIDD != "GUAN" & z < 10)
# check colinearity among variables

cor(mode.data1[, 4:19])
# mapping the correlation
ggcorrplot(cor(mode.data1[, 4:19]), hc.order = TRUE, type = "lower", lab = TRUE) #
# bold was related with soil C, and MAT;cec was related with soil C; I decided to exclude cec and bold
# build model with plot included as the random effect
# standardized the data.

range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

mode.data1[, 4:19] <- apply(mode.data1[, 4:19], 2, range01) %>% data.frame()

# site and plot are nested
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand +  bio1+bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich  + (1 | siteIDD / plotID), data = mode.data1)

# plotID as the only random effect (483 plots)

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = mode.data1)

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data1)


ggpredict(mod, terms = c("funrich", "siteIDD"), type = "re") %>% 
  plot()



effect_all_clima_soil <- summary(mod)
effect_all_clima_soil <- effect_all_clima_soil$coefficients
effect_all_clima_soil <- data.frame(effect_all_clima_soil)
sig <- effect_all_clima_soil$Pr...t..
sig <- round(sig, 3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
effect_all_clima_soil <- cbind(effect_all_clima_soil, sig)

# make the predictions, with the neon data only

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | site), data = d2)

pred.mm <- ggpredict(mod, terms = c("funrich","site"),type="re")  # this gives overall predictions for the model

ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +         # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = d2, aes(x = funrich, y = z, colour = site),alpha=0.5) + 
  labs(x = "fungal diversity", y = "z", 
       title = "") + 
  theme_minimal()


ggpredict(mod, terms = c("funrich", "site"), type = "re") %>% 
  ggplot() +
  theme_minimal()


modinter <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 | plotID), data = mode.data1)



# the random slope model seems better
ranslope <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (0 + c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 | plotID), data = mode.data1)


# what if c is not included in the model?
mod <- lmer(z ~ organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = mode.data1)
summary(mod)
# site as the only random effect

mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data1)
step(mod)

# best model
mod <- lmer(z ~ c + organicCPercent + nitrogen + sand + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data1)
summary(mod)


## climate, plant diversity and soil model: only neon sites were included
mode.data2 <- model_data[, c(1:20)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data2 <- subset(mode.data2, siteIDD != "GUAN" & z < 10 & richness > 0) # some sites do have plant data
# check colinearity among variables
cor(mode.data2[, 4:20])
ggcorrplot(cor(mode.data2[, 4:20]), hc.order = TRUE, type = "lower", lab = TRUE) #
mode.data2[, 4:20] <- apply(mode.data2[, 4:20], 2, range01) %>% data.frame()

# when site and plot are nested,  c was included, only bio12 and soil C was significant, and c and funrich still affect the z values
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + bio12 + bio15 + bio4 + spei + funrich + bio1 + (1 | siteIDD / plotID), data = mode.data2)
# extract the fixed effect
effect_neon_clima_soil <- summary(mod)
effect_neon_clima_soil <- effect_neon_clima_soil$coefficients
effect_neon_clima_soil <- data.frame(effect_neon_clima_soil)
sig <- effect_neon_clima_soil$Pr...t..
sig <- round(sig, 3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
effect_neon_clima_soil <- cbind(effect_neon_clima_soil, sig)

# if we just include the plotID as the random effect

mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + bio12 + bio15 + bio4 + spei + funrich + bio1 + (1 | plotID), data = mode.data2)

# site as the random effect
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + bio12 + bio15 + bio4 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data2)

# extract the fixed effect
effect_neon_clima_soil_plotran <- summary(mod)
effect_neon_clima_soil_plotran <- effect_neon_clima_soil_plotran$coefficients
effect_neon_clima_soil_plotran <- data.frame(effect_neon_clima_soil_plotran)
sig <- effect_neon_clima_soil_plotran$Pr...t..
sig <- round(sig, 3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
effect_neon_clima_soil_plotran <- cbind(effect_neon_clima_soil_plotran, sig)

pred.mm <- ggpredict(mod, terms = c("richness"))  # this gives overall predictions for the model

#population level

d3$predm2=predict(mod,re.form=NA)
d3$predm2fit=predict(mod)

ggplot(d3,aes(richness,z))+
  geom_point()+
  geom_line(colour="red",aes(y=predm2))+
  geom_line(colour="dark grey",aes(y=predm2fit,group=siteIDD)) 

# the relationship between the tree diversity and 

ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted,colour=group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = mode.data2, aes(x = richness, y = z, colour = siteIDD),alpha=0.1) + 
  labs(x = "fungal diversity", y = "z", 
       title = "Body length does not affect intelligence in dragons") + 
  theme_minimal()+guides(colour=FALSE)




# random slope model for the neon data?

model_neon_slope <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (0 + organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 | siteIDD / plotID), data = mode.data2)

# adding the code, substr the plotID generates the same results
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD / code), data = mode.data2)

# if we included only the plotID as the random effect, soil N(-),sand (-),bio8(-),bio4(+),bio12(+),bio15(+),bio1(-)
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = mode.data2)
# if we included only the site as the random effect, soil N(-),rich (+),sand (-),bio8(+),bio18(-),bio4(+),bio12(+),funrich(+)

mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data2)

# best model
mod <- lmer(z ~ organicCPercent + c + sand + bio12 + bio15 + funrich + bio1 + (1 | siteIDD / plotID), data = mode.data2)
## only site as random
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data2)

# only plot as random
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = mode.data2)

# climate, soil, plant and root model# (104 plots for the NEON sites)

mode.data3 <- model_data[, c(2:28)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data3 <- subset(mode.data3, siteIDD != "GUAN" & z < 10 & rich > 0 & fine > 0 & rootc > 0) # some sites do have plant data
# check colinearity among variables
cor(mode.data3[, 3:27])
# coarse and fine roots are correlated;bio1 and bio4, bio1 and bold; organic and bold; organic and cec; root n and root cn;
# i therefor excluded coarse, bio4,bold,cec and root n
ggcorrplot(cor(mode.data3[, 3:27]), hc.order = TRUE, type = "lower", lab = TRUE) #
mode.data3[, 3:27] <- apply(mode.data3[, 3:27], 2, range01) %>% data.frame()
# bold-bio1 related with r>0/75
# bio1-bio4(was excluded)
# bold-soil oc
# coarse-fine
# cec-bold
# root n-root cn

# site and plot are nested, only funrich showed a significant effect with plant rich showing a marginal effect
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = mode.data3)

effect_neon_clima_soil_root <- summary(mod)
effect_neon_clima_soil_root <- effect_neon_clima_soil_root$coefficients
effect_neon_clima_soil_root <- data.frame(effect_neon_clima_soil_root)
sig <- effect_neon_clima_soil_root$Pr...t..
sig <- round(sig, 3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
effect_neon_clima_soil_root <- cbind(effect_neon_clima_soil_root, sig)

# when only plotID was included as the random effect model
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | plotID), data = mode.data3)

effect_neon_clima_soil_root_plotran <- summary(mod)
effect_neon_clima_soil_root_plotran <- effect_neon_clima_soil_root_plotran$coefficients
effect_neon_clima_soil_root_plotran <- data.frame(effect_neon_clima_soil_root_plotran)
sig <- effect_neon_clima_soil_root_plotran$Pr...t..
sig <- round(sig, 3)
sig[sig > 0.05] <- "no"
sig[sig < 0.05] <- "sig"
effect_neon_clima_soil_root_plotran <- cbind(effect_neon_clima_soil_root_plotran, sig)



# random slope model

rand_slop_root <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = mode.data3)


# site as the only random effect
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD), data = mode.data3)
step(mod)
# best model
mod <- lmer(z ~ c + rich + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD), data = mode.data3)

# plot as the random
# site as the only random effect
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + rich + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | plotID), data = mode.data3)
step(mod)
# best model
mod <- lmer(z ~ c + rich + sand + bio8 + bio15 + funrich + rootc + (1 | plotID), data = mode.data3)
#
#
## for different guilds
####
