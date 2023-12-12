# model fitting with the npp data
library(ggcorrplot)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggeffects)
library(MuMIn)

# function

range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

# model include climate, soil, plant, fungal richness and npp
mode.data2 <- model_data[, c(1:20,28)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)

mode.data2 <- subset(mode.data2, siteIDD != "GUAN" & z < 10 & richness > 0&npp>0) # some sites do have plant data

# check colinearity among variables
cor(mode.data2[, 4:21])

ggcorrplot(cor(mode.data2[, 4:21]), hc.order = TRUE, type = "lower", lab = TRUE) #
#bold was related with organic c and was excluded

mode.data2[, c(4:6,8:21)] <- apply(mode.data2[, c(4:6,8:21)], 2, range01) %>% data.frame()

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

mod <- lmer(z ~ organicCPercent + c +ph+ nitrogen +sand+ cec+richness +npp +  bio1+bio2 + bio8 + bio18 + bio12 + bio15 + bio4 + spei + funrich + (1 | siteIDD/plotID), data = mode.data2)

step(mod)

mod=lmer(z ~ c + richness + ph+npp + bio12 + bio4 + funrich + (1 | siteIDD/plotID),data=mode.data2)


