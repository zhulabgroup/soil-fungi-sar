# model selection for each model
# for the prestented figures, site and plot are nested in the random effect setting

library(sjPlot)
library(ggplot2)
library(effects)
library(lme4)
library(lmerTest)

# function for data standardization
range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

# plot the results of the fixed effects
set_theme(
  base = theme_classic(), # To remove the background color and the grids
  theme.font = "Arial", # To change the font type
  axis.title.size = 2.0, # To change axis title size
  axis.textsize.x = 1, # To change x axis text size
  axis.textsize.y = 1,
  title.size = 2,
  title.align = "center"
) # To change y axis text size

0. # load the [model_data] from data folder

1. # climate and soil model:both dob and neon sites were included, here only the plotID was treated as a random effect
mode.data1 <- model_data[, c(1:18, 20)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data1 <- subset(mode.data1, siteIDD != "GUAN" & z < 10)
# check colinearity among variables

cor(mode.data1[, 4:19])
# mapping the correlation
ggcorrplot(cor(mode.data1[, 4:19]), hc.order = TRUE, type = "lower", lab = TRUE) #
# bold was related with soil C, and MAT;cec was related with soil C; I decided to exclude cec and bold
# build model with plot included as the random effect

mode.data1[, 4:19] <- apply(mode.data1[, 4:19], 2, range01) %>% data.frame()

# plot and site included as nested random effects, here some sites have only samples from the same plot
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD / plotID), data = mode.data1)
step(mod)

# best model
mod_nest <- lmer(z ~ c + nitrogen + sand + bio2 + bio12 + bio15 + funrich + (1 | siteIDD / plotID), data = mode.data1)
p1 <- plot_model(mod_nest, axis.labels = c("Fun.rich", "Pre.seas.", "MAP", "MDR", "Sand", "SoilN", "c"), color = c("blue", "red"), rm.terms = "c", title = "Climate+Soil \n(N=483)", dot.size = 4)

2. # climate, soil and plant diversity model

mode.data2 <- model_data[, c(1:20)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data2 <- subset(mode.data2, siteIDD != "GUAN" & z < 10 & richness > 0) # some sites do have plant data
# check colinearity among variables
cor(mode.data2[, 4:20])
ggcorrplot(cor(mode.data2[, 4:20]), hc.order = TRUE, type = "lower", lab = TRUE) #
mode.data2[, 4:20] <- apply(mode.data2[, 4:20], 2, range01) %>% data.frame()

#
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + bio12 + bio15 + bio4 + spei + funrich + bio1 + (1 | siteIDD / plotID), data = mode.data2)

step(mod)

mod_nestCSP <- lmer(z ~ c + nitrogen + richness + sand + bio2 + bio12 + bio15 + funrich + (1 | siteIDD / plotID), data = mode.data2)

summary(mod_nestCSP)

p2 <- plot_model(mod_nestCSP, axis.labels = c("Fun.rich", "Pre.seas.", "MAP", "MDR", "Sand", "Pla.rich", "SoilN", "c"), colors = c("blue", "red"), rm.terms = "c", title = "Climate+Soil+Plant \n(N=438)", dot.size = 4)


3. # climate, soil, plant and root model

mode.data3 <- model_data[, c(1:28)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data3 <- subset(mode.data3, siteIDD != "GUAN" & z < 10 & richness > 0 & fine > 0 & rootc > 0) # some sites do have plant data
# check colinearity among variables
cor(mode.data3[, 4:27])
# coarse and fine roots are correlated;bio1 and bio4, bio1 and bold; organic and bold; organic and cec; root n and root cn;
# i therefor excluded coarse, bio4,bold,cec and root n
ggcorrplot(cor(mode.data3[, 3:26]), hc.order = TRUE, type = "lower", lab = TRUE) #
mode.data3[, 3:26] <- apply(mode.data3[, 3:26], 2, range01) %>% data.frame()
# bold-bio1 related with r>0/75
# bio1-bio4(was excluded)
# bold-soil oc
# coarse-fine
# cec-bold
# root n-root cn

# site and plot are nested, only funrich showed a significant effect with plant rich showing a marginal effect
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = mode.data3)
step(mod)
mod_nestCSPR <- lmer(z ~ c + richness + funrich + bio1 + (1 | siteIDD / plotID), data = mode.data3)
# create the effect size plot
p3 <- plot_model(mod_nestCSPR, axis.labels = c("MAT", "Fun.rich", "Pla.rich"), color = c("blue", "red"), rm.terms = "c", title = "Climate+Soil+Plant+Root \n(N=104)", dot.size = 5)

plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO")
