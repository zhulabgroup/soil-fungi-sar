# model selection for each model

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

# site included as the random effect, here some sites have only samples from the same plot
mod <- lmer(z ~ c + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD), data = mode.data1)
step(mod)
#best model
mod=lmer(z ~ c + organicCPercent + nitrogen + sand + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD),data=mode.data1)
r.squaredGLMM(mod)

2. # climate, soil and plant diversity model

mode.data2 <- model_data[, c(1:20)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data2 <- subset(mode.data2, siteIDD != "GUAN" & z < 10 & richness > 0) # some sites do have plant data
# check colinearity among variables
cor(mode.data2[, 4:20])
ggcorrplot(cor(mode.data2[, 4:20]), hc.order = TRUE, type = "lower", lab = TRUE) #
mode.data2[, 4:20] <- apply(mode.data2[, 4:20], 2, range01) %>% data.frame()

# when site and plot are nested,  c was included, only bio12 and soil C was significant, and c and funrich still affect the z values
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + bio12 + bio15 + bio4 + spei + funrich + bio1 + (1 | siteIDD ), data = mode.data2)

step(mod)
mod=lmer(z ~ organicCPercent + c + nitrogen + richness + sand + bio18 + bio12 + bio15 + bio4 + spei + funrich + bio1 + (1 | siteIDD),data=mode.data2)
summary(mod)
r.squaredGLMM(mod)
3.# climate, soil, plant and root model

mode.data3 <- model_data[, c(2:27)] # GUAN don't have soil variables and will be excluded(possibly this site is out of place)
mode.data3 <- subset(mode.data3, siteIDD != "GUAN" & z < 10 & richness > 0 & fine > 0 & rootc > 0) # some sites do have plant data
# check colinearity among variables
cor(mode.data3[, 3:26])
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
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD), data = mode.data3)

step(mod)
mod=lmer(z ~ c + richness + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD),data=mode.data3)
summary(mod)



