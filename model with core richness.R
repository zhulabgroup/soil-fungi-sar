# when core level richness was considered while c was included in the model while
library(lmerTest)
library(lme4)
library(MuMIn)
1.# determine core level fungal richness

# test the relationship between the number of soil cores and the estimated z value
load("rare_all.Rdata") # load the all rarefied data
d <- sample_data(rare_all)
table(d$Project) # the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1 <- data.frame(d[, c("geneticSampleID", "Site")])
plotID <- substr(d1$geneticSampleID, 1, 8)
d1 <- cbind(d1, plotID)
iddob <- d1$Site[1:908] # a site corresponds to a plot
idneon <- d1$plotID[909:6378] # an unique plotID corresponds to a plot
plotIDM <- data.frame(c(iddob, idneon))
names(plotIDM) <- "plotIDM" # the plot id used for the SAR
row.names(plotIDM) <- row.names(d)
plotIDM <- sample_data(plotIDM)
d <- merge_phyloseq(rare_all, plotIDM) # merge the new plotid with the initial data
# select an unique plot and build a SAR within the plot
a1 <- sample_data(d) # the unique plotID, we have 476 plots
a1 <- unique(a1$plotIDM)
rare_all <- d

## the mean richness per core for a site(alpha diversity)

rich<- numeric()

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  d1=estimate_richness(data_sub,measures="Observed")
  rich[i]=mean(d1$Observed)
}

core_rich=cbind(a1,rich)
core_rich=data.frame(core_rich)
names(core_rich)=c("plotID","corich")

# save the data
save(core_rich,file="core_rich.RData")

data_corich=merge(model_data[,1:20],core_rich,by="plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10)
data_corich[, 5:21] <- apply(data_corich[, 5:21], 2, range01) %>% data.frame()

2.# climate and soil model

mod <- lmer(z ~ corich+ organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | siteIDD / plotID), data = data_corich)
mod_corich=lmer(z ~ corich + bio4 + bio12 + bio15 + (1 | siteIDD/plotID),data=data_corich)

3.# climate, soil and plant richness model
data_corich=merge(model_data[,1:20],core_rich,by="plotID")
data_corich <- subset(data_corich, siteIDD != "GUAN" & z < 10&richness>0)
data_corich[, 5:21] <- apply(data_corich[, 5:21], 2, range01) %>% data.frame()
mod <- lmer(z ~ corich+ organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei +richness+ funrich + bio1 + (1 | siteIDD / plotID), data = data_corich)


4.# # climate, soil,plant richness and root traits model 

data_corich=merge(model_data[,1:27],core_rich,by="plotID")
data_corich=subset(data_corich,siteIDD != "GUAN" & z < 10 & richness > 0 & fine > 0 & rootc > 0)
data_corich[, 5:27] <- apply(data_corich[, 5:27], 2, range01) %>% data.frame()
mod <- lmer(z ~ organicCPercent + c + ph + nitrogen + richness + sand + bio2 + bio8 + bio18 + +bio12 + bio15 + spei + funrich + bio1 + fine + d13C + rootc + rootcn + (1 | siteIDD / plotID), data = data_corich)
step(mod)
mod <- lmer(z ~ organicCPercent + corich + richness + funrich + rootc + (1 | siteIDD/plotID), data = data_corich)
