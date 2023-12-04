# look at the relationship between the in-situ soil variables and the rastered ones.
# this was tested for three soil variables: soil ph, soil N and soil organic content.

library(tidyverse)
library(phyloseq)

# the soil variables used in the model
# get the mean vlaues of each site
##
1. # read the data
neon_dob <- readRDS("phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)

2. # get the sample ID
d <- sample_data(neon_dob)
colnames(d)
plotid <- substr(d$geneticSampleID, 1, 8)
project <- d$Project
plotid <- cbind(plotid, project) %>% data.frame()

# note: the dob and neon data are not rbind in order so we cant subset different projects based on the row number

# add the plot ID to the data
plotID <- vector()
for (i in 1:dim(plotid)[1])
{
  d2 <- plotid[i, ]

  if (d2[, 2] == "DoB") {
    plotID[i] <- substr(d2$plotid, 1, 3)
  } else {
    plotID[i] <- substr(d2$plotid, 1, 8)
  }
}

# alternatively, the below chunk of codes work

plotID <- vector()
for (i in 1:dim(plotid)[1])
{
  d2 <- plotid[i, ]

  if (str_detect(d2[, 2], "DoB")) # this funcitona can be used to detect specific string
    {
      plotID[i] <- substr(d2$plotid, 1, 3)
    } else {
    plotID[i] <- substr(d2$plotid, 1, 8)
  }
}

# add the plotid to the initial sample data

plotID <- data.frame(plotID)
rownames(plotID) <- rownames(d)
plotID <- sample_data(plotID)
d <- merge_phyloseq(d, plotID)

# get the mean values of the soil variables, ignoring the differences among soil layers
sam <- sample_data(d)
sam <- sam[, c(5, 6:10, 24)]
sam <- data.frame(sam)


# look at the relationship betwwen the rastered and the in-situ measured variables for (N)
sam_n <- subset(sam, nitrogenPercent > 0)
sam_n <- aggregate(nitrogenPercent ~ plotID, data = sam_n, FUN = mean)
#
model_soil <- aggregate(model_data[, c(6:11)], by = list(model_data$plotID), FUN = mean)

names(model_soil)[1] <- "plotID"

# merge the two datasets

soiln <- merge(sam_n, model_soil, by = "plotID") # (r2=0.18, not well related,n=223)

# # look at the relationship betwwen the rastered and the in-situ measured variables for (soil organic carbon,n=235)
sam_c <- subset(sam, organicCPercent > 0)
sam_c <- aggregate(organicCPercent ~ plotID, data = sam_c, FUN = mean)

soilc <- merge(sam_c, model_soil, by = "plotID") # (r2=0.48, moderately related)

# look at the relationship betwwen the rastered and the in-situ measured variables for (soil pH,n=501)
sam_ph <- subset(sam, soilInCaClpH > 0)
sam_ph <- aggregate(soilInCaClpH ~ plotID, data = sam_ph, FUN = mean)

soilph <- merge(sam_ph, model_soil, by = "plotID") # (r2=0.76, highly correlated,n=235)
