setwd("/../../../sar")
library(doParallel)
library(raster)
library(sp)
install.packages("devtools")
require(devtools)
devtools::install_github("seschaub/getSpei")
require(getSpei)
require(ncdf4)
require(chron)
# read in the data
neon_dob <- readRDS("/Users/luowenqi/Desktop/sar/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project") == "NEON")
rm(neon_dob)
# neon <- subset_samples(neon, horizon == "O")
neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0) # both soil horizons were included
# rm(neon_dob)

# neon dataset
d <- sample_data(neon) # sample data data frame
d <- data.frame(d)
plotID <- substr(d$geneticSampleID, 1, 8) # get the id for each plot so that we can estimate the SAR based on the plot
d <- sample_data(neon)
plotID <- data.frame(plotID)
row.names(plotID) <- row.names(d)
plotID <- sample_data(plotID)

d <- merge_phyloseq(neon, plotID) # adding the plotID to the sample data

# the coordinates of the plot
plot_loca <- sample_data(d)[, c("lon", "lat")]
# extract the climate variables for each NEON plot
head(plot_loca)
plot_loca <- data.frame(plot_loca)
plot_loca <- cbind(plotID, plot_loca)
plot_loca <- unique(plot_loca)
## get the 18 climate variables
r <- getData("worldclim", var = "bio", res = 10)
points <- SpatialPoints(plot_loca[, 2:3], proj4string = r@crs)
values <- extract(r, points)
df <- cbind.data.frame(coordinates(points), values)
# select the seven climate variables that do not show colinearity
##
r <- getData("worldclim", var = "bio", res = 10)
points <- SpatialPoints(plot_loca[, 2:3], proj4string = r@crs)
values <- extract(r, points)
df <- cbind.data.frame(coordinates(points), values)
# select the seven climate variables
neon_climate <- df[, c("bio2", "bio8", "bio18", "bio4", "bio12", "bio15", "bio1")]
# get the land aridity index
names(plot_loca) <- c("site", "longitude", "latitude")
site_spei <- spec_spei(spei_files = c("spei01"), start_y = 2010, end_y = 2018, locations = plot_loca)
aridity.mean <- aggregate(spei01 ~ location_id, FUN = mean, data = site_spei)
neon_climate <- cbind(neon_climate, plot_loca)
names(aridity.mean) <- c("plotID", "aridity")
# combine the climate and aridity index at the plot level
neon_climate_aridity <- merge(neon_climate, aridity.mean, by = "plotID")
write.csv(neon_climate_aridity, "neon_climate_aridity.csv")
