library(vegan)
library(phyloseq)

0. # an example showing how to determine the coordinates of a new point based on a reference point, the distance and bearing angle between them with Haversine formula
calculate_coordinates <- function(lat1, lon1, distance_km, bearing_deg) {
  # Radius of the Earth in kilometers
  earth_radius_km <- 6371
  # Convert latitude and longitude from degrees to radians
  lat1 <- deg2rad(lat1)
  lon1 <- deg2rad(lon1)
  bearing_rad <- deg2rad(bearing_deg)
  # Calculate the new latitude
  lat2 <- asin(sin(lat1) * cos(distance_km / earth_radius_km) + cos(lat1) * sin(distance_km / earth_radius_km) * cos(bearing_rad))
  # Calculate the new longitude
  lon2 <- lon1 + atan2(sin(bearing_rad) * sin(distance_km / earth_radius_km) * cos(lat1), cos(distance_km / earth_radius_km) - sin(lat1) * sin(lat2))

  # Convert latitude and longitude back to degrees
  lat2 <- rad2deg(lat2)
  lon2 <- rad2deg(lon2)

  return(c(lat2, lon2))
}

# Example usage
lat1 <- 40.7128 # Latitude of the starting point (e.g., New York City)
lon1 <- -74.0060 # Longitude of the starting point
distance_km <- c(100, 200, 112) # Distance in kilometers
bearing_deg <- c(45, 36, 15) # Bearing in degrees

# Calculate the coordinates of the new point
new_coordinates <- calculate_coordinates(lat1, lon1, distance_km, bearing_deg)

cat("New Latitude:", new_coordinates[1], "\n")
cat("New Longitude:", new_coordinates[2], "\n")



1. # adding the gx and gy for each core within the 40 x 40 m plot for the neon sites
load("rare_all.Rdata") # load the all rarefied data
d <- sample_data(rare_all)
table(d$Project) # the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1 <- data.frame(d[, c("geneticSampleID", "Site")])
plotID <- substr(d1$geneticSampleID, 1, 8)
d1 <- cbind(d1, plotID)
iddob <- d1$Site[1:908] # a site correspondes to a plot
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

# subset the neon data
sub_neon <- subset_samples(rare_all, get_variable(rare_all, "Project") == "NEON")
sub_neon <- subset_samples(sub_neon, get_variable(sub_neon, "horizon") != "AH")
sub_neon <- subset_samples(sub_neon, get_variable(sub_neon, "horizon") != "OH")

# the location of each core (gx and gy) within the 40 by 40 m plot
loca <- data.frame(sample_data(sub_neon))["geneticSampleID"] # need to be
loca <- substr(loca$geneticSampleID, 11, 20)
loca <- data.frame(loca)
d <- strsplit(loca$loca, "-")
a <- matrix(nrow = dim(loca)[1], ncol = 2) # the location of each core within the 40 by 40 plot
for (i in 1:dim(loca)[1])
{
  a[i, ] <- as.numeric(d[[i]][2:3])
}
# need to cbind the location data with the plotIDM
a <- data.frame(a)
names(a) <- c("gx", "gy")
row.names(a) <- row.names(sample_data(sub_neon))
a <- sample_data(a)
sub_neon <- merge_phyloseq(sub_neon, a) # adding the location data to the full data.
a1 <- sample_data(sub_neon) # the unique plotID, we have 476 plots
a1 <- unique(a1$plotIDM) # 4


2. # get the bearing angle for each soil core, relative to the origin. let's do this based on the plotID.

bear0 <- list()
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  angle_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples for each plot
  k <- sample_data(angle_sub)[, 25] / sample_data(angle_sub)[, 26]
  bear0[[i]] <- data.frame(atan(k) * 180 / pi)
}
# the bearing angles for all cores based on plotID
bear1 <- bear0[[1]]
for (i in 2:length(a1)) {
  bear1 <- rbind(bear1, bear0[[i]])
}

bear1$gx <- gsub("NaN", "90", bear1$gx)
bear1$gx <- as.numeric(bear1$gx)
# the distance from each core to the origin


3. # compute the distance of each core to the origin, let's do this based on the plotID

pair <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples for a given plotID
  location <- sample_data(data_sub)[, c("gx", "gy")]
  dis_tance <- sqrt((location$gx)^2 + (location$gy)^2)
  pair[[i]] <- matrix(dis_tance)
}
# bind all the distance output
pair_cb <- pair[[1]]
for (i in 2:length(a1))
{
  pair_cb <- rbind(pair_cb, pair[[i]])
}

pair_cb <- data.frame(pair_cb)

# add the distance and bearing information to the full dataset

bear_dis <- cbind(bear1, pair_cb)
sam_names <- data.frame(rownames(sample_data(sub_neon)))
bear_dis <- bear_dis[sam_names$rownames.sample_data.sub_neon.., ] # reorder the rows of the data so that it can binded with the full data
names(bear_dis) <- c("angle", "distance")

bear_dis <- sample_data(bear_dis)

sub_neon <- merge_phyloseq(sub_neon, bear_dis)

4. # for a given plotID, determined the coordinates of each soil core

cod <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples in a given plotID
  location <- sample_data(data_sub)[, c("lon", "lat", "distance", "angle")]
  location <- data.frame(location)
  new_coordinates <- calculate_coordinates(location[1, 1], location[1, 2], location[, 3] / 1000, location[, 4])
  cod[[i]] <- matrix(new_coordinates, nrow = dim(location)[1], ncol = 2, byrow = FALSE)
  rownames(cod[[i]]) <- rownames(location)
}
# bind all the plots
cod_cb <- cod[[1]]
for (i in 2:length(a1))
{
  cod_cb <- rbind(cod_cb, cod[[i]])
}
cod_cb <- data.frame(cod_cb)

cod_cb <- cod_cb[sam_names$rownames.sample_data.sub_neon.., ] # reorder the rows of the data to make it consistent with the full data
names(cod_cb) <- c("flon", "flat") # the final lon and lat

cod_cb <- sample_data(cod_cb)

neon_cod <- merge_phyloseq(sub_neon, cod_cb) # the neon dataset with each core having an unique coordiate
save(neon_cod, file = "neon_cod.Rdata")

# look at the diversity and distance

beta_diversity_neon <- vegdist(otu_table(neon_cod), "jaccard") # beta diversity between all pairwise soil cores
d <- sample_data(dob)
dist_mat <- distm(d[, 1:2], fun = distHaversine)
dist_mat <- as.dist(dist_mat)
dist_mat <- matrix(dist_mat)
beta_diversity_mat <- matrix(beta_diversity_mat)

dob_dis <- cbind(beta_diversity_mat, dist_mat)
dob_dis <- data.frame(dob_dis)
names(dob_dis) <- c("beta", "distance")
