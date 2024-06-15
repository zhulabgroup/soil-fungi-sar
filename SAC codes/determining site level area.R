# these codes determine the species richness at the site level
#note, the data was unrarefied. For rarefied data, please see the "1-..." file
setwd("/Users/luowenqi/Desktop/sar")
neon_dob <- readRDS("/Users/luowenqi/Desktop/sar/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project") == "NEON")
rm(neon_dob)

neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0)
d <- sample_data(neon) #
d <- data.frame(d)
plotID <- substr(d$geneticSampleID, 6, 8)
d <- cbind(d, plotID)
a <- unique(d$Site)

site <- numeric() # look at the number of plots per site
for (i in 1:length(a))
{
  a1 <- subset(d, Site == a[i])
  site[i] <- length(unique(a1$plotID)) # the number of plots in one site
}

pp <- vector("list", 45) # location of the plots within a site
for (i in seq_along(pp))
{
  a1 <- subset(d, Site == a[i])[, c("lon", "lat")]
  a1 <- unique(a1)
  pp[[i]] <- ggplot() +
    geom_point(data = a1, aes(x = lon, y = lat), color = "red") # the number of
}

dtlon <- numeric() # the lon and lat range of the plots within a site
for (i in 1:length(a))
{
  a1 <- subset(d, Site == a[i])[, c("lon", "lat")]
  dtlon[i] <- max(a1$lon) - min(a1$lon)
}

dtlat <- numeric()
for (i in 1:length(a))
{
  a1 <- subset(d, Site == a[i])[, c("lon", "lat")]
  dtlat[i] <- max(a1$lat) - min(a1$lat)
}
site_area <- data.frame(cbind(a, dtlon, dtlat))
site_area$dtlon <- as.numeric(site_area$dtlon)
site_area$dtlat <- as.numeric(site_area$dtlat)

# estimate an area that can cover all the plots in any investigated sites
# 1 degree of lon==54.6 miles
# 1 degree of lat=69 miles
# https://www.usgs.gov/faqs/how-much-distance-does-a-degree-minute-and-second-cover-your-maps#:~:text=One%2Ddegree%20of%20longitude%20equals,one%20second%20equals%2080%20feet.

dit1 <- site_area$dtlon * 54.6
dit2 <- site_area$dtlat * 69
site_area <- cbind(site_area, dit1, dit2) #
# the site for the STEI site covers a quite large area, the maxium sampling area is 27.61029*22.67657(28*23 miles)
