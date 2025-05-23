library(vegan)
library(geosphere)
neon_dob <- readRDS("../sar-project/data/phylo_V3.1.RDS")
dob <- subset_samples(neon_dob, Project == "DoB")
neon <- subset_samples(neon_dob, Project == "NEON")
rm(neon_dob)

# 1.dob 
dob <- subset_samples(dob, horizon %in% c("O", "M"))
d <- sample_data(dob)

# extract the site and core info from sample name
site_names <- unique(d$Site)
sample_names <- matrix(unlist(strsplit(d$geneticSampleID, "[.]")), ncol = 3, byrow = T)
colnames(sample_names) <- c("Site1", "Location", "Horizon")
d <- cbind(d, sample_names)
sample_data(dob) <- d

# see the distribution of otu in dob
dob <- subset_taxa(dob, taxa_sums(dob) > 0)
dob <- subset_taxa(dob, sample_sums(dob) > 0)
hist((sample_sums(otu_table(dob))))
quantile((sample_sums(otu_table(dob))), 0.1)
set.seed(1000)

# rarefy
dob_rare <- rarefy_even_depth(dob, 3000)# as we are treating the two data sets as a whole, i think we should rarefy the full data set rather than separately?
beta_diversity_mat <- vegdist(otu_table(dob_rare), "jaccard")
d <- sample_data(dob_rare)

# calculate distance
dist_mat <- distm(d[, 1:2], fun = distHaversine)
dist_mat <- as.dist(dist_mat)

plot(dist_mat, beta_diversity_mat)
plot(dist_mat[which(dist_mat <= 41 & dist_mat > 0)],
  beta_diversity_mat[which(dist_mat <= 41 & dist_mat > 0)],
  xlab = "distance", ylab = "Jaccard Distance", main = "DoB"
)
plot(dist_mat[which(dist_mat > 41)],
  beta_diversity_mat[which(dist_mat > 41)],
  xlab = "distance", ylab = "Jaccard Distance", main = "DoB"
)
plot(log(dist_mat),
  beta_diversity_mat,
  xlab = "log(distance)", ylab = "Jaccard Distance", main = "DoB"
)
rownames(dist_mat) <- d$geneticSampleID
colnames(dist_mat) <- d$geneticSampleID

library(RColorBrewer)
mypalatte <- brewer.pal(7, "Set1")
mp <- ggplot() +
  mapworld +
  ylim(25, 90) +
  xlim(-150, -50)

# calculate the distance between dob sites 
# with respect to sample 000 in each site
head_point <- d[which(d$Location == "000"), c("lon", "lat", "Site")]
head_point <- head_point[duplicated(head_point), ]
dist_mat <- distm(head_point[, 1:2], fun = distHaversine)
row.names(dist_mat) <- head_point$Site
colnames(dist_mat) <- head_point$Site
hist(dist_mat[which(dist_mat > 0)],
  breaks = 1000, xlab = "distance between sites",
  main = "histogram of the distance between sites"
)

# plot the distance-decay pattern within each site/plot
for (site_i in site_names) {
  dob_sub <- subset_samples(dob, Site == site_i)
  d_sub <- sample_data(dob_sub)
  dist_mat_sub <- as.dist(distm(d_sub[, 1:2], fun = distHaversine))
  beta_mat_sub <- vegdist(otu_table(dob_sub), method = "jaccard")
  summary(lm(beta_mat_sub ~ dist_mat_sub))
  # png(paste("Fig/", site_i, ".png", sep = ""))
  plot(dist_mat_sub, beta_mat_sub,
    xlab = "distance", ylab = "jaccard distance",
    main = site_i
  )
  # dev.off()
}

# 2. NEON
beta_neon <- readRDS("E:/beta_dist.RDS") # the calculated beta diversity matrix between each 2 samples in all NEON dataset
d <- sample_data(neon)
d <- d[which(row.names(d) %in% attributes(beta_neon)[["Labels"]]),]
dist_neon <- distm(d[, 1:2], fun = distHaversine)
dist_neon <- as.dist(dist_neon)
plot(dist_neon[which(dist_neon > 0)], beta_neon[which(dist_neon > 0)], xlab = "distance", ylab = "jaccard distance")

# distance-decay pattern within plot
neon <- subset_samples(neon, sample_sums(neon) >= 3000)
neon <- subset_samples(neon, taxa_sums(neon) > 0)
d <- sample_data(neon)
plot_names <- unique(substr(d$geneticSampleID, 1, 8))
d$plot <- substr(d$geneticSampleID, 1, 8)
sample_data(neon) <- d

for (site_i in plot_names) {
  neon_sub <- subset_samples(neon, plot == site_i)
  neon_sub <- rarefy_even_depth(neon_sub, rngseed = 10, sample.size = 3000, replace = F)
  d_sub <- sample_data(neon_sub)
  
  if (nrow(d_sub) <= 3) next
  
  loc_within_plot <- matrix(as.numeric(unlist(strsplit(gsub(".*[OM]-(\\d+\\.?\\d*-\\d+\\.?\\d*)-.*", "\\1", d_sub$geneticSampleID), "-"))), ncol = 2, byrow = T)
  row.names(loc_within_plot) <- row.names(d_sub)
  dist_within_plot <- dist(loc_within_plot)
  beta_within_plot <- vegdist(otu_table(neon_sub), method = "jaccard")
  
  #png(paste("Fig/", site_i, ".png", sep = ""))
  plot(dist_within_plot, beta_within_plot,
    xlab = "distance", ylab = "jaccard distance",
    main = site_i
  )
  #dev.off()
}

d$loc_within_plot <- matrix(as.numeric(unlist(strsplit(gsub(".*[OM]-(\\d+\\.?\\d*-\\d+\\.?\\d*)-.*", "\\1", d$geneticSampleID), "-"))), ncol = 2, byrow = T)
d$coord_x <- NA
d$coord_y <- NA
for(plot_i in plot_names){  
  d_sub <- d[which(d$plot == plot_i),]
  site_coord <- unique(d_sub[, 1:2])  # For instance, New York City coordinates
  
  # Relative distances in x and y axes from the known plot
  relative_distance_x <- d_sub$loc_within_plot[,1]  # Hypothetical relative distance on x-axis
  relative_distance_y <- d_sub$loc_within_plot[,2]  # Hypothetical relative distance on y-axis
  
  # Calculate direct distance using Pythagorean theorem
  direct_distance <- sqrt(relative_distance_x^2 + relative_distance_y^2)
  
  # Calculate angle using inverse tangent (arctan) function
  angle <- atan2(relative_distance_y, relative_distance_x) * (180 / pi)  # Convert radians to degrees
  
  # Calculate the new plot's coordinates based on the known plot
  coord <- t(apply(cbind(angle, direct_distance), 1, function(x) destPoint(p = site_coord, b = x[1], d = x[2])))
  d$coord_x[which(d$plot == plot_i)] <- coord[,1]
  d$coord_y[which(d$plot == plot_i)] <- coord[,2]
}