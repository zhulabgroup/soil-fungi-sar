library(vegan)
library(geosphere)
neon_dob <- readRDS("../sar-project/data/phylo_V3.1.RDS")
dob <- subset_samples(neon_dob, Project == "DoB")
neon <- subset_samples(neon_dob, Project == "NEON")
rm(neon_dob)
dob <- subset_samples(dob, horizon %in% c("O", "M"))
d <- sample_data(dob)

site_names <- unique(d$Site)
sample_names <- matrix(unlist(strsplit(d$geneticSampleID, "[.]")), ncol = 3, byrow = T)
colnames(sample_names) <- c("Site1", "Location", "Horizon")
d <- cbind(d, sample_names)
sample_data(dob) <- d

dob <- subset_taxa(dob, taxa_sums(dob) > 0)
dob <- subset_taxa(dob, sample_sums(dob) > 0)
hist((sample_sums(otu_table(dob))))
quantile((sample_sums(otu_table(dob))), 0.1)
set.seed(1000)
dob_rare <- rarefy_even_depth(dob, 3000)
beta_diversity_mat <- vegdist(otu_table(dob_rare), "jaccard")
d <- sample_data(dob_rare)
dist_mat <- distm(d[, 1:2], fun = distHaversine)
dist_mat <- as.dist(dist_mat)
plot(dist_mat, beta_diversity_mat)
plot(dist_mat[which(dist_mat <= 41 & dist_mat > 0)], 
     beta_diversity_mat[which(dist_mat <= 41 & dist_mat > 0)],
     xlab = "distance", ylab = "Jaccard Distance", main = "DoB")
plot(dist_mat[which(dist_mat > 41)], 
     beta_diversity_mat[which(dist_mat > 41)],
     xlab = "distance", ylab = "Jaccard Distance", main = "DoB")
plot(log(dist_mat), 
     beta_diversity_mat,
     xlab = "log(distance)", ylab = "Jaccard Distance", main = "DoB")
rownames(dist_mat) <- d$geneticSampleID
colnames(dist_mat) <- d$geneticSampleID
library(RColorBrewer)
# dist_clusters <- kmeans(d[,1:2],centers = 7)[["cluster"]]
# d$cluster <- dist_clusters
mypalatte <- brewer.pal(7, "Set1")
mp <- ggplot()+
  mapworld+
  ylim(25, 90) +
  xlim(-150, -50)
mp2 <- mp +
  geom_point(aes(x = d$lon, y = d$lat)) +
  xlab("longitude") + 
  ylab("latitude")
mp2

head_point <- d[which(d$Location == "000"),c("lon", "lat", "Site")]
head_point <- head_point[duplicated(head_point),]
dist_mat <- distm(head_point[, 1:2], fun = distHaversine)
row.names(dist_mat) <- head_point$Site
colnames(dist_mat) <- head_point$Site
hist(dist_mat[which(dist_mat > 0)], breaks = 1000, xlab = "distance between sites",
     main = "histogram of the distance between sites")

for (site_i in site_names){
  dob_sub <- subset_samples(dob, Site == site_i)
  d_sub <- sample_data(dob_sub)
  dist_mat_sub <- as.dist(distm(d_sub[, 1:2], fun = distHaversine))
  beta_mat_sub <- vegdist(otu_table(dob_sub), method = "jaccard")
  summary(lm(beta_mat_sub ~ dist_mat_sub))
  
  png(paste("Fig/",site_i,".png", sep = ""))
  plot(dist_mat_sub, beta_mat_sub, xlab = "distance", ylab = "jaccard distance",
       main = site_i)
  dev.off()
}

beta_neon <- readRDS("E:/beta_dist.RDS")
d <- read.table("E:/sample_data_neon_dob.txt")
d <- d[which(d$Project == "NEON"),]
d <- d[which(rownames(d) %in% attributes(beta_neon)[["Labels"]]),]
dist_neon <- distm(d[, 1:2], fun = distHaversine)
dist_neon <- as.dist(dist_neon)
plot(dist_neon[which(dist_neon > 0)], beta_neon[which(dist_neon > 0)], xlab = "distance", ylab = "jaccard distance")

loc_within_plot <- matrix(as.numeric(unlist(strsplit(gsub(".*[OM]-(\\d+\\.?\\d*-\\d+\\.?\\d*)-.*", "\\1", d$geneticSampleID), "-"))), ncol = 2, byrow = T)
row.names(loc_within_plot) <- row.names(d)
neon <- subset_samples(neon, sample_sums(neon) >= 3000)
neon <- subset_samples(neon, taxa_sums(neon) > 0)

d <- sample_data(neon)
plot_names <- unique(substr(d$geneticSampleID, 1, 8))
d$plot <- substr(d$geneticSampleID, 1, 8)
sample_data(neon) <- d
for (site_i in plot_names){
  neon_sub <- subset_samples(neon, plot == site_i)
  neon_sub <- rarefy_even_depth(neon_sub, rngseed = 10, sample.size = 3000, replace = F)
  d_sub <- sample_data(neon_sub)
  if (nrow(d_sub) <= 3) next
  loc_within_plot <- matrix(as.numeric(unlist(strsplit(gsub(".*[OM]-(\\d+\\.?\\d*-\\d+\\.?\\d*)-.*", "\\1", d_sub$geneticSampleID), "-"))), ncol = 2, byrow = T)
  row.names(loc_within_plot) <- row.names(d_sub)
  dist_within_plot <- dist(loc_within_plot)
  beta_within_plot <- vegdist(otu_table(neon_sub), method = "jaccard")
  png(paste("Fig/",site_i,".png", sep = ""))
  plot(dist_within_plot, beta_within_plot, xlab = "distance", ylab = "jaccard distance",
       main = site_i)
  dev.off()
}
