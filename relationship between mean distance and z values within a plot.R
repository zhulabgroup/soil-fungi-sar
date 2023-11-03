# determine the mean distance among pairwise soil cores within the 40 x 40 m plot
# we can look at how this mean distance impact the simulated z values
rm(list = ls())
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon") != "AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon") != "OH") # the data only include the O and M soil horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all <- rarefy_even_depth(neon_dob, rngseed = 10, sample.size = 3000, replace = F) # 764 samples were removed
head(sample_data(rare_all))

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

# for the neon data
sub_neon <- subset_samples(rare_all, get_variable(rare_all, "Project") == "NEON")
sub_neon <- subset_samples(sub_neon, get_variable(sub_neon, "horizon") != "AH")
sub_neon <- subset_samples(sub_neon, get_variable(sub_neon, "horizon") != "OH")

# the location within the 40 x 40 m plot
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
a1 <- unique(a1$plotIDM) # 472 plotID


ddt <- numeric()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples in a given plotID
  location <- sample_data(data_sub)[, c("gx", "gy")]
  dis_tance <- dist(location, diag = TRUE, upper = TRUE)
  ddt[i] <- mean(matrix(dis_tance))
}

plot_mean_dist <- data.frame(cbind(plotID = a1, distance = ddt))
plot_mean_dist$distance <- as.numeric(plot_mean_dist$distance)

# does the mean distance affect the estimated z values?

d <- aggregate(z ~ plotID, data = model_data, FUN = mean) # the full dataset includes many variables including the z
d <- merge(d, plot_mean_dist, by = "plotID")
plot(z ~ distance, data = subset(d, z < 10))
summary(lm(z ~ distance, data = subset(d, z < 10)))
# weak while significant positive relationship between the z and the z values
# slope:0.002642, and p-value of 0.0495 *
