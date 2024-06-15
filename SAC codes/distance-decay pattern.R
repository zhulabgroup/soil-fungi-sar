## look at the distance-decay pattern within the 40 by 40 m plot based on the relative location of the soil cores
rm(list=ls())
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="OH")# the data only include the O and M soil horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob<- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob<- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all=rarefy_even_depth(neon_dob, rngseed=10,sample.size = 3000, replace = F)#764 samples were removed
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

# the location within the 40 by 40 m plot
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
a1 <- unique(a1$plotIDM) # 472 plotID for the neon site



pair <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples in a given plotID
  location <- sample_data(data_sub)[, c("gx", "gy")]
  dis_tance <- dist(location, diag = TRUE, upper = TRUE)
  pair[[i]] <- matrix(dis_tance)
}

# the beta diversity between all pairwise cores within a 40 by 40 plot

beta <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples in a given plotID
  com <- otu_table(data_sub)
  ddis <- vegdist(com, method = "jaccard")
  beta[[i]] <- matrix(ddis) # the spatial distance between two locations
}

# combine the community distance and beta diversity

decay <- list()
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) #
  decay[[i]] <- cbind(pair[[i]], beta[[i]])
}

decay_com <- decay[[1]]
for (i in 2:length(a1)) {
  decay_com <- rbind(decay_com, decay[[i]])
}

decay_com <- data.frame(decay_com)
names(decay_com) <- c("distance", "beta")

# look at how the relationship changes as increasing distance

r <- numeric()
for (i in 1:53)
{
  d <- summary(lm(beta ~ distance, data = subset(decay_com, distance < i)))
  r[i] <- d$coefficients[2, 1]
}
