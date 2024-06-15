# step 4 is to estimate the observed fungal richness for each plot. The resulting output was saved locally.

neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon") != "AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon") != "OH") # the data only include the O and M soil horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all <- rarefy_even_depth(neon_dob, rngseed = 10, sample.size = 3000, replace = F) # 764 samples were removed
save(rare_all, file = "rare_all.Rdata") # save the rarefied data

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

rich <- numeric()
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(d, plotIDM == a1[i])
  rb <- taxa_sums(data_sub)
  rich[i] <- length(rb[rb > 0]) # the observed species richness per plot
}

rich_plot <- cbind(a1, rich)
write.csv(rich_plot, "rich_plot.csv")
