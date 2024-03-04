# to show the variation in the estimated z values, taking a plot for an example
library(doParallel)
library(permute)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(cowplot)
# rarefy all the full data
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon") != "AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon") != "OH") # the data only include the O and M soil horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all <- rarefy_even_depth(neon_dob, rngseed = 10, sample.size = 3000, replace = F) # 764 samples were removed
#
d <- sample_data(rare_all)
table(d$Project) # the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1 <- data.frame(d[, c("geneticSampleID", "Site")])
plotID <- substr(d1$geneticSampleID, 1, 8)
d1 <- cbind(d1, plotID)
iddob <- d1$Site[1:908] # a site corresponded to a plot
idneon <- d1$plotID[909:6378] # an unique plotID corresponds to a plot
plotIDM <- data.frame(c(iddob, idneon))
names(plotIDM) <- "plotIDM" # the plot id used for the SAR
row.names(plotIDM) <- row.names(d)
plotIDM <- sample_data(plotIDM)
d <- merge_phyloseq(rare_all, plotIDM) # merge the new plotid with the initial data
# select an unique plot and build a SAR within the plot
a1 <- sample_data(d) # the unique plotID, we have 515 plots
a1 <- unique(a1$plotIDM)

### when i=1, we can get the z values for a dob site. i stopped computing after running i=1
times <- 30
power.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  neon_sub <- subset_samples(d, plotIDM == a1[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  if (dim1[1] >= 3) # to construct a linear regression, at least three cores are required
    {
      cl <- makeCluster(3)
      registerDoParallel(cl)
      power.z[[i]] <- foreach(i = 1:times, .combine = "rbind", .packages = c("phyloseq", "permute", "doParallel")) %dopar% {
        species <- vector(length = dim1[1]) # create a vector to save diversity
        for (j in 1:dim1[1]) #
        {
          sample_seq <- shuffle(c(1:dim1[1]))

          # otu table for each site
          otu_tab <- otu_table(neon_sub)
          otu_tab <- matrix(otu_tab, nrow = dim1[1], ncol = 67769, byrow = FALSE) # 67769 taxa

          species[1] <- sum(otu_tab[sample_seq[1], ] > 0)

          for (j in 2:(dim1[1]))
          {
            # take out samples as the sequence in sample_seq
            temp <- colSums(otu_tab[c(sample_seq[1:j]), ])
            # count species
            species[j] <- sum(temp > 0)
          }
        }

        ex <- as.data.frame(cbind(species, "A" = c(1:dim1[1])))

        return(ex)
      }
      stopCluster(cl)
    } else {
    power.z[[i]] <- matrix(NA, nrow = 2, ncol = dim1[1]) ### the 10 is randomly selected, to creat a matrix for the sites with less 3 cores to avoid NULL output
  }
}

# select a site, when i =1, there are 24 cores within the plot

dd <- rep(1:30, each = 24) #
ex <- cbind(power.z[[1]], dd)
ex_dob <- ex
# to get the mean richness of each area interval
aggregate(species ~ A, data = ex_dob, FUN = mean)
# to get the sd of species richness of each area interval
aggregate(species ~ A, data = ex_dob, FUN = sd)

a <- ggplot(data = ex_dob, aes(x = log(A), y = log(species), color = as.factor(dd))) +
  geom_point(size = 3, alpha = 0.5) +
  ggtitle("a DoB plot with 24 cores") +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
  xlab("Log(Area)") +
  ylab("Log(Species)") +
  guides(color = "none")
geom_smooth(method = "lm", se = FALSE)


b <- ggplot(data = ex_dob, aes(x = log(A), y = log(species), color = as.factor(dd))) +
  geom_point(size = 3, alpha = 0.5) +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
  xlab("Log(Area)") +
  ylab("Log(Species)") +
  guides(color = "none") +
  geom_smooth(method = "lm", se = FALSE)



c <- ggplot(ex_dob, aes(x = as.factor(A), y = species)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2, color = "red") +
  geom_boxplot() +
  xlab("Area") +
  ylab("Species") +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 10, hjust = 0.5), axis.text.y = element_text(hjust = 0, angle = 90), axis.text.x = element_text(hjust = 1, size = 10, angle = 90), axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))


# for the neon site

times <- 30
power.z <- vector("list", length(a1))

for (i in 45:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  neon_sub <- subset_samples(d, plotIDM == a1[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  if (dim1[1] >= 3) # to construct a linear regression, at least three cores are required
    {
      cl <- makeCluster(3)
      registerDoParallel(cl)
      power.z[[i]] <- foreach(i = 1:times, .combine = "rbind", .packages = c("phyloseq", "permute", "doParallel")) %dopar% {
        species <- vector(length = dim1[1]) # create a vector to save diversity
        for (j in 1:dim1[1]) #
        {
          sample_seq <- shuffle(c(1:dim1[1]))

          # otu table for each site
          otu_tab <- otu_table(neon_sub)
          otu_tab <- matrix(otu_tab, nrow = dim1[1], ncol = 67769, byrow = FALSE) # 67769 taxa

          species[1] <- sum(otu_tab[sample_seq[1], ] > 0)

          for (j in 2:(dim1[1]))
          {
            # take out samples as the sequence in sample_seq
            temp <- colSums(otu_tab[c(sample_seq[1:j]), ])
            # count species
            species[j] <- sum(temp > 0)
          }
        }

        ex <- as.data.frame(cbind(species, "A" = c(1:dim1[1])))

        return(ex)
      }
      stopCluster(cl)
    } else {
    power.z[[i]] <- matrix(NA, nrow = 2, ncol = dim1[1]) ### the 10 is randomly selected, to creat a matrix for the sites with less 3 cores to avoid NULL output
  }
}


dd <- rep(1:30, each = 12) #

ex <- cbind(power.z[[45]], dd)
ex_neon <- ex

d <- ggplot(data = ex_neon, aes(x = log(A), y = log(species), color = as.factor(dd))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
  xlab("Log(Area)") +
  ylab("Log(Species)") +
  guides(color = "none") +
  geom_smooth(method = "lm", se = FALSE)


# to get the mean richness of each area interval
aggregate(species ~ A, data = ex_neon, FUN = mean)
# to get the sd of species richness of each area interval
aggregate(species ~ A, data = ex_neon, FUN = sd)
# for a neon site

e <- ggplot(data = ex_neon, aes(x = log(A), y = log(species), color = as.factor(dd))) +
  geom_point(size = 3, alpha = 0.5) +
  ggtitle("a NEON plot with 12 cores") +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA)) +
  xlab("Log(Area)") +
  ylab("Log(Species)") +
  guides(color = "none")

geom_smooth(method = "lm", se = FALSE)

f <- ggplot(ex_neon, aes(x = as.factor(A), y = species)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2, color = "red") +
  geom_boxplot() +
  xlab("Area") +
  ylab("Species") +
  theme(legend.position = "bottom", text = element_text(size = 18), plot.title = element_text(size = 15, hjust = 0.5), axis.text.y = element_text(hjust = 0, angle = 90), axis.text.x = element_text(hjust = 1), axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.ticks.x = element_blank(), panel.background = element_rect(fill = "NA"), panel.border = element_rect(color = "black", size = 1.5, fill = NA))

a <- ggplotGrob(a)

b <- ggplotGrob(b)
c <- ggplotGrob(c)
d <- ggplotGrob(d)
e <- ggplotGrob(e)
f <- ggplotGrob(f)

c$widths <- a$widths

plot_grid(a, e, b, d, c, f, ncol = 2, labels = c("A", "D", "B", "E", "C", "F"), label_size = 20)
