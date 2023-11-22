# test the relationship between the number of soil cores and the estimated z value
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

2 # the number of soil cores per plot

ncore <- numeric()
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(d, plotIDM == a1[i])
  ncore[i] <- dim(otu_table(data_sub))[1]
}

ncore <- cbind(ncore, plotID = a1) %>% data.frame()

ncore <- ncore[, -3]
# 3 load the simulated 1000 data:

z <- aggregate(z ~ plotID, data = data_sen_mod, FUN = mean)

names(ncore) <- c("plotID", "z")
ncore$ncore <- as.numeric(ncore$ncore)
z <- merge(z, ncore, by = "plotID")
# create a plot
ggplot(data = z, aes(x = log(ncore), y = z, color = project)) +
  geom_point(data = z, aes(x = log(ncore), y = z, color = project), size = 2) +
  theme(legend.key.size = unit(0.15, "inches"), legend.position = c(0.5, 0.85), legend.text.align = 0, panel.background = element_blank(), panel.border = element_rect(fill = NA, size = 1, color = "black"), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15), plot.title = element_text(hjust = 0.5, face = "bold", size = 18), axis.text.x = element_text(size = 15)) +
  xlab("Number of cores") +
  ylab(expression(italic(z) * " value")) +
  geom_smooth(method = "lm") +
  scale_color_manual("", breaks = c("FALSE", "TRUE"), labels = c("DoB", "NEON"), values = c("black", "mediumpurple"))

# does this is because different projects

project <- str_detect(z$plotID, "_")
z <- cbind(z, project)
