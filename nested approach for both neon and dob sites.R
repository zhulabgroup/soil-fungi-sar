# load the rarefied data first
library(doParallel)
library(permute)

# rarefy all the data
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="OH")# the data only include the O and M soil horizon

neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob<- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob<- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all=rarefy_even_depth(neon_dob, rngseed=10,sample.size = 3000, replace = F)#764 samples were removed

# 
d=sample_data(rare_all)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site correspondes to a plot
idneon=d1$plotID[909:6378]# an unique plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))
names(plotIDM)="plotIDM"# the plot id used for the SAR
row.names(plotIDM)=row.names(d)
plotIDM=sample_data(plotIDM)
d<- merge_phyloseq(rare_all, plotIDM)# merge the new plotid with the initial data 
# select an unique plot and build a SAR within the plot
a1= sample_data(d)# the unique plotID, we have 476 plots
a1=unique(a1$plotIDM)

###
times <- 30
power.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  neon_sub <- subset_samples(d, plotIDM == a1[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  if (dim1[1] >= 3) # we can construct a linear regression with at least three sites
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq", "permute", "doParallel")) %dopar% {
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) #
      {
        sample_seq <- shuffle(c(1:dim1[1]))
        
        # take out otu_tab for each site
        otu_tab <- otu_table(neon_sub)
        otu_tab <- matrix(otu_tab, nrow = dim1[1], ncol = 67769, byrow = FALSE) #
        
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
      
      temp <- summary(lm(log(species) ~ log(A), ex))[["coefficients"]] # to get the estimated c and z values
      return(c(temp[2, 1], temp[1, 1], temp[2, 4], temp[1, 4]))
    }
    stopCluster(cl)
  } else {
    power.z[[i]] <- matrix(NA, nrow = 4, ncol = dim1[1]) ### the 10 is randomly selected, to creat a matrix for the sites with less 3 cores to avoid NULL output
  }
}



neon_z <- matrix(nrow = length(a1), ncol = 30) # get the 30 simulated z values for the neon site
for (i in 1:length(a1)) {
  neon_z[i, ] <- power.z[[i]][1, ]
}
# 1:12，12#


