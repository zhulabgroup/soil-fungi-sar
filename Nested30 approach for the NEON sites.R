## Nested method for the NEON site,which means that the prior sites were included when adding the new ones

setwd("/../../../sar")
library(doParallel)
# read in the data
neon_dob <- readRDS("/Users/luowenqi/Desktop/sar/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
rm(neon_dob)
# neon <- subset_samples(neon, horizon == "O")
neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0)# both soil horizons were included 
#rm(neon_dob)

dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")

# neon dataset
d <- sample_data(neon) # sample data data frame
d=data.frame(d)
plotID=substr(d$geneticSampleID,1,8)# get the id for each plot so that we can estimate the SAR based on the plot
d <- sample_data(neon)
plotID=data.frame(plotID)
row.names(plotID)=row.names(d)
plotID=sample_data(plotID)

d<- merge_phyloseq(neon, plotID) # adding the plotID to the sample data

a1= sample_data(d)# the unique of the plotID, we have 476 plots
a1=unique(a1$plotID)

times=30
power.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(d, plotID==a1[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression with at least three sites
  {
    cl <- makeCluster(4)
    registerDoParallel(cl)
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","permute","doParallel")) %dopar%{
      
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) # 
      { 
        
        sample_seq <- shuffle(c(1:dim1[1]))
        
        # take out otu_tab for each site
        otu_tab <- otu_table(neon_sub)
        otu_tab <- matrix(otu_tab, nrow = dim1[1], ncol=63685,byrow  = TRUE)#
        
        species[1] <- sum(otu_tab[sample_seq[1],] > 0)
        
        for (j in 2:(dim1[1]))
        { 
          # take out samples as the sequence in sample_seq
          temp <- colSums(otu_tab[c(sample_seq[1:j]),])
          # count species
          species[j] <- sum(temp > 0)
        }
      }
      
      ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
      
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]# to get the estimated c and z values
      return(c(temp[2,1], temp[1,1],temp[2,4],temp[1,4]))
    }
    stopCluster(cl)
  }
  else
  {
    power.z[[i]]=matrix(NA,nrow=4,ncol = dim1[1])###the 10 is randomly selected, to creat a matrix for the sites with less 3 cores to avoid NULL output
  }
}


neon_z=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for the neon site
for(i in 1:length(a1)){
  neon_z[i,]=power.z[[i]][1,]
}

neon_z=data.frame(a1,neon_z)
write.csv(neon_z,"neon_z.nest30.csv")

neon_c=matrix(nrow=length(a1),ncol=30)# get the 30 simulated c values for the neon site
for(i in 1:length(a1)){
  neon_c[i,]=power.z[[i]][2,]
}

neon_c=data.frame(a1,neon_c)
write.csv(neon_c,"neon_c.nest30.csv")

neon_z_mean_nest=apply(neon_z[,2:31],1,mean)# get the mean value of the estimated z and c
neon_c_mean_nest=apply(neon_c[,2:31],1,mean)# get the mean value of the estimated z and c

neon_z_mean_nest=data.frame(neon_z["a1"],neon_z_mean_nest)
neon_c_mean_nest=data.frame(neon_c["a1"],neon_c_mean_nest)

plot_level_zc_nest=cbind(neon_z_mean_nest,neon_c_mean_nest)

write.csv(plot_level_zc,"plot_level_zc_nest.csv")