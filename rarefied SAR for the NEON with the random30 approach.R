rm(list=ls())
# rarefy data#
neon_dob <- readRDS("/Users/luowenqi/Desktop/sar/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
rm(neon_dob)
# neon <- subset_samples(neon, horizon == "O")
neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0)# both soil horizons were included 
#rm(neon_dob)
# suggestions for the sampling depth for rarefy. I choose 1000 reads, which will lead to the reoval of 10% samples
#https://forum.qiime2.org/t/how-to-select-a-sampling-depth/4265/5
#https://forum.qiime2.org/t/what-is-the-threshold-above-that-i-select-from-the-sequence-count-for-diversity-analysis-in-qiime-2/6897
#dim(subset(k,k<3000))[1]/6222
rare_neon=rarefy_even_depth(neon, rngseed=10,sample.size = 3000, replace = F)#749 samples were removed
d <- sample_data(neon_rare) # sample data data frame
d=data.frame(d)
plotID=substr(d$geneticSampleID,1,8)# get the id for each plot so that we can estimate the SAR based on the plot
d <- sample_data(neon_rare)
plotID=data.frame(plotID)
row.names(plotID)=row.names(d)
plotID=sample_data(plotID)
d<- merge_phyloseq(neon_rare, plotID) # adding the plotID to the sample data
a1= sample_data(d)# the unique of the plotID, we have 476 plots
a1=unique(a1$plotID)
times=30
power.z <- vector("list", length(a1))
for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_rare_sub <- subset_samples(d, plotID==a1[i])
  dim1 <- dim(otu_table(neon_rare_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression with at least three sites
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) # Random All ?# yes, because we do not need to control the sample numbers within each plot
      { 
        # randomly sample j samples in the site 
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], j)] <- TRUE
        temp <- merge_samples(neon_rare_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
      }
      ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]# to get the estimated c and z values
      return(c(temp[2,1], temp[1,1]))
    }
    stopCluster(cl)
  }
  else
  {
    power.z[[i]]=matrix(10,nrow=2,ncol = dim1[1])###the 10 is randomly selected, to creat a matrix for the sites with less 3 cores to avoid NULL output
  }
}

rare_neon_z=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for the rare_neon site
for(i in 1:length(a1)){
  rare_neon_z[i,]=power.z[[i]][1,]
}
rare_neon_z=data.frame(a1,rare_neon_z)
rare_neon_c=matrix(nrow=length(a1),ncol=30)# get the 30 simulated c values for the rare_neon site
for(i in 1:length(a1)){
  rare_neon_c[i,]=power.z[[i]][2,]
}
rare_neon_c=data.frame(a1,rare_neon_c)
rare_neon_z_mean=apply(rare_neon_z[,2:31],1,mean)# get the mean value of the estimated z and c
rare_neon_c_mean=apply(rare_neon_c[,2:31],1,mean)# get the mean value of the estimated z and c
rare_neon_z_mean=data.frame(rare_neon_z["a1"],rare_neon_z_mean)
rare_neon_c_mean=data.frame(rare_neon_c["a1"],rare_neon_c_mean)
plot_level_zc=cbind(rare_neon_z_mean,rare_neon_c_mean)
write.csv(plot_level_zc,"plot_level_zc_rare.csv")
