# rarefy all the data
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="OH")# the data only include the O and M soil horizon

neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob<- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob<- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all=rarefy_even_depth(neon_dob, rngseed=10,sample.size = 3000, replace = F)#764 samples were removed
save(rare_all,file="rare_all.Rdata")# save the rarefied data

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

set.seed(1010)
times=30
power.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression with at least three cores
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], j)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
      }
      ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]# extract the estimated c and z values
      return(c(temp[2,1], temp[1,1]))
    }
    stopCluster(cl)
  }
  
  else
  {
    power.z[[i]]=matrix(10,nrow=2,ncol = dim1[1])#the number 10 is randomly selected, to create a matrix for the plots with < 3 cores to avoid NULL output
  }
}

all_z=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for all the plots
for(i in 1:length(a1)){
  all_z[i,]=power.z[[i]][1,]
}

all_z=data.frame(a1,all_z)
write.csv(all_z,"all_z.ranall.csv")

all_c=matrix(nrow=length(a1),ncol=30)# get the 30 simulated c values for the neon site
for(i in 1:length(a1)){
  all_c[i,]=power.z[[i]][2,]
}

all_c=data.frame(a1,all_c)
write.csv(all_c,"all_c.ranall.csv")
