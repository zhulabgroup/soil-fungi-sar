## guild specifc for the dob site

tax_table(dob_permutation)%>%data.frame()%>%dplyr::select(ta2)%>%table()%>%as.data.frame()%>%arrange(Freq)
#
#           mycoparasite   403
#         root_endophyte   431
#        animal_parasite   497
#  arbuscular_mycorrhizal  786
#               epiphyte   930
#          plant_pathogen  1675
#  unspecified_saprotroph  1838
#        wood_saprotroph   2220
#       litter_saprotroph  3792
#    animal_endosymbiont   6062
#  nectar/tap_saprotroph   6287
#         soil_saprotroph  8808
#        ectomycorrhizal   9143
######

# to select eight fungal guilds

data_EM <- subset_taxa(dob_permutation, ta2 == "ectomycorrhizal")
data_AM <- subset_taxa(dob_permutation, ta2 == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(dob_permutation, ta2 == "soil_saprotroph")
data_littersap <- subset_taxa(dob_permutation, ta2 == "litter_saprotroph")
data_plapat <- subset_taxa(dob_permutation, ta2 == "plant_pathogen")
data_woodsap <- subset_taxa(dob_permutation, ta2 == "wood_saprotroph")
data_para <- subset_taxa(dob_permutation, ta2%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(dob_permutation, ta2 == "epiphyte")

## for the 

set.seed=(989)
times=30
a4=sample_data(data_EM)
a4=unique(a4$subplotID5)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, subplotID5==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  if(dim1>=2)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:1) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 1)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}

## at the 5x5 m2 scale

richness_subplot5_dob_EM=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot5_dob_EM[i,1]=richness[[i]] 
    richness_subplot5_dob_EM[i,2]=a4[i]
  }
  else
  {
    richness_subplot5_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot5_dob_EM[i,2]=a4[i]
  }
}


# get the mean value for each subplot

richness_subplot5_dob_EM%>%mutate(plotid=substr(richness_subplot5_dob_EM$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness),sd_value = sd(richness))->richness_subplot5_dob_EM
save(richness_subplot5_dob_EM,file="richness_subplot5_dob_EM.RData")



### at the 10m x 10m scale

set.seed=(2202)
times=30
a4=sample_data(data_EM)

a4=unique(a4$subplotID10)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, subplotID10==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>3)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:3) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 3)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<3)
  {
    richness[[i]]=NA
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_subplot10_dob_EM=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_dob_EM[i,1]=richness[[i]] 
    richness_subplot10_dob_EM[i,2]=a4[i]
  }
  else
  {
    richness_subplot10_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot10_dob_EM[i,2]=a4[i]
  }
}


# get the mean value for each subplot

richness_subplot10_dob_EM%>%mutate(plotid=substr(richness_subplot10_dob_EM$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot10_dob_EM

save(richness_subplot10_dob_EM,file="richness_subplot10_dob_EM.RData")

### at the 20 m by 20 m scale

set.seed=(2204)
times=30
a4=sample_data(data_EM)
a4=unique(a4$subplotID20)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, subplotID20==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  if(dim1>4)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:4) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 4)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<4)
  {
    richness[[i]]=NA
  }
  else{
    
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}

richness_subplot20_dob_EM=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot20_dob_EM[i,1]=richness[[i]] 
    richness_subplot20_dob_EM[i,2]=a4[i]
  }
  else
  {
    richness_subplot20_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot20_dob_EM[i,2]=a4[i]
  }
}


# get the mean value for each subplot

richness_subplot20_dob_EM%>%mutate(plotid=substr(richness_subplot20_dob_EM$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot20_dob_EM
save(richness_subplot20_dob_EM,file="richness_subplot20_dob_EM.RData")

## at the 30m x 30 m scale

set.seed=(5566)
times=30
a4=sample_data(data_EM)

a4=unique(a4$point_assign1)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, point_assign1==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>6)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:6) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 6)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<6)
  {
    richness[[i]]=NA
  }
  else{
    
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}
### bind the data 

richness_subplot30A_dob_EM=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30A_dob_EM[i,1]=richness[[i]] 
    richness_subplot30A_dob_EM[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30A_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot30A_dob_EM[i,2]=a4[i]
  }
}


richness_subplot30A_dob_EM%>%mutate(plotid=substr(richness_subplot30A_dob_EM$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30A_dob_EM

save(richness_subplot30A_dob_EM,file="richness_subplot30A_dob_EM.RData")


### for the second approach 

set.seed=(5566)
times=30
a4=sample_data(data_EM)
a4=unique(a4$point_assign2)
mm=str_detect(a4," * 1")
mm[!mm]=NA
which(!is.na(mm))
richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, point_assign2==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>6)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:6) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 6)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<6)
  {
    richness[[i]]=NA
  }
  else{
    
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}

### bind the data 

richness_subplot30B_dob_EM=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30B_dob_EM[i,1]=richness[[i]] 
    richness_subplot30B_dob_EM[i,2]=a4[i]
  }
  else
  {
    richness_subplot30B_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot30B_dob_EM[i,2]=a4[i]
  }
}


richness_subplot30B_dob_EM%>%mutate(plotid=substr(richness_subplot30B_dob_EM$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30B_dob_EM

save(richness_subplot30B_dob_EM,file="richness_subplot30B_dob_EM.RData")


###
# in this case, we have only "BC5 * 4" and we do not have "BC5 * 1", suggesting that no cores were assigned within the zone
set.seed=(5566)
times=30
a4=sample_data(data_EM)
a4=unique(a4$point_assign3)
mm=str_detect(a4," * 1")
mm[!mm]=NA
which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, point_assign3==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>6)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:6) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 6)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<6)
  {
    richness[[i]]=NA
  }
  else{
    
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


### bind the data 

richness_subplot30C_dob_EM=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30C_dob_EM[i,1]=richness[[i]] 
    richness_subplot30C_dob_EM[i,2]=a4[i]
  }
  else
  {
    richness_subplot30C_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot30C_dob_EM[i,2]=a4[i]
  }
}


richness_subplot30C_dob_EM%>%mutate(plotid=substr(richness_subplot30C_dob_EM$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30C_dob_EM

save(richness_subplot30C_dob_EM,file="richness_subplot30C_dob_EM.RData")

###
set.seed=(5566)
times=30
a4=sample_data(data_EM)

a4=unique(a4$point_assign4)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, point_assign4==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>6)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:6) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 6)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<6)
  {
    richness[[i]]=NA
  }
  else{
    
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}

### bind the data 

richness_subplot30D_dob_EM=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30D_dob_EM[i,1]=richness[[i]] 
    richness_subplot30D_dob_EM[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30D_dob_EM[i,1]=mean(richness[[i]])
    richness_subplot30D_dob_EM[i,2]=a4[i]
  }
}


richness_subplot30D_dob_EM%>%mutate(plotid=substr(richness_subplot30D_dob_EM$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30D_dob_EM

save(richness_subplot30D_dob_EM,file="richness_subplot30D_dob_EM.RData")

# at the 40 by 40 m scale

set.seed=(5679)

times=30
a4=sample_data(data_EM)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, plotIDM==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>=10)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:10) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 10)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob_EM=matrix(ncol=1,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob_EM[i,]=mean(richness[[i]])
}

richness_mean_subplot40_dob_EM=cbind(plotid=a4,richness=richness_mean_subplot40_dob_EM,rep(1600,43))%>%data.frame()%>%rename_all(~paste0(c("plotid","richness","area")))

## 
save(richness_mean_subplot40_dob_EM,file="richness_mean_subplot40_dob_EM.RData")





