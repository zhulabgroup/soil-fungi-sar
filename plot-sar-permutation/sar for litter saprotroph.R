set.seed=(989)
times=30
a4=sample_data(data_littersap)
a4=unique(a4$subplotID5)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplotID5==a4[i])
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}

##
richness_subplot5_neon_littersap=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot5_neon_littersap[i,1]=richness[[i]] 
    richness_subplot5_neon_littersap[i,2]=a4[i]
  }
  else
  {
    richness_subplot5_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot5_neon_littersap[i,2]=a4[i]
  }
}

# get the mean

# get the mean value for each subplot

richness_subplot5_neon_littersap%>%mutate(plotid=substr(richness_subplot5_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness),sd_value = sd(richness))->richness_subplot5_neon_littersap

save(richness_subplot5_neon_littersap,file="richness_subplot5_neon_littersap.RData")

# for the 10 by 10 m scale

set.seed=(3202)
times=30
a4=sample_data(data_littersap)
a4=unique(a4$subplotID10)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplotID10==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>3)# we select four soil core in each subplot
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
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

richness_subplot10_neon_littersap=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_neon_littersap[i,1]=richness[[i]] 
    richness_subplot10_neon_littersap[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot10_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot10_neon_littersap[i,2]=a4[i]
  }
}
# get the mean

richness_subplot10_neon_littersap%>%mutate(plotid=substr(richness_subplot10_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot10_neon_littersap

save(richness_subplot10_neon_littersap,file="richness_subplot10_neon_littersap.RData")

## at the 20 m by 20 m

set.seed=(1203)
times=30
a4=sample_data(data_littersap)

a4=unique(a4$subplotID20)
length(a4)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplotID20==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>3)
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
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


richness_subplot20_neon_littersap=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot20_neon_littersap[i,1]=richness[[i]] 
    richness_subplot20_neon_littersap[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot20_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot20_neon_littersap[i,2]=a4[i]
  }
}

richness_subplot20_neon_littersap%>%mutate(plotid=substr(richness_subplot20_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot20_neon_littersap

save(richness_subplot20_neon_littersap,file="richness_subplot20_neon_littersap.RData")

## for the 30 by 30 m scale

set.seed=(5204)
times=30
a4=sample_data(data_littersap)
a4=unique(a4$subplot_ID_30_A)
length(a4)

mm=str_detect(a4," * 1")# just selected the 1 subplots

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplot_ID_30_A==a4[i])
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
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

richness_subplot30A_neon_littersap=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30A_neon_littersap[i,1]=richness[[i]] 
    richness_subplot30A_neon_littersap[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30A_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot30A_neon_littersap[i,2]=a4[i]
  }
}

richness_subplot30A_neon_littersap%>%mutate(plotid=substr(richness_subplot30A_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30A_neon_littersap

save(richness_subplot30A_neon_littersap,file="richness_subplot30A_neon_littersap.RData")

# for the second approach

set.seed=(5208)
times=30
a4=sample_data(data_littersap)
a4=unique(a4$subplot_ID_30_B)
mm=str_detect(a4," * 1")
mm[!mm]=NA
which(!is.na(mm))
richness <- vector("list", length(which(!is.na(mm))))
for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplot_ID_30_B==a4[i])
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
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
###

richness_subplot30B_neon_littersap=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30B_neon_littersap[i,1]=richness[[i]] 
    richness_subplot30B_neon_littersap[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30B_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot30B_neon_littersap[i,2]=a4[i]
  }
}

richness_subplot30B_neon_littersap%>%mutate(plotid=substr(richness_subplot30B_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30B_neon_littersap

save(richness_subplot30B_neon_littersap,file="richness_subplot30B_neon_littersap.RData")

##

set.seed=(4206)
times=30
a4=sample_data(data_littersap)

a4=unique(a4$subplot_ID_30_C)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplot_ID_30_C==a4[i])
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
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
###########
richness_subplot30C_neon_littersap=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30C_neon_littersap[i,1]=richness[[i]] 
    richness_subplot30C_neon_littersap[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30C_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot30C_neon_littersap[i,2]=a4[i]
  }
}


richness_subplot30C_neon_littersap%>%mutate(plotid=substr(richness_subplot30C_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30C_neon_littersap

save(richness_subplot30C_neon_littersap,file="richness_subplot30C_neon_littersap.RData")

####

set.seed=(6206)
times=30
a4=sample_data(data_littersap)

a4=unique(a4$subplot_ID_30_D)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplot_ID_30_D==a4[i])
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
        tlittersapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tlittersapp)["TRUE"] > 0)
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

richness_subplot30D_neon_littersap=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30D_neon_littersap[i,1]=richness[[i]] 
    richness_subplot30D_neon_littersap[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30D_neon_littersap[i,1]=mean(richness[[i]])
    richness_subplot30D_neon_littersap[i,2]=a4[i]
  }
}

richness_subplot30D_neon_littersap%>%mutate(plotid=substr(richness_subplot30D_neon_littersap$ncol,1,8))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot30D_neon_littersap

save(richness_subplot30D_neon_littersap,file="richness_subplot30D_neon_littersap.RData")

richness_mean_subplot30_neon_littersap=rbind(richness_subplot30A_neon_littersap,richness_subplot30B_neon_littersap,richness_subplot30C_neon_littersap,richness_subplot30D_neon_littersap)%>%group_by(plotid)%>%
  summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(mean_value,na.rm = TRUE))%>%filter(!is.na(mean_value))%>%mutate(area=rep(900,length.out = n()))

save(richness_mean_subplot30_neon_littersap,file="richness_mean_subplot30_neon_littersap.RData")

# 
