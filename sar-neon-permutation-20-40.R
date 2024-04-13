# at the 20 by 20 m spatial scale

set.seed=(1203)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot20ID)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot20ID==a4[i])
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

##  get the richness for each subplot at the 10 by 10 spatial scale

richness_mean_subplot20=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_mean_subplot20[i,1]=richness[[i]] 
    richness_mean_subplot20[i,2]=a4[i]
  }
  
  else
  {
    richness_mean_subplot20[i,1]=mean(richness[[i]])
    richness_mean_subplot20[i,2]=a4[i]
  }
}

richness_mean_subplot20=mutate(richness_mean_subplot20,plotid=substr(richness_mean_subplot20$ncol,1,8))

names(richness_mean_subplot20)[1]="richness"

# the mean richness for each 10 x 10 subplot

richness_mean_subplot20=aggregate(richness~plotid,data=richness_mean_subplot20,FUN=mean,na.rm=TRUE)

richness_sd_subplo20=aggregate(richness~plotid,data=richness_mean_subplot20,FUN=sd,na.rm=TRUE)

## at the 30 by 30 m scale

set.seed=(1204)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_A)

mm=str_detect(a4," - 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_A==a4[i])
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

## get the richness in each 30 by 30 subplot estimated with the first approach

richness_mean_subplot30A=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_mean_subplot30A[i,1]=richness[[i]] 
    richness_mean_subplot30A[i,2]=a4[i]
  }
  
  else
  {
    richness_mean_subplot30A[i,1]=mean(richness[[i]])
    richness_mean_subplot30A[i,2]=a4[i]
  }
}

richness_mean_subplot30A=mutate(richness_mean_subplot30A,plotid=substr(richness_mean_subplot30A$ncol,1,8))

names(richness_mean_subplot30A)[1]="richness"

# the mean richness for each 10 x 10 subplot

richness_sd_subplo30A=aggregate(richness~plotid,data=richness_mean_subplot30A,FUN=sd,na.rm=TRUE)

richness_mean_subplot30A=aggregate(richness~plotid,data=richness_mean_subplot30A,FUN=mean,na.rm=TRUE)


# for the second approach

set.seed=(1205)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_B)

mm=str_detect(a4," - 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_B==a4[i])
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

## get the richness in each 30 by 30 subplot estimated with the first approach

richness_mean_subplot30B=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_mean_subplot30B[i,1]=richness[[i]] 
    richness_mean_subplot30B[i,2]=a4[i]
  }
  
  else
  {
    richness_mean_subplot30B[i,1]=mean(richness[[i]])
    richness_mean_subplot30B[i,2]=a4[i]
  }
}

richness_mean_subplot30B=mutate(richness_mean_subplot30B,plotid=substr(richness_mean_subplot30B$ncol,1,8))
names(richness_mean_subplot30B)[1]="richness"

# the mean richness for each 10 x 10 subplot
richness_sd_subplo30B=aggregate(richness~plotid,data=richness_mean_subplot30B,FUN=sd,na.rm=TRUE)
richness_mean_subplot30B=aggregate(richness~plotid,data=richness_mean_subplot30B,FUN=mean,na.rm=TRUE)

###


set.seed=(1206)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_C)

mm=str_detect(a4," - 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_C==a4[i])
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

## get the richness in each 30 by 30 subplot estimated with the first approach

richness_mean_subplot30C=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_mean_subplot30C[i,1]=richness[[i]] 
    richness_mean_subplot30C[i,2]=a4[i]
  }
  
  else
  {
    richness_mean_subplot30C[i,1]=mean(richness[[i]])
    richness_mean_subplot30C[i,2]=a4[i]
  }
}

richness_mean_subplot30C=mutate(richness_mean_subplot30C,plotid=substr(richness_mean_subplot30C$ncol,1,8))
names(richness_mean_subplot30C)[1]="richness"

# the mean richness for each 10 x 10 subplot
richness_sd_subplo30C=aggregate(richness~plotid,data=richness_mean_subplot30C,FUN=sd,na.rm=TRUE)
richness_mean_subplot30C=aggregate(richness~plotid,data=richness_mean_subplot30C,FUN=mean,na.rm=TRUE)

## for the fourth method

set.seed=(1206)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_D)

mm=str_detect(a4," - 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_D==a4[i])
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

## get the richness in each 30 by 30 subplot estimated with the first approach

richness_mean_subplot30D=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_mean_subplot30D[i,1]=richness[[i]] 
    richness_mean_subplot30D[i,2]=a4[i]
  }
  
  else
  {
    richness_mean_subplot30D[i,1]=mean(richness[[i]])
    richness_mean_subplot30D[i,2]=a4[i]
  }
}

richness_mean_subplot30D=mutate(richness_mean_subplot30D,plotid=substr(richness_mean_subplot30D$ncol,1,8))
names(richness_mean_subplot30D)[1]="richness"

# the mean richness for each 30 x 30 subplot
richness_sd_subplo30D=aggregate(richness~plotid,data=richness_mean_subplot30D,FUN=sd,na.rm=TRUE)
richness_mean_subplot30D=aggregate(richness~plotid,data=richness_mean_subplot30D,FUN=mean,na.rm=TRUE)

