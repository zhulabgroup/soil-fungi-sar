library(doParallel)
library(stringr)
library(dplyr)
library(phyloseq)

# at the 20 by 20 m spatial scale

# Define the dimensions of each grid cell
grid_cell_size <- 20

# Create grid boundaries
x_breaks <- seq(0, plot_width, by = grid_cell_size)
y_breaks <- seq(0, plot_height, by = grid_cell_size)

# Assign a ID to each point based on their coordinates
subplot_ID$gx=as.numeric(subplot_ID$gx)
subplot_ID$gy=as.numeric(subplot_ID$gy)
subplot_ID$gx <- cut(subplot_ID$gx, breaks = x_breaks, labels = FALSE,include.lowest = TRUE)
subplot_ID$gy <- cut(subplot_ID$gy, breaks = y_breaks, labels = FALSE,include.lowest = TRUE)

# Display the result
head(subplot_ID)

###

3# add the subplotID to the phyloseq object

subplotID20=paste(subplot_ID$gx,"-",subplot_ID$gy)%>%data.frame()
plotIDM=sample_data(a1)%>%data.frame()
plotIDM=plotIDM$plotIDM
subplotID20=paste(plotIDM,"*",subplotID20$.)
subplotID20=data.frame(subplotID20)
names(subplotID20)="subplotID20"
row.names(subplotID20)=row.names(sample_data(a1))
subplotID20=sample_data(subplotID20)
a1=merge_phyloseq(a1,subplotID20)

###

set.seed=(1203)
times=30
a4=sample_data(a1)

a4=unique(a4$subplotID20)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplotID20==a4[i])
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

##  get the richness for each subplot at the 20 by 20 spatial scale


richness_subplot20_neon=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot20_neon[i,1]=richness[[i]] 
    richness_subplot20_neon[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot20_neon[i,1]=mean(richness[[i]])
    richness_subplot20_neon[i,2]=a4[i]
  }
}

richness_subplot20_neon=mutate(richness_subplot20_neon,plotid=substr(richness_subplot20_neon$ncol,1,8))

names(richness_subplot20_neon)[1]="richness"

# the mean richness for each 20 x 20 subplot

richness_mean_subplot20_neon=aggregate(richness~plotid,data=richness_subplot20_neon,FUN=mean,na.rm=TRUE)

richness_sd_subplo20_neon=aggregate(richness~plotid,data=richness_subplot20_neon,FUN=sd,na.rm=TRUE)


## at the 30 by 30 m scale

subplot_ID=sample_data(a1)%>%data.frame()

subplot_ID=subplot_ID[,c("gx","gy")]


point_assign1=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=0 && rowid$gx<=30 && rowid$gy>=0 && rowid$gy<=30)
  {
    point_assign1[i]=1
  }
  else
  {
    point_assign1[i]=4
  }
}

# the second approach
point_assign2=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=0 && rowid$gx<=30 && rowid$gy>=10 && rowid$gy<=40)
  {
    point_assign2[i]=1
  }
  else
  {
    point_assign2[i]=4
  }
}

# the third approach
point_assign3=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=10 && rowid$gx<=40 && rowid$gy>=10 && rowid$gy<=40)
  {
    point_assign3[i]=1
  }
  else
  {
    point_assign3[i]=4
  }
}

# the fourth approach

point_assign4=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=10 && rowid$gx<=40 && rowid$gy>=0 && rowid$gy<=30)
  {
    point_assign4[i]=1
  }
  else
  {
    point_assign4[i]=4
  }
}

point_assign1=paste(plotIDM,"*",point_assign1)
point_assign2=paste(plotIDM,"*",point_assign2)
point_assign3=paste(plotIDM,"*",point_assign3)
point_assign4=paste(plotIDM,"*",point_assign4)

four_assign=cbind(point_assign1,point_assign2,point_assign3,point_assign4)%>%data.frame()
names(four_assign)=c("subplot_ID_30_A","subplot_ID_30_B","subplot_ID_30_C","subplot_ID_30_D")
row.names(four_assign)=row.names(sample_data(a1))
four_assign=sample_data(four_assign)

a1=merge_phyloseq(a1,four_assign)


## need to assign the plotID
#for the first approach
set.seed=(5204)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_A)

mm=str_detect(a4," * 1")

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

richness_subplot30A_neon=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30A_neon[i,1]=richness[[i]] 
    richness_subplot30A_neon[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30A_neon[i,1]=mean(richness[[i]])
    richness_subplot30A_neon[i,2]=a4[i]
  }
}

richness_subplot30A_neon=mutate(richness_subplot30A_neon,plotid=substr(richness_subplot30A_neon$ncol,1,8))

names(richness_subplot30A_neon)[1]="richness"

# the mean richness for each 30 x 30 subplot

richness_sd_subplo30A_neon=aggregate(richness~plotid,data=richness_subplot30A_neon,FUN=sd,na.rm=TRUE)

richness_mean_subplot30A_neon=aggregate(richness~plotid,data=richness_subplot30A_neon,FUN=mean,na.rm=TRUE)


# for the second approach

set.seed=(5208)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_B)

mm=str_detect(a4," * 1")

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

richness_subplot30B_neon=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30B_neon[i,1]=richness[[i]] 
    richness_subplot30B_neon[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30B_neon[i,1]=mean(richness[[i]])
    richness_subplot30B_neon[i,2]=a4[i]
  }
}

richness_subplot30B_neon=mutate(richness_subplot30B_neon,plotid=substr(richness_subplot30B_neon$ncol,1,8))

names(richness_subplot30B_neon)[1]="richness"

# the mean richness for each 30 x 30 subplot

richness_mean_subplot30B_neon=aggregate(richness~plotid,data=richness_subplot30B_neon,FUN=mean,na.rm=TRUE)
richness_sd_subplo30B_neon=aggregate(richness~plotid,data=richness_subplot30B_neon,FUN=sd,na.rm=TRUE)

###

# for the third approach



set.seed=(4206)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_C)

mm=str_detect(a4," * 1")

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

richness_subplot30C_neon=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30C_neon[i,1]=richness[[i]] 
    richness_subplot30C_neon[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30C_neon[i,1]=mean(richness[[i]])
    richness_subplot30C_neon[i,2]=a4[i]
  }
}

richness_subplot30C_neon=mutate(richness_subplot30C_neon,plotid=substr(richness_subplot30C_neon$ncol,1,8))
names(richness_subplot30C_neon)[1]="richness"

# the mean richness for each 10 x 10 subplot
richness_sd_subplo30C_neon=aggregate(richness~plotid,data=richness_subplot30C_neon,FUN=sd,na.rm=TRUE)
richness_mean_subplot30C_neon=aggregate(richness~plotid,data=richness_subplot30C_neon,FUN=mean,na.rm=TRUE)

## for the fourth method

set.seed=(6206)
times=30
a4=sample_data(a1)

a4=unique(a4$subplot_ID_30_D)

mm=str_detect(a4," * 1")

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

richness_subplot30D_neon=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30D_neon[i,1]=richness[[i]] 
    richness_subplot30D_neon[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30D_neon[i,1]=mean(richness[[i]])
    richness_subplot30D_neon[i,2]=a4[i]
  }
}

richness_subplot30D_neon=mutate(richness_subplot30D_neon,plotid=substr(richness_subplot30D_neon$ncol,1,8))
names(richness_subplot30D_neon)[1]="richness"

# the mean richness for each 10 x 10 subplot
richness_sd_subplo30D_neon=aggregate(richness~plotid,data=richness_subplot30D_neon,FUN=sd,na.rm=TRUE)
richness_mean_subplot30D_neon=aggregate(richness~plotid,data=richness_subplot30D_neon,FUN=mean,na.rm=TRUE)

# bind the richness data determined with three different approaches

richness_subplot30_neon=rbind(richness_mean_subplot30A_neon,richness_mean_subplot30B_neon,richness_mean_subplot30C_neon,richness_mean_subplot30D_neon)

richness_mean_subplot30_neon=aggregate(richness~plotid,data=richness_subplot30_neon,FUN=mean)
richness_sd_subplot30_neon=aggregate(richness~plotid,data=richness_subplot30_neon,FUN=sd)

richness_mean_subplot30_neon=mutate(richness_mean_subplot30_neon,area=rep(900,dim(richness_mean_subplot30_neon)[1]))


### at the 40 by 40 plot

6.# get the plot level richness, should we include all the richness all a subset of cores

a4=sample_data(a1)
a4=unique(a4$plotIDM)

richness=numeric()
con_40=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
  con_40[i]=dim(sample_data(data_sub))[1]
  total_richness=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  richness[i]=total_richness
}

richness_40plot=data.frame(a4,richness)%>%data.frame()

names(richness_40plot)[1]="plotID"

# if we included a number of cores in each plot

set.seed=(5677)
times=30
a4=sample_data(a1)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
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
    richness[[i]]=NA
  }
}

# to see the number of soil core for each neon plot

a4=sample_data(a1)
a4=unique(a4$plotIDM)

con_40 <-numeric()

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
  con_40[i]=dim(sample_data(data_sub))[1]
  
}





## to see how many cores are there for one site


a4=sample_data(a1)
a4=unique(a4$siteID)

con_40 <-numeric()

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, siteID==a4[i])
  con_40[i]=dim(sample_data(data_sub))[1]
  
}


# how many plot per site

a4=sample_data(a1)
a4=unique(a4$siteID)

plot_number_persite <-numeric()

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, siteID==a4[i])
  
  plot_number_persite[i]=dim(data.frame(sample_data(data_sub))["plotIDM"]%>%unique())[1]
  
}
plot_number_persite=data.frame(plot_number_persite,a4)

names(plot_number_persite)=c("plotN","a")


richness_subplot40_neon=matrix(ncol=1,nrow=472)
for (i in 1:472)
{
  
  richness_subplot40_neon[i,]=mean(richness[[i]])
}

richness_mean_subplot40_neon=cbind(plotid=a4,richness=richness_subplot40_neon,rep(1600,472))%>%data.frame()

richness_mean_subplot40_neon$V3=as.numeric(richness_mean_subplot40_neon$V3)

names(richness_mean_subplot40_neon)=c("plotid","richness","area")

# combing the richness at different scales

richness_mean_subplot5_neon=mutate(richness_mean_subplot5_neon,area=rep(25,dim(richness_mean_subplot5_neon)[1]))
richness_mean_subplot10_neon=mutate(richness_mean_subplot10_neon,area=rep(100,dim(richness_mean_subplot10_neon)[1]))
richness_mean_subplot20_neon=mutate(richness_mean_subplot20_neon,area=rep(400,dim(richness_mean_subplot20_neon)[1]))
richness_mean_subplot30_neon=mutate(richness_mean_subplot30_neon,area=rep(900,dim(richness_mean_subplot30_neon)[1]))

richness_mean_subplot30_neon=mutate(richness_mean_subplot30_neon,area=rep(900,dim(richness_mean_subplot30_neon)[1]))


sar_neon_permutation=rbind(richness_mean_subplot5_neon,richness_mean_subplot10_neon,richness_mean_subplot20_neon,richness_mean_subplot30_neon,richness_mean_subplot40_neon)


save(sar_neon_permutation,file="sar_neon_permutation.RData")



# get the site level richness
# for each site we randomly select 100 cores

set.seed=(5800)
times=30
a4=sample_data(a1)
a4=unique(a4$siteID)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, siteID==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>=100)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:100) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 100)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=NA
  }
}


richness_site_neon=data.frame(nrow=length(a4),ncol=3)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_site_neon[i,1]=richness[[i]] 
    richness_site_neon[i,2]=a4[i]
    richness_site_neon[i,3]=0
  }
  
  else
  {
    richness_site_neon[i,1]=mean(richness[[i]])
    richness_site_neon[i,2]=a4[i]
    richness_site_neon[i,3]=sd(richness[[i]])
  }
}

names(richness_site_neon)=c("richness_mean","a","richness_sd")

# get the richness predicted for these sites based on the SAR function


zvalue=numeric()
pvalue=numeric()
c=numeric()
a5=unique(point_number_three_neon$Var1)
for (i in 1:length(a5))
{
  df1=subset(sar_neon_permutation,plotid==a5[i])
  if(dim(df1)[1]<3){
    zvalue[i]=NA
    pvalue[i]=NA
    c[i]=NA
  }
  else{
    ft=lm(log(richness)~log(area),data=df1)%>%summary()
    zvalue[i]=ft$coefficients[2,1]
    pvalue[i]=ft$coefficients[1,4]
    c[i]=ft$coefficients[1]
  }
  
}

df=cbind(data.frame(a5),zvalue,c)

df$zvalue=as.numeric(df$zvalue)

names(df)=c("plotid","zvalue","logc")

df=mutate(df,a=substr(plotid,1,4))

# all the z values determined with more than three dots 

names(point_number_neon)[1]="plotid"
df=merge(df,point_number_neon,by="plotid")
z_neon=df
# to get the mean value of each log and z for each site

site_pred_rich_sd=aggregate( pre_rich~a,data=z_neon,FUN=sd)
site_pred_rich_mean=aggregate( pre_rich~a,data=z_neon,FUN=mean)

data.frame(table(z_neon$a))

site_area_10=site_area_10[,-c(11:14)]

# add the observed richness to the data

site_area_10=merge(site_area_10,richness_site_neon,by="a")

site_pred_rich_mean=cbind(site_pred_rich_mean,site_pred_rich_sd)
site_pred_rich_mean=site_pred_rich_mean[,-3]
names(site_pred_rich_mean)=c("a","pred_richness","pred_sd")

site_area_10=merge(site_area_10,site_pred_rich_mean,by="a")

site_area_10[,-c("richness")]

save(site_area_10,file="site_area_10.RData")

p1=ggplot()+
  geom_segment(data=site_area_10,color="gray",aes(x=pred_richness,y=log(richness_mean)-log(richness_sd/10),yend=log(richness_mean)+log(richness_sd/10),xend=pred_richness))+

geom_segment(data=site_area_10,color="gray",aes(y=log(richness_mean),yend=log(richness_mean),x=pred_richness-pred_sd, xend=pred_richness+pred_sd))+
  xlab("Predicted site-level richness")+
  ylab("Observed site-level richness")+
  theme(legend.position ="right", 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  geom_smooth(data=site_area_10,aes(x=pred_richness,y=log(richness_mean)),method="lm")+
  geom_point(data=site_area_10,pch=21,alpha=0.8,color="black",aes(x=pred_richness,y=log(richness_mean),fill=a),size=4)+
  scale_fill_manual("Sites",breaks=site_area_10$a,values=custom_colors)+
  annotate("text",x=7.25,y=11.5,label=expression(italic(R^2)*"=77%***"),size=6)+
  guides(fill="none")

  
  # if we just use the mean of the logc and the mean z value for a site
meanz=aggregate(zvalue~a,data=z_neon,FUN=mean)
meanc=aggregate(logc~a,data=z_neon,FUN=mean)

mean_z_c=cbind(meanz,meanc)[,-3]
names(mean_z_c)=c("a","meanz","meanc")


p2=ggplot()+
  geom_point(data=site_area_10,pch=21,aes(x=pred_rich2,y=log(richness_mean),fill=a),size=4)+
  geom_segment(data=site_area_10,color="gray",aes(x=pred_rich2,y=log(richness_mean)-log(richness_sd/10),yend=log(richness_mean)+log(richness_sd/10),xend=pred_rich2))+
  scale_fill_manual("Sites",breaks=site_area_10$a,values=custom_colors)+
  geom_smooth(data=site_area_10,pch=21,aes(x=pred_rich2,y=log(richness_mean)),method="lm")+
  
  xlab("Predicted site-level richness")+
  ylab("")+
  theme(legend.position ="right", 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  annotate("text",x=9.5,y=11.5,label=expression(italic(R^2)*"=47%**"),size=6)+
  guides(fill="none")
  
  

###
library(RColorBrewer)

custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                                   "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                                   "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                                   "#E5C494", "#B3B3B3", "#8DD3C7")

p3=ggplot()+
  geom_point(data=site_area_10,pch=21,aes(x=pred_rich2,y=log(obrich),fill=a),size=4)+
  scale_fill_manual("Sites",breaks=site_area_10$a,values=custom_colors)+
  geom_smooth(data=site_area_10,pch=21,aes(x=pred_rich2,y=log(richness_mean)),method="lm")+
  
  xlab("Predicted site-level richness")+
  ylab("")+
  theme(legend.position ="right", 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  annotate("text",x=10,y=9.5,label=expression(italic(R^2)*"=42%**"),size=6)




# to see how many plots are there in a site

# to get the richness at the site level and all the cores were selected

set.seed=(6900)
times=30

a4=sample_data(a1)
a4=unique(a4$siteID)

richness <- matrix(nrow=length(a4),ncol=2)
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, siteID==a4[i])
  
  df=data.frame(t(colSums(otu_table(data_sub))))
  # Step 2: Estimate Chao1 diversity
  abundance_matrix <- as.matrix(df[,1:ncol(df)])
  
  chao1_estimate <- estimateR(abundance_matrix, method = "chao1")
  
  richness[i,1]=chao1_estimate[1,]
  richness[i,2]=chao1_estimate[2,]
  }
  # if we convert all the value into 0-1

rar_richness=cbind(a4,richness)%>%data.frame()

names(rar_richness)=c("a","obrich","chao1")


a4=sample_data(a1)
a4=unique(a4$siteID)

richness <- matrix(nrow=length(a4),ncol=2)
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, siteID==a4[i])
  
  df=data.frame(t(colSums(otu_table(data_sub))))
  df[df>0]=1
  # Step 2: Estimate Chao1 diversity
  abundance_matrix <- as.matrix(df[,1:ncol(df)])
  
  chao1_estimate <- estimateR(abundance_matrix, method = "chao1")
  
  richness[i,1]=chao1_estimate[1,]
  richness[i,2]=chao1_estimate[2,]
}
  
 
  
  
  
  
  dmm=data.frame(t(dmm))
  rownames(dmm)="try"
  
  tax_table1=tax_table(dmm)
  
  samp=matrix(ncol=5,nrow=1,1:5)%>%data.frame()
  
  sampdf=sample_data(samp)
  
  abundance_matrix <- as.matrix(dmm[,1:ncol(dmm)])
  
  physeq <- phyloseq(otu_table(abundance_matrix, taxa_are_rows = TRUE),
                     tax_table1,
                     sampdf)
    

    richness[[i]]=estimate_richness(data_sub,measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
  }



richness_site_neon_rarify=data.frame(nrow=length(a4),ncol=3)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_site_neon[i,1]=richness[[i]] 
    richness_site_neon[i,2]=a4[i]
    richness_site_neon[i,3]=0
  }
  
  else
  {
    richness_site_neon[i,1]=mean(richness[[i]])
    richness_site_neon[i,2]=a4[i]
    richness_site_neon[i,3]=sd(richness[[i]])
  }
}

names(richness_site_neon)=c("richness_mean","a","richness_sd")










# Nested method for the NEON site,which means that the prior sites were included when adding the new ones




a2= sample_data(a1)# the unique of the plotID, we have 476 plots
a2=unique(a2$siteID)


species <- vector("list", length(a2))

for (i in 1:length(a2))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(a1, siteID==a2[i])
  dim1 <- dim(sample_data(neon_sub )) # the number of samples in one site
  
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) # 
      { 
      
        sample_seq <- shuffle(c(1:dim1[1]))
        # take out otu_tab for each site
        otu_tab <- otu_table(neon_sub)
        otu_tab <- matrix(otu_tab, nrow = dim1[1], ncol=67769,byrow  = TRUE)#
        species[1] <- sum(otu_tab[sample_seq[1],] > 0)
        
        for (j in 2:(dim1[1]))
        { 
          cat('\r',paste(paste0(rep("*", round(j/ 1, 0)), collapse = ''), j, collapse = ''))# inf
          # take out samples as the sequence in sample_seq
          temp <- colSums(otu_tab[c(sample_seq[1:j]),])
          # count species
          species[j] <- sum(temp > 0)
        }
      }
      
      return(c(temp[2,1], temp[1,1],temp[2,4],temp[1,4]))
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




# create a plot for the HARV_037 plot for an example
#
save(d_HARV_037,file="d_HARV_037.RData")

a=ggplot()+
  geom_point(data=d_HARV_037,aes(x=log(area),y=log(richness)),size=3)+
  geom_segment(data=d_HARV_037,aes(x=log(area),y=log(richness)-log(sdd)/sqrt(30),xend=log(area),yend=log(richness)+log(sdd)/sqrt(30)))+
  xlab("log(Area)")+
  ylab("Log(Richness)")+
  theme(legend.position = c(0.5,0.8), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))

sar_dob_permutation=subset(sar_dob_permutation,plotid!="CO2")

b=ggplot(sar_dob_permutation)+geom_point(data=sar_dob_permutation,aes(x=log(area),y=log(richness),color=plotid))+
  geom_smooth(data=sar_dob_permutation,aes(x=log(area),y=log(richness),color=plotid),se=FALSE,method="lm")+
  guides(color="none")+
  xlab("Log(Area)")+
  ylab("Log(Richness)")+
  theme(legend.position = c(0.5,0.8), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))

###

ggplot(sar_dob_permutation)+geom_point(data=sar_dob_permutation,aes(x=log(area),y=log(richness),color=plotid))+
  geom_smooth(data=sar_dob_permutation,aes(x=log(area),y=log(richness),color=plotid),se=FALSE,method="lm")+
  guides(color="none")+
  xlab("Log(Area)")+
  ylab("Log(Richness)")+
  theme(legend.position = c(0.5,0.8), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))

  
  

