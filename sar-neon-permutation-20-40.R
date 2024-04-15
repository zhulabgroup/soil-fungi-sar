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
  
  

