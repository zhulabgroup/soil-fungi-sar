
# get the dob data

library(stringr)

dob= subset_samples(rare_all,Project=="DoB")# 

########
subplot_ID=sample_data(dob)%>%data.frame()
sample_ID=substr(subplot_ID$geneticSampleID,1,7)

assign_ID=list()
for (i in 1:908)
{
  sub_sample=sample_ID[[i]]
  
  if(str_detect(sub_sample,"00"))
  {
    assign_ID[[i]]=c(0,0)
  }
  else if(str_detect(sub_sample,"A05|A5"))
  {
    assign_ID[[i]]=c(0,5)
    
  }
  
  else if(str_detect(sub_sample,"A10"))
  {
    assign_ID[[i]]=c(0,10)
    
  }
  
  else if(str_detect(sub_sample,"A20"))
  {
    assign_ID[[i]]=c(0,20)
  }
  else if(str_detect(sub_sample,"A40"))
  {
    assign_ID[[i]]=c(0,40)
  }
  else if(str_detect(sub_sample,"B05|B5"))
  {
    assign_ID[[i]]=c(5,5)
    
  }
  else if(str_detect(sub_sample,"B10"))
  {
    assign_ID[[i]]=c(10,10)
  }
  else if(str_detect(sub_sample,"B20"))
  {
    assign_ID[[i]]=c(20,20)
  }
  
  else if(str_detect(sub_sample,"B40"))
  {
    assign_ID[[i]]=c(40,40)
  }
  
  else if(str_detect(sub_sample,"C05|C5")&!str_detect(sub_sample,"BC5"))# done
  {
    assign_ID[[i]]=c(5,0)
  }
  
  else if(str_detect(sub_sample,"C10"))
  {
    assign_ID[[i]]=c(10,0)
  }
  
  else if(str_detect(sub_sample,"C20"))
  {
    assign_ID[[i]]=c(20,0)
  }
  else if(str_detect(sub_sample,"C40"))
  {
    assign_ID[[i]]=c(40,0)
  }
  else if(str_detect(sub_sample,"BC5.C5."))# done
    {
    assign_ID[[i]]=c(5,5)
  }
  else
  {
    assign_ID[[i]]=c(15,15)
  }
  }

# rbind all the coordinates 

coord_dob=matrix(nrow=908,ncol=2)
   for (i in 1:908)  {
     
     coord_dob[i,1]=assign_ID[[i]][1]
     coord_dob[i,2]=assign_ID[[i]][2]
   }            
                 
# to add a ID to the data set for the dob sites

coord_dob=data.frame(coord_dob)
row.names(coord_dob)=row.names(sample_data(dob))
coord_dob=sample_data(coord_dob)

dob=merge_phyloseq(dob,coord_dob)

# divide the plot into subplots at the 5 by 5 m spatial scales

subplot_ID=sample_data(dob)%>%data.frame()
subplot_ID=subplot_ID[,c("X1","X2")]

# Define the dimensions of the plot area
plot_width <- 40
plot_height <- 40

# at the 5 by 5 m spatial scales

# Define the dimensions of each grid cell
grid_cell_size <- 5

# Create grid boundaries
x_breaks <- seq(0, plot_width, by = grid_cell_size)
y_breaks <- seq(0, plot_height, by = grid_cell_size)

# Assign a ID to each point based on their coordinates

subplot_ID$X1 <- cut(subplot_ID$X1, breaks = x_breaks, labels = FALSE,include.lowest = TRUE)
subplot_ID$X2 <- cut(subplot_ID$X2, breaks = y_breaks, labels = FALSE,include.lowest = TRUE)

# Display the result
head(subplot_ID)

names(subplot_ID)=c("gx","gy")


# add subplot ID to the phyloseq data

plotIDM=sample_data(dob)
plotIDM=plotIDM$plotIDM

subplotID5=paste(plotIDM,"*",paste(subplot_ID$gx,"*",subplot_ID$gy))%>%data.frame()

# add the ID to the phyloseq

names(subplotID5)="subplotID5"
row.names(subplotID5)=row.names(sample_data(dob))
subplotID5=sample_data(subplotID5)
dob=merge_phyloseq(dob,subplotID5)

## get the samples at the 5 by 5 scale for the dob site

set.seed=(2201)
times=30
a4=sample_data(dob)

a4=unique(a4$subplotID5)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, subplotID5==a4[i])
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
# get the number of species within each 5 by 5 subplots

richness_subplot5_dob=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot5_dob[i,1]=richness[[i]] 
    richness_subplot5_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot5_dob[i,1]=mean(richness[[i]])
    richness_subplot5_dob[i,2]=a4[i]
  }
}

# get the mean richness of the subplots

richness_subplot5_dob=mutate(richness_subplot5_dob,plotid=substr(richness_subplot5_dob$ncol,1,3))

names(richness_subplot5_dob)[1]="richness"

richness_sd_subplot5_dob=aggregate(richness~plotid,data=richness_subplot5_dob,FUN=sd)

richness_mean_subplot5_dob=aggregate(richness~plotid,data=richness_subplot5_dob,FUN=mean)

# for the 10 by 10 subplots
# divide the plot into several subplots at the 10 by 10 scal

subplot_ID=sample_data(dob)%>%data.frame()
subplot_ID=subplot_ID[,c("X1","X2")]

# Define the dimensions of the plot area
plot_width <- 40
plot_height <- 40
grid_cell_size <- 10

# Create grid boundaries
x_breaks <- seq(0, plot_width, by = grid_cell_size)
y_breaks <- seq(0, plot_height, by = grid_cell_size)

# Assign a ID to each point based on their coordinates

subplot_ID$X1 <- cut(subplot_ID$X1, breaks = x_breaks, labels = FALSE,include.lowest = TRUE)
subplot_ID$X2 <- cut(subplot_ID$X2, breaks = y_breaks, labels = FALSE,include.lowest = TRUE)
names(subplot_ID)=c("gx","gy")

# add the plot ID to the subplot ID

plotIDM=sample_data(dob)
plotIDM=plotIDM$plotIDM

subplotID10=paste(plotIDM,"*",paste(subplot_ID$gx,"*",subplot_ID$gy))%>%data.frame()

# add the ID to the phyloseq
names(subplotID10)="subplotID10"
row.names(subplotID10)=row.names(sample_data(dob))

subplotID10=sample_data(subplotID10)
dob=merge_phyloseq(dob,subplotID10)

###



## get the sample at the 10 by 10 scale for the dob site

set.seed=(2202)
times=30
a4=sample_data(dob)

a4=unique(a4$subplotID10new)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, subplotID10new==a4[i])
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

# get the richness at the 10 by 10

richness_subplot10_dob=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_dob[i,1]=richness[[i]] 
    richness_subplot10_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot10_dob[i,1]=mean(richness[[i]])
    richness_subplot10_dob[i,2]=a4[i]
  }
}


richness_subplot10_dob=mutate(richness_subplot10_dob,plotid=substr(richness_subplot10_dob$ncol,1,3))

names(richness_subplot10_dob)[1]="richness"

richness_sd_subplot10_dob=aggregate(richness~plotid,data=richness_subplot10_dob,FUN=sd)
richness_mean_subplot10_dob=aggregate(richness~plotid,data=richness_subplot10_dob,FUN=mean)

# at the 20 by 20 spatial scale

subplot_ID=sample_data(dob)%>%data.frame()
subplot_ID=subplot_ID[,c("X1","X2")]

# Define the dimensions of the plot area
plot_width <- 40
plot_height <- 40

# Define the dimensions of each grid cell
grid_cell_size <- 20

# Create grid boundaries
x_breaks <- seq(0, plot_width, by = grid_cell_size)
y_breaks <- seq(0, plot_height, by = grid_cell_size)

# Assign a ID to each point based on their coordinates

subplot_ID$X1 <- cut(subplot_ID$X1, breaks = x_breaks, labels = FALSE,include.lowest = TRUE)
subplot_ID$X2 <- cut(subplot_ID$X2, breaks = y_breaks, labels = FALSE,include.lowest = TRUE)

# Display the result
head(subplot_ID)

names(subplot_ID)=c("gx","gy")

# add the plot ID to the subplot ID

plotIDM=sample_data(dob)
plotIDM=plotIDM$plotIDM

subplotID20new=paste(plotIDM,"-",paste(subplot_ID$gx,"-",subplot_ID$gy))
subplotID20new=data.frame(subplotID20new)
# add thi ID to the phyloseq
names(subplotID20new)="subplotID20new"
row.names(subplotID20new)=row.names(sample_data(dob))

subplotID20new=sample_data(subplotID20new)
dob=merge_phyloseq(dob,subplotID20new)
## get the sample at the 20 by 20 scale for the dob site

set.seed=(2204)
times=30
a4=sample_data(dob)

a4=unique(a4$subplotID20new)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, subplotID20new==a4[i])
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



# get the richness of the 20 by 20 plot
richness_subplot20_dob=data.frame(nrow=158,ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot20_dob[i,1]=richness[[i]] 
    richness_subplot20_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot20_dob[i,1]=mean(richness[[i]])
    richness_subplot20_dob[i,2]=a4[i]
  }
}


richness_subplot20_dob=mutate(richness_subplot20_dob,plotid=substr(richness_subplot20_dob$ncol,1,3))

names(richness_subplot20_dob)[1]="richness"

richness_sd_subplot20_dob=aggregate(richness~plotid,data=richness_subplot20_dob,FUN=sd)
richness_mean_subplot20_dob=aggregate(richness~plotid,data=richness_subplot20_dob,FUN=mean)

# for the 30 by 30 plot level 
# there four approaches to define this scale

subplot_ID=sample_data(dob)%>%data.frame()
subplot_ID=subplot_ID[,c("X1","X2")]


# the first approach
names(subplot_ID)=c("gx","gy")

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

row.names(four_assign)=row.names(sample_data(dob))
four_assign=sample_data(four_assign)

dob=merge_phyloseq(dob,four_assign)

## get the richness at the 30 by 30 m spatial scales
# 

set.seed=(5566)
times=30
a4=sample_data(dob)

a4=unique(a4$point_assign1)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, point_assign1==a4[i])
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

richness_subplot30A_dob=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30A_dob[i,1]=richness[[i]] 
    richness_subplot30A_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30A_dob[i,1]=mean(richness[[i]])
    richness_subplot30A_dob[i,2]=a4[i]
  }
}

richness_subplot30A_dob=mutate(richness_subplot30A_dob,plotid=substr(richness_subplot30A_dob$ncol,1,3))

names(richness_subplot30A_dob)[1]="richness"

# the mean richness for each 10 x 10 subplot


richness_mean_subplot30A_dob=aggregate(richness~plotid,data=richness_subplot30A_dob,FUN=mean,na.rm=TRUE)
richness_sd_subplo30A_dob=aggregate(richness~plotid,data=richness_subplot30A_dob,FUN=sd,na.rm=TRUE)


# for the second approach

set.seed=(5567)
times=30
a4=sample_data(dob)

a4=unique(a4$point_assign2)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, point_assign2==a4[i])
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

richness_subplot30B_dob=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30B_dob[i,1]=richness[[i]] 
    richness_subplot30B_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30B_dob[i,1]=mean(richness[[i]])
    richness_subplot30B_dob[i,2]=a4[i]
  }
}

richness_subplot30B_dob=mutate(richness_subplot30B_dob,plotid=substr(richness_subplot30B_dob$ncol,1,3))

names(richness_subplot30B_dob)[1]="richness"

# the mean richness for each 10 x 10 subplot


richness_mean_subplot30B_dob=aggregate(richness~plotid,data=richness_subplot30B_dob,FUN=mean,na.rm=TRUE)
richness_sd_subplo30B_dob=aggregate(richness~plotid,data=richness_subplot30B_dob,FUN=sd,na.rm=TRUE)


## for the third method

set.seed=(5568)
times=30
a4=sample_data(dob)

a4=unique(a4$point_assign3)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, point_assign3==a4[i])
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

richness_subplot30C_dob=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30C_dob[i,1]=richness[[i]] 
    richness_subplot30C_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30C_dob[i,1]=mean(richness[[i]])
    richness_subplot30C_dob[i,2]=a4[i]
  }
}

richness_subplot30C_dob=mutate(richness_subplot30C_dob,plotid=substr(richness_subplot30C_dob$ncol,1,3))

names(richness_subplot30C_dob)[1]="richness"

# the mean richness for each 10 x 10 subplot


richness_mean_subplot30C_dob=aggregate(richness~plotid,data=richness_subplot30C_dob,FUN=mean,na.rm=TRUE)
richness_sd_subplo30C_dob=aggregate(richness~plotid,data=richness_subplot30C_dob,FUN=sd,na.rm=TRUE)


####

set.seed=(5568)
times=30
a4=sample_data(dob)

a4=unique(a4$point_assign4)

mm=str_detect(a4," * 1")

mm[!mm]=NA

which(!is.na(mm))

richness <- vector("list", length(which(!is.na(mm))))

for (i in c(1,which(!is.na(mm))))# only select the subplot coded with "-1"
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, point_assign4==a4[i])
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

richness_subplot30D_dob=data.frame(nrow=length(which(!is.na(mm))),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in c(1,which(!is.na(mm))))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot30D_dob[i,1]=richness[[i]] 
    richness_subplot30D_dob[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot30D_dob[i,1]=mean(richness[[i]])
    richness_subplot30D_dob[i,2]=a4[i]
  }
}

richness_subplot30D_dob=mutate(richness_subplot30D_dob,plotid=substr(richness_subplot30D_dob$ncol,1,3))

names(richness_subplot30D_dob)[1]="richness"

# the mean richness for each 10 x 10 subplot


richness_mean_subplot30D_dob=aggregate(richness~plotid,data=richness_subplot30D_dob,FUN=mean,na.rm=TRUE)
richness_sd_subplo30D_dob=aggregate(richness~plotid,data=richness_subplot30D_dob,FUN=sd,na.rm=TRUE)

### rbind all the estimates based on the four different approaches

bind_rows(richness_mean_subplot30A_dob,richness_mean_subplot30B_dob,richness_mean_subplot30C_dob,richness_mean_subplot30D_dob)%>%group_by(plotid)%>%dplyr::summarise(mean_value=mean(mean_value,na.rm = TRUE),sd_value = sd(sd_value,na.rm = TRUE))->richness_mean_subplotABCD_dob

# at the 40 by 40 scale

set.seed=(5679)

times=30
a4=sample_data(dob)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, plotIDM==a4[i])
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


richness_mean_subplot40_dob=matrix(ncol=1,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,]=mean(richness[[i]])
}

richness_mean_subplot40_dob=cbind(plotid=a4,richness=richness_mean_subplot40_dob,rep(1600,43))%>%data.frame()

names(richness_mean_subplot40_dob)=c("plotid","richness","area")

## rbind all the different spatial scales
richness_mean_subplot5_dob=mutate(richness_mean_subplot5_dob,area=rep(25,43))

richness_mean_subplot10_dob=mutate(richness_mean_subplot10_dob,area=rep(100,41))

richness_mean_subplot20_dob=mutate(richness_mean_subplot20_dob,area=rep(400,42))

richness_mean_subplot30_dob=mutate(richness_mean_subplotABCD_dob,area=rep(900,40))

richness_mean_subplot40_dob=mutate(richness_mean_subplot40_dob,area=rep(1600,43))
## 

sar_dob_permutation=rbind(richness_mean_subplot5_dob,richness_mean_subplot10_dob,richness_mean_subplot20_dob,richness_mean_subplot30_dob,richness_mean_subplot40_dob)

sar_dob_permutation$richness=as.numeric(sar_dob_permutation$richness)
sar_dob_permutation$area=as.numeric(sar_dob_permutation$area)

save(sar_dob_permutation,file="sar_dob_permutation.RData")
dob_permutation=dob
save(dob_permutation,file="dob_permutation.RData")

save(sar_dob_permutation,file="sar_dob_permutation.RData")



