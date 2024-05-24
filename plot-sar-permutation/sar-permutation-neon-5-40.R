
0# Generate an subplot ID for each subplot at different scales
load("rare_all.RData")

a1= subset_samples(rare_all,Project=="NEON")# the unique plotID, we have 476 plots
data_sub <- sample_data(a1)%>%data.frame()



0. # add the subplot ID to the phyloseq object

subplot=strsplit(data_sub$geneticSampleID,"-")
subplot_matrix=matrix(ncol=6,nrow=dim(data_sub)[1])
for(i in 1:dim(data_sub)[1]){
  subplot_matrix[i,]=subplot[[i]]
}
subplot_matrix%>%data.frame()->subplot_ID

subplot_ID=subplot_ID[,3:4]


names(subplot_ID)=c("gx","gy")

row.names(subplot_ID)=row.names(sample_data(a1))# need to be the same with the sample data
subplot_ID=sample_data(subplot_ID)
a1=merge_phyloseq(a1,subplot_ID)

# Define the dimensions of the plot area
plot_width <- 40
plot_height <- 40

# at the 5 by 5 m spatial scales

# Define the dimensions of each grid cell
grid_cell_size <- 5

# Create grid boundaries
x_breaks <- seq(0, plot_width, by = grid_cell_size)
y_breaks <- seq(0, plot_height, by = grid_cell_size)

# Assign a ID to each soil core based on their coordinates
subplot_ID$gx=as.numeric(subplot_ID$gx)
subplot_ID$gy=as.numeric(subplot_ID$gy)

subplot_ID$gx <- cut(subplot_ID$gx, breaks = x_breaks, labels = FALSE,include.lowest = TRUE)
subplot_ID$gy <- cut(subplot_ID$gy, breaks = y_breaks, labels = FALSE,include.lowest = TRUE)

# Display the result
head(subplot_ID)

###
3# add the subplotID to the phyloseq object

subplotID5=paste(subplot_ID$gx,"-",subplot_ID$gy)%>%data.frame()

plotIDM=sample_data(a1)%>%data.frame()

plotIDM=plotIDM$plotIDM

subplotID5=paste(plotIDM,"*",subplotID5$.)

subplotID5=data.frame(subplotID5)

names(subplotID5)="subplotID5"

row.names(subplotID5)=row.names(sample_data(a1))

subplotID5=sample_data(subplotID5)

a1=merge_phyloseq(a1,subplotID5)

a1=subset_samples(a1,select=-subplotID5)# if we want to delete the column

2# to see how many cores are there within each of the 5 by 5 m subplot

a4=sample_data(a1)
a4=unique(a4$subplotID5)
con_5=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  
  data_sub <- subset_samples(a1, subplotID5==a4[i])
  
  con_5[i]=dim(sample_data(data_sub))[1]# as long as be assigned 
}
# an error occurred for i=1997

# to get the richness at each 5 x 5 m subplots, and if there were more than two cores we just randomly select one of them, permutating 30 times 

set.seed=(1201)
times=30
a4=sample_data(a1)
a4=unique(a4$subplotID5)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplotID5==a4[i])
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

# get the richness at the 5 x 5 m level for the NEON sites

richness_subplot5_neon=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
  {
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
    {
    richness_subplot5_neon[i,1]=richness[[i]] 
    richness_subplot5_neon[i,2]=a4[i]
  }
 
else
  {
  richness_subplot5_neon[i,1]=mean(richness[[i]])
  richness_subplot5_neon[i,2]=a4[i]
}
}
# get the mean value for each subplot

richness_subplot5_neon=mutate(richness_subplot5_neon,plotid=substr(richness_subplot5_neon$ncol,1,8))

names(richness_subplot5_neon)[1]="richness"

# the mean richness for each 5 x 5 subplot

richness_mean_subplot5_neon=aggregate(richness~plotid,data=richness_subplot5_neon,FUN=mean)

richness_sd_subplot5_neon=aggregate(richness~plotid,data=richness_subplot5_neon,FUN=sd)

3. #at the 10 by 10 plot ,we can also do a permutation.

data_sub <- sample_data(a1)%>%data.frame()
subplot=strsplit(data_sub$geneticSampleID,"-")
subplot_matrix=matrix(ncol=6,nrow=dim(data_sub)[1])
for(i in 1:dim(data_sub)[1]){
  subplot_matrix[i,]=subplot[[i]]
}
subplot_matrix%>%data.frame()->subplot_ID
subplot_ID=subplot_ID[,3:4]
names(subplot_ID)=c("gx","gy")


# Define the dimensions of the plot area
plot_width <- 40
plot_height <- 40

# at the 10 by 10 m spatial scales

# Define the dimensions of each grid cell
grid_cell_size <- 10

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

subplotID10=paste(subplot_ID$gx,"-",subplot_ID$gy)%>%data.frame()

plotIDM=sample_data(a1)%>%data.frame()

plotIDM=plotIDM$plotIDM

subplotID10=paste(plotIDM,"*",subplotID10$.)

subplotID10=data.frame(subplotID10)

names(subplotID10)="subplotID10"

row.names(subplotID10)=row.names(sample_data(a1))

subplotID10=sample_data(subplotID10)

a1=merge_phyloseq(a1,subplotID10)

#

set.seed=(3202)
times=30
a4=sample_data(a1)

a4=unique(a4$subplotID10)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplotID10==a4[i])
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

##  get the richness for each subplot at the 10 by 10 spatial scale

richness_subplot10_neon=data.frame(nrow=length(a4),ncol=2)# the mean indiactes the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_neon[i,1]=richness[[i]] 
    richness_subplot10_neon[i,2]=a4[i]
  }
  
  else
  {
    richness_subplot10_neon[i,1]=mean(richness[[i]])
    richness_subplot10_neon[i,2]=a4[i]
  }
}

richness_subplot10_neon=mutate(richness_subplot10_neon,plotid=substr(richness_subplot10_neon$ncol,1,8))

names(richness_subplot10_neon)[1]="richness"

# the mean richness for each 10 x 10 subplot

richness_mean_subplot10_neon=aggregate(richness~plotid,data=richness_subplot10_neon,FUN=mean,na.rm=TRUE)

richness_sd_subplot10_neon=aggregate(richness~plotid,data=richness_subplot10_neon,FUN=sd,na.rm=TRUE)


