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
    assign_ID[[i]]=c(5,0)
    
  }
  else if(str_detect(sub_sample,"B10"))
  {
    assign_ID[[i]]=c(10,0)
  }
  else if(str_detect(sub_sample,"B20"))
  {
    assign_ID[[i]]=c(20,0)
  }
  
  else if(str_detect(sub_sample,"B40"))
  {
    assign_ID[[i]]=c(40,0)
  }
  
  else if(str_detect(sub_sample,"C05|C5"))
  {
    assign_ID[[i]]=c(5,5)
  }
  
  else if(str_detect(sub_sample,"C10"))
  {
    assign_ID[[i]]=c(10,10)
  }
  
  else if(str_detect(sub_sample,"C20"))
  {
    assign_ID[[i]]=c(20,20)
  }
  else if(str_detect(sub_sample,"C40"))
  {
    assign_ID[[i]]=c(40,40)
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

# devide the plot into subplots at the 5 by 5 m spatial scales

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

sample_data(dob)

# add the plot ID to the subplot ID

plotIDM=sample_data(dob)
plotIDM=plotIDM$plotIDM

subplotID5=paste(plotIDM,"*",paste(subplot_ID$gx,"*",subplot_ID$gy))%>%data.frame()

# add thi ID to the phyloseq
names(subplotID5)="subplotID5"
row.names(subplotID5)=row.names(sample_data(dob))

subplotID5=sample_data(subplotID5)
dob=merge_phyloseq(dob,subplotID5)
## get the sample at the 5 by 5 scale for the dob site

set.seed=(2201)
times=30
a4=sample_data(dob)

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