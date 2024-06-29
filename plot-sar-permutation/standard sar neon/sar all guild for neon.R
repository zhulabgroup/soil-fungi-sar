
p=sample_data(rare_all_assign)%>%data.frame()%>%select(gx,gy)

neon_data=sample_data(rare_all_assign)%>%data.frame()
neon_data$gx=NULL
neon_data$gy=NULL
sample_data(rare_all_assign)=sample_data(neon_data)


p$gx=as.numeric(p$gx)
p$gy=as.numeric(p$gy)
row.names(p)=row.names(sample_data(rare_all_assign))
p=sample_data(p)
rare_all_assign=merge_phyloseq(rare_all_assign,p)
#####

pp=sample_data(rare_all_assign)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))
plot_ID=sample_data(rare_all_assign)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

###at the 10 by 10 m scale

set.seed=(1201)
times=30
a4=sample_data(rare_all_assign)
a4=unique(a4$subplotID10)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all_assign, subplotID10==a4[i])
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

# get the richness for each 5 x 5 m subplot for the NEON sites

richness_subplot10_neon=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_neon[i,1]=richness[[i]] 
    richness_subplot10_neon[i,2]=richness[[i]] 
  }
  else
  {
    richness_subplot10_neon[i,1]=mean(richness[[i]])
    richness_subplot10_neon[i,2]=sd(richness[[i]])
  }
}

#####

richness_subplot10_neon_standar=richness_subplot10_neon%>%mutate(plotid=substr(a4,1,8))%>%group_by(plotid)%>%summarise(mean_value=mean(nrow,na.rm=TRUE),sd_value=sd(nrow,na.rm=TRUE),area=rep(100,length.out=n))
save(richness_subplot10_neon_standar,file="richness_subplot10_neon_standar.RData")



## at the 30 by 30 m scale
plot_ID=sample_data(rare_all_assign)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

set.seed(520)
  richness_plot_neon=numeric()
  for(k in 1:length(plot_ID))
  {
  cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
    
    core_number=numeric()
    for(i in 1:dim(ori)[1])
    {
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      d=sample_data(rare_all_assign)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
      core_number[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
    }
    core_number_8=core_number>8
    core_number_8[!core_number_8]=NA
    
    richness=numeric()
    for(j in c(which(!is.na(core_number_8))))# only select the locations that over 9 cores can be sampled/81 cores
    {
      # informs the processing
      sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
      sub_samp=sample(rownames(sub_samp),9)# select only nine cores
      kk=subset_samples(rare_all_assign,sample_names(rare_all_assign)%in%sub_samp)
      richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
    } 
    
    richness_plot_neon[k]=richness%>%mean(na.rm = TRUE)
  }
  
 ###
  
#richness_subplot40_neon_standar=richness_plot_neon%>%bind_cols(a4)%>%data.frame()%>%mutate(area=rep(900,length.out=n))%>%rename_all(~paste0(c("mean_value","plotid","area")))%>%select(plotid,mean_value,area)

#save(richness_subplot40_neon_standar,file="richness_subplot40_neon_standar.RData")


## at the 40 by 40 scale ###
set.seed=(5679)
times=30
a4=sample_data(rare_all_assign)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all_assign, plotIDM==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>=16)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:16) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 16)] <- TRUE
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


richness_mean_subplot40_neon=matrix(ncol=2,nrow=length(a4))
for (i in 1:length(a4))
{
  
  richness_mean_subplot40_neon[i,1]=mean(richness[[i]])
  richness_mean_subplot40_neon[i,2]=sd(richness[[i]])
}


richness_subplot40_neon_standar=cbind(plotid=a4,mean_value=richness_mean_subplot40_neon[,1],sd_value=richness_mean_subplot40_neon[,2],rep(1600,length.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value","area")))

save(richness_subplot40_neon_standar,file="richness_subplot40_neon_standar.RData")
#





### at the 20 

set.seed=(1203)
times=30
a4=sample_data(rare_all_assign)

a4=unique(a4$subplotID20)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all_assign, subplotID20==a4[i])
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
    richness_subplot20_neon[i,2]=richness[[i]] 
  }
  
  else
  {
    richness_subplot20_neon[i,1]=mean(richness[[i]])
    richness_subplot20_neon[i,2]=sd(richness[[i]])
  }
}

richness_subplot20_neon_standar=richness_subplot20_neon%>%mutate(plotid=substr(a4,1,8),area=rep(400,length.out=n))%>%group_by(plotid)%>%summarise(mean_value=mean(nrow,na.rm = TRUE),sd_value = sd(nrow,na.rm = TRUE))%>%mutate(area=rep(400,length.out=n))%>%rename_all(~paste0(c("plotid","mean_value","sd_value","area")))

save(richness_subplot20_neon_standar,file="richness_subplot20_neon_standar.RData")

####




