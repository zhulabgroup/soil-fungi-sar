##use the standardized approach to quantify the z values

dob=dob_permutation
## get the sample at the 10 by 10 scale for the dob site

set.seed=(2202)
times=30
a4=sample_data(dob)

a4=unique(a4$subplotID10)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, subplotID10==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  if(dim1>1)
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    richness[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq","dplyr")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:1) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], 1)] <- TRUE# select 3 cores 
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else if (dim1<1)
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

richness_subplot10_dob%>%mutate(plotid=substr(richness_subplot10_dob$ncol,1,3))%>%dplyr::rename(richness=nrow)%>%
  group_by(plotid)%>%summarise(mean_value=mean(richness,na.rm = TRUE),sd_value = sd(richness,na.rm = TRUE))->richness_subplot10_dob

richness_subplot10_dob_standar=richness_subplot10_dob%>%mutate(area=rep(100,43))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")

save(richness_subplot10_dob_standar,file="richness_subplot10_dob_standar.RData")

#at the 20 m scale, we have already quantified this

set.seed=(1203)
times=30
a4=sample_data(dob)

a4=unique(a4$subplotID20)

richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob, subplotID20==a4[i])
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

richness_subplot20_dob=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  
  if(is.null(dim(richness[[i]])[1]))
  {
    richness_subplot20_dob[i,1]=richness[[i]] 
    richness_subplot20_dob[i,2]=richness[[i]] 
  }
  
  else
  {
    richness_subplot20_dob[i,1]=mean(richness[[i]])
    richness_subplot20_dob[i,2]=sd(richness[[i]])
  }
}


richness_subplot20_dob_standar=richness_subplot20_dob%>%mutate(plotid=substr(a4,1,3))%>%group_by(plotid)%>%summarise(mean_value=mean(nrow,na.rm = TRUE),sd_value = sd(nrow,na.rm = TRUE))%>%mutate(area=rep(400,length.out=n))

save(richness_subplot20_dob_standar,file="richness_subplot20_dob_standar.RData")


##at the 30 by 30 scale

richness_plot_permu=list()
for (p in 1:30)
{
  cat("\r", paste(paste0(rep("*", round(p / 1, 0)), collapse = ""), p, collapse = "")) # informs the processing
  
  richness_plot_dob=numeric()
  for(k in 1:length(plot_ID))
  {
    
    core_number=numeric()
    for(i in 1:dim(ori)[1])
    {
      d=sample_data(dob_permutation)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
      core_number[i]=subset(d,X1>=ori$x[i]&X1<=ori$x[i]+30 &X2<=ori$y[i]+30&X2>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
    }
    core_number_8=core_number>8
    core_number_8[!core_number_8]=NA
    richness=numeric()
    for(j in c(which(!is.na(core_number_8))))# only select the locations that over 9 cores can be sampled/81 cores
    {
      # informs the processing
      sub_samp=subset(d,X1>=ori$x[j]&X1<=ori$x[j]+30 &X2<=ori$y[j]+30&X2>=ori$y[j])
      sub_samp=sample(rownames(sub_samp),9)# select only nine cores
      kk=subset_samples(dob_permutation,sample_names(dob_permutation)%in%sub_samp)
      richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
    } 
    
    richness_plot_dob[k]=richness%>%mean(na.rm = TRUE)
  }
  richness_plot_permu[[p]]=richness_plot_dob
}

# get the mean for the 30 permutations

mean_permu=matrix(ncol=length(plot_ID),nrow=30)
for(p in 1:30)
{
  mean_permu[p,] =richness_plot_permu[[p]]
}

richness_subplot30_dob_standar=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar,file="richness_subplot30_dob_standar.RData")



#at the 40 by 40 scale
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
    
    if(dim1>=2)
    {
      data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
      ot=t(otu_table(data_sub))%>%data.frame()
      n=dim(ot)[2]# the number of "sites" for each plot
      ot1=rowSums(ot)
      out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 16, length.out=16)), se=FALSE)
      richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    }
    
    else{
      richness[[i]]=NA
    }
  }
}


richness_mean_subplot40_dob=matrix(ncol=1,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,]=mean(richness[[i]])
}

richness_subplot40_dob_standar=cbind(plotid=a4,richness=richness_mean_subplot40_dob,rep(1600,43))%>%data.frame()%>%rename_all(~paste0(c("plotid","richness","area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")

save(richness_subplot40_dob_standar,file="richness_subplot40_dob_standar.RData")

