## for the guild specific for the dob site
# the coordinates should be in the format of numeric


tax_table(dob)%>%data.frame()%>%dplyr::select(ta2)%>%table()%>%as.data.frame()%>%arrange(Freq)
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

# to select eight fungal guilds

data_EM <- subset_taxa(dob, ta2 == "ectomycorrhizal")
data_AM <- subset_taxa(dob, ta2 == "arbuscular_mycorrhizal")
data_soilsap <- subset_taxa(dob, ta2 == "soil_saprotroph")
data_littersap <- subset_taxa(dob, ta2 == "litter_saprotroph")
data_plapat <- subset_taxa(dob, ta2 == "plant_pathogen")
data_woodsap <- subset_taxa(dob, ta2 == "wood_saprotroph")
data_para <- subset_taxa(dob, ta2%in%c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite"))
data_epiphy <- subset_taxa(dob, ta2 == "epiphyte")

### at the 30 m scale for the dob sites
set.seed(123)
pp=sample_data(data_EM)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_EM)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

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
      d=sample_data(data_EM)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_EM,sample_names(data_EM)%in%sub_samp)
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

richness_subplot30_dob_standar_EM=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,sd_value,area)
save(richness_subplot30_dob_standar_EM,file="richness_subplot30_dob_standar_EM.RData")

########################at the 40 by 40 scale###################################
set.seed=(567)
times=30
a4=sample_data(data_EM)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_EM, plotIDM==a4[i])
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


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_EM=bind_cols(plotid=a4)%>%mutate(richness=richness_mean_subplot40_dob,rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","richness","area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_EM,file="richness_subplot40_dob_standar_EM.RData")




### for the AM

set.seed=(989)
times=30
a4=sample_data(data_AM)
a4=unique(a4$subplotID10)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_AM, subplotID10==a4[i])
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
        tEMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tEMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}



richness_subplot10_dob_AM=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_dob_AM[i,1]=a4[i]
    richness_subplot10_dob_AM[i,2]=richness[[i]] 
    
  }
  else
  {
    
    richness_subplot10_dob_AM[i,1]=a4[i]
    richness_subplot10_dob_AM[i,2]=mean(richness[[i]])
  }
}


# get the mean value for each subplot

richness_subplot10_dob_AM%>%mutate(plotid=substr(richness_subplot10_dob_AM$nrow,1,3))%>%dplyr::rename(richness=ncol)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm=TRUE),sd_value = sd(richness,na.rm=TRUE))->richness_subplot10_dob_standar_AM

save(richness_subplot10_dob_standar_AM,file="richness_subplot10_dob_standar_AM.RData")

############

pp=sample_data(data_AM)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_AM)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

set.seed(456)
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
      d=sample_data(data_AM)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_AM,sample_names(data_AM)%in%sub_samp)
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

richness_subplot30_dob_standar_AM=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,sd_value,area)


save(richness_subplot30_dob_standar_AM,file="richness_subplot30_dob_standar_AM.RData")

####################at the 40 by 40 m scale###################

set.seed=(567)
times=30
a4=sample_data(data_AM)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_AM, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_AM=richness_mean_subplot40_dob%>%bind_cols(plotid=a4)%>%mutate(area=rep(1600,lengh.out=n))%>%rename_all(~paste0(c("mean_value","sd_value","plotid" ,"area")))%>%select(plotid, mean_value, sd_value,  area)

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_AM,file="richness_subplot40_dob_standar_AM.RData")
####




########################## for the soil saprophytic fungi#######################################



set.seed=(989)
times=30
a4=sample_data(data_soilsap)
a4=unique(a4$subplotID10)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_soilsap, subplotID10==a4[i])
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
        tsoilsapp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tsoilsapp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}



richness_subplot10_dob_soilsap=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_dob_soilsap[i,1]=a4[i]
    richness_subplot10_dob_soilsap[i,2]=richness[[i]] 
    
  }
  else
  {
    
    richness_subplot10_dob_soilsap[i,1]=a4[i]
    richness_subplot10_dob_soilsap[i,2]=mean(richness[[i]])
  }
}


# get the mean value for each subplot

richness_subplot10_dob_soilsap%>%mutate(plotid=substr(richness_subplot10_dob_soilsap$nrow,1,3))%>%dplyr::rename(richness=ncol)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm=TRUE),sd_value = sd(richness,na.rm=TRUE))->richness_subplot10_dob_soilsap

save(richness_subplot10_dob_soilsap,file="richness_subplot10_dob_soilsap.RData")# should be named



##############################

pp=sample_data(data_soilsap)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_soilsap)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

set.seed(456)
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
      d=sample_data(data_soilsap)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_soilsap,sample_names(data_soilsap)%in%sub_samp)
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

richness_subplot30_dob_standar_soilsap=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar_soilsap,file="richness_subplot30_dob_standar_soilsap.RData")



set.seed=(567)
times=30
a4=sample_data(data_soilsap)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_soilsap, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_soilsap=bind_cols(plotid=a4,richness_mean_subplot40_dob)%>%mutate(area=rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value", "area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_soilsap,file="richness_subplot40_dob_standar_soilsap.RData")


####################################for the litter saprophytic fungi################################################


set.seed=(989)
times=30
a4=sample_data(data_littersap)
a4=unique(a4$subplotID10)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, subplotID10==a4[i])
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
        tEMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tEMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}



richness_subplot10_dob_littersap=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_dob_littersap[i,1]=a4[i]
    richness_subplot10_dob_littersap[i,2]=richness[[i]] 
    
  }
  else
  {
    
    richness_subplot10_dob_littersap[i,1]=a4[i]
    richness_subplot10_dob_littersap[i,2]=mean(richness[[i]])
  }
}


# get the mean value for each subplot

richness_subplot10_dob_littersap%>%mutate(plotid=substr(richness_subplot10_dob_littersap$nrow,1,3))%>%dplyr::rename(richness=ncol)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm=TRUE),sd_value = sd(richness,na.rm=TRUE))->richness_subplot10_dob_standar_littersap

save(richness_subplot10_dob_standar_littersap,file="richness_subplot10_dob_standar_littersap.RData")


set.seed(456)
pp=sample_data(data_littersap)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_littersap)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)


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
      d=sample_data(data_littersap)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_littersap,sample_names(data_littersap)%in%sub_samp)
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

richness_subplot30_dob_standar_littersap=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar_littersap,file="richness_subplot30_dob_standar_littersap.RData")


### at the 40 by 40m scale

set.seed=(567)
times=30
a4=sample_data(data_littersap)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_littersap, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_littersap=bind_cols(plotid=a4,richness_mean_subplot40_dob)%>%mutate(area=rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value", "area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_littersap,file="richness_subplot40_dob_standar_littersap.RData")


## for the plant pathogens

##############
set.seed=(989)
times=30
a4=sample_data(data_plapat)
a4=unique(a4$subplotID10)
richness <- vector("list", length(a4))
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_plapat, subplotID10==a4[i])
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
        tEMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tEMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}

richness_subplot10_dob_plapat=data.frame(nrow=length(a4),ncol=2)# the mean indicates the mean value for the cores within the same subplot
for (i in 1:length(a4))
{
  ak=dim(richness[[i]])[1]
  
  if(is.null(ak))
  {
    richness_subplot10_dob_plapat[i,1]=a4[i]
    richness_subplot10_dob_plapat[i,2]=richness[[i]] 
    
  }
  else
  {
    
    richness_subplot10_dob_plapat[i,1]=a4[i]
    richness_subplot10_dob_plapat[i,2]=mean(richness[[i]])
  }
}


# get the mean value for each subplot

richness_subplot10_dob_plapat%>%mutate(plotid=substr(richness_subplot10_dob_plapat$nrow,1,3))%>%dplyr::rename(richness=ncol)%>%dplyr::group_by(plotid)%>%dplyr::summarise(mean_value=mean(richness,na.rm=TRUE),sd_value = sd(richness,na.rm=TRUE))->richness_subplot10_dob_standar_plapat

save(richness_subplot10_dob_standar_plapat,file="richness_subplot10_dob_standar_plapat.RData")

###################plapat at the 30 m scale###############
pp=sample_data(data_plapat)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_plapat)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

set.seed(456)
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
      d=sample_data(data_plapat)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_plapat,sample_names(data_plapat)%in%sub_samp)
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

richness_subplot30_dob_standar_plapat=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar_plapat,file="richness_subplot30_dob_standar_plapat.RData")

### at the 40 by 40 m scale
set.seed=(567)
times=30
a4=sample_data(data_plapat)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_plapat, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_plapat=bind_cols(plotid=a4,richness_mean_subplot40_dob)%>%mutate(area=rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value", "area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_plapat,file="richness_subplot40_dob_standar_plapat.RData")


### for the wood sap

set.seed(456)
pp=sample_data(data_woodsap)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_woodsap)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)


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
      d=sample_data(data_woodsap)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_woodsap,sample_names(data_woodsap)%in%sub_samp)
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

richness_subplot30_dob_standar_woodsap=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar_woodsap,file="richness_subplot30_dob_standar_woodsap.RData")

###at the 40 by 40 scale##########

set.seed=(567)
times=30
a4=sample_data(data_woodsap)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_woodsap, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_woodsap=bind_cols(plotid=a4,richness_mean_subplot40_dob)%>%mutate(area=rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value", "area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_woodsap,file="richness_subplot40_dob_standar_woodsap.RData")


## for different guilds

pp=sample_data(data_para)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_para)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

set.seed(456)
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
      d=sample_data(data_para)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_para,sample_names(data_para)%in%sub_samp)
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

richness_subplot30_dob_standar_para=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar_para,file="richness_subplot30_dob_standar_para.RData")

###at the 40 by 40m scale


set.seed=(567)
times=30
a4=sample_data(data_para)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_para, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_para=bind_cols(plotid=a4,richness_mean_subplot40_dob)%>%mutate(area=rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value", "area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_para,file="richness_subplot40_dob_standar_para.RData")


## for different guilds

pp=sample_data(data_epiphy)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

plot_ID=sample_data(data_epiphy)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)

set.seed=(123)

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
      d=sample_data(data_epiphy)%>%data.frame()%>%select(X1,X2,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
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
      kk=subset_samples(data_epiphy,sample_names(data_epiphy)%in%sub_samp)
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

richness_subplot30_dob_standar_epiphy=mean_permu%>%data.frame()%>%t()%>%data.frame()%>%rowwise()%>% summarise(mean_value = mean(c_across(everything()), na.rm = TRUE),sd_value = sd(c_across(everything()), na.rm = TRUE))%>%mutate(plotID=plot_ID,area=rep(900,length.out=n))%>%select( plotID,mean_value,area)

save(richness_subplot30_dob_standar_epiphy,file="richness_subplot30_dob_standar_epiphy.RData")

### also at the 40 by 40 m scale
set.seed=(567)
times=30
a4=sample_data(data_epiphy)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_epiphy, plotIDM==a4[i])
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
        tAMp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(tAMp)["TRUE"] > 0)
        return(species[j])
      }
    }
    stopCluster(cl)
  }
  else{
    richness[[i]]=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  }
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=43)
for (i in 1:43)
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])# also we can have the sd
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])#
  
}

richness_subplot40_dob_standar_epiphy=bind_cols(plotid=a4,richness_mean_subplot40_dob)%>%mutate(area=rep(1600,lengh.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value", "area")))

setwd("/Users/luowenqi/soil-sar/plot-sar-permutation/standardize sar")
save(richness_subplot40_dob_standar_epiphy,file="richness_subplot40_dob_standar_epiphy.RData")



