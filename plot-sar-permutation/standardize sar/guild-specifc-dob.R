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

### for the AM

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

########################## for the soil saprophytic fungi#######################################

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

####################################for the litter saprophytic fungi################################################

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

## for the plant pathogens


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

### for the wood sap


## for different guilds
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

###

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

###


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

### also at the 40 scale



