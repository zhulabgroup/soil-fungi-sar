
## use the extrapolation approach to determine the 40 m scale richness for the neon site

##need to confirm that the gx and gy are in the format of numeric

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
# replace the initial object
#####
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))
#

# for different guilds
# to add a new column called guild to the intial data
# these several groups were grouped into one broad group of "para"
taxa_df <- as.data.frame(tax_table(rare_all_assign))
guild=tax_table(rare_all_assign)%>%data.frame()%>%select(ta2)
term=c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite")

guild%>%mutate(ta2=ifelse(ta2%in%term,"para",ta2))%>%rename(guild=ta2)->guild
taxa_df$guild=guild
new_tax_table <- as.matrix(taxa_df)
tax_table(rare_all_assign) <- new_tax_table

guild_select=c("ectomycorrhizal","arbuscular_mycorrhizal","soil_saprotroph","litter_saprotroph","plant_pathogen","wood_saprotroph",
               "para","epiphyte")

## estimating fungal diversity at the 10 m by 10 m scale for the NEON plots

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

# get the diversity for each 10 x 10 m subplot for the NEON sites

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



richness_subplot10_neon_standar=richness_subplot10_neon%>%mutate(plotid=substr(a4,1,8))%>%group_by(plotid)%>%summarise(mean_value=mean(nrow,na.rm=TRUE),sd_value=sd(nrow,na.rm=TRUE),area=rep(100,length.out=n))
save(richness_subplot10_neon_standar,file="richness_subplot10_neon_standar.RData")




## at the 20 m by 20 m scale for all the guilds

set.seed=(1203)
times=30
a4=sample_data(rare_all_assign)
a4=unique(a4$subplotID20)
richness <- vector("list", length(a4))

for (i in 1:length(a4))# note that the i here is large
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all_assign, subplotID20==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  tryCatch({if(dim1>=4)
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
    else if (dim1>=2&dim1<4)
    {
      sub_samp=rownames(sample_data(data_sub))
      kk=subset_samples(rare_all_assign,sample_names(rare_all_assign)%in%sub_samp)
      data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
      ot=t(otu_table(data_sub))%>%data.frame()
      n=dim(ot)[2]# the number of "sites" for each plot
      ot1=rowSums(ot)
      out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 4, length.out=4)), se=FALSE)
      richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    }
    else{
      
      richness[[i]]=NA
    }
  },
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))}
  )
}


richness_mean_subplot20_neon=matrix(ncol=2,nrow=length(a4))
for (i in 1:length(a4))
{
  richness_mean_subplot20_neon[i,1]=mean(richness[[i]])
  richness_mean_subplot20_neon[i,2]=sd(richness[[i]])
}

richness_mean_subplot20_neon%>%data.frame()%>%bind_cols(a4)%>%mutate(plotid=substr(a4,1,8))%>%
  group_by(plotid)%>%summarise(mean_value=mean(X1,na.rm=TRUE),sd_value=sd(X1,na.rm=TRUE))%>%
  mutate(area=rep(400,n()),guild=rep("all",n()))->data_neon_20_all

### for individual guilds at the 20 by 20 m scale for the neon sites

data=list()
for (k in 1:8){
  
  data_select=subset_taxa(rare_all_assign,guild==guild_select[k])
  
  set.seed=(1203)
  times=30
  a4=sample_data(data_select)
  a4=unique(a4$subplotID20)
  richness <- vector("list", length(a4))
  
  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, subplotID20==a4[i])
    dim1=dim(sample_data(data_sub))[1]
    
    tryCatch({if(dim1>=4)
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
      else if (dim1>=2&dim1<4)
      {
        sub_samp=rownames(sample_data(data_sub))
        kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
        data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
        ot=t(otu_table(data_sub))%>%data.frame()
        n=dim(ot)[2]# the number of "sites" for each plot
        ot1=rowSums(ot)
        out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 4, length.out=4)), se=FALSE)
        richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
      }
      else{
        richness[[i]]=NA
      }
    },
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", i, ":", e))
    },
    finally ={
      print(paste("Processed element", i))}
    )
  }
  
  richness_mean_subplot20_neon=matrix(ncol=2,nrow=length(a4))
  for (i in 1:length(a4))
  {
    richness_mean_subplot20_neon[i,1]=mean(richness[[i]])
    richness_mean_subplot20_neon[i,2]=sd(richness[[i]])
  }
  
  data[[k]]=richness_mean_subplot20_neon%>%data.frame()%>%bind_cols(a4)%>%mutate(plotid=substr(a4,1,8))%>%
    group_by(plotid)%>%summarise(mean_value=mean(X1,na.rm=TRUE),sd_value=sd(X1,na.rm=TRUE))%>%
    mutate(area=rep(400,n()),guild=rep(guild_select[k],n()))
}

data_neon_20_all%>%bind_rows(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],data[[6]],data[[7]],data[[8]])->data_neon_20m_scale

save(data_neon_20m_scale,file="data_neon_20m_scale.RData")



###whole-community fungal diversity at the 30m by 30 m scale

plot_ID=sample_data(rare_all_assign)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)
set.seed(520)
richness_plot_neon=matrix(ncol=3,nrow=length(plot_ID))
for(k in 1:length(plot_ID))
{
  cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
  
  tryCatch({core_number=numeric()
  for(i in 1:dim(ori)[1])
  {
    #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    d=sample_data(rare_all_assign)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==plot_ID[k])
    
    core_number[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
  }
  core_number_8=core_number>8
  core_number_8[!core_number_8]=NA
  
  richness=numeric()
  for(j in 1:length(core_number))# only select the locations that over 9 cores can be sampled/81 cores
  {
    cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
    tryCatch({if (j%in%c(which(!is.na(core_number_8)))){
      # informs the processing
      sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
      sub_samp=sample(rownames(sub_samp),9)# select only nine cores
      kk=subset_samples(rare_all_assign,sample_names(rare_all_assign)%in%sub_samp)
      richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
    } 
      else{
        sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
        if(dim(sub_samp)[1]<2)
        {
          richness[j]=NA 
        }
        else
        {
          sub_samp=rownames(sub_samp)# get the sample names
          kk=subset_samples(rare_all_assign,sample_names(rare_all_assign)%in%sub_samp)#based on this to extrapolate the richness
          data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
          ot=t(otu_table(data_sub))%>%data.frame()
          n=dim(ot)[2]# the number of "sites" for each plot
          ot1=rowSums(ot)
          out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 9, length.out=9)), se=FALSE)
          richness[j]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
        }
      }
    },
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", j, ":", e))
    },
    finally ={
      print(paste("Processed element", j))}
    )
  }
  
  richness_plot_neon[k,1]=plot_ID[k]
  richness_plot_neon[k,2]=richness%>%mean(na.rm = TRUE)
  richness_plot_neon[k,3]=richness%>%sd(na.rm = TRUE)},
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", k, ":", e))
  },
  finally ={
    print(paste("plot finished", k))}
  )
}


# data for the 30 by 30 m scale was computed on the great lakes
# used to search for a 30 m^2 square in the plot
# for individual guilds
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

data=list()
for (m in 1:8)
{
  data_select=subset_taxa(rare_all_assign,guild==guild_select[m])
  plot_ID=sample_data(data_select)%>%data.frame()
  plot_ID=unique(plot_ID$plotIDM)
  set.seed(520)
  richness_plot_neon=matrix(ncol=3,nrow=length(plot_ID))
  for(k in 1:length(plot_ID))
  {
    cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
    
    tryCatch({core_number=numeric()
    for(i in 1:dim(ori)[1])
    {
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      d=sample_data(data_select)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
      core_number[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
    }
    core_number_8=core_number>8
    core_number_8[!core_number_8]=NA
    
    richness=numeric()
    for(j in 1:length(core_number))# only select the locations that over 9 cores can be sampled/81 cores
    {
      cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
      tryCatch({if (j%in%c(which(!is.na(core_number_8)))){
        # informs the processing
        sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
        sub_samp=sample(rownames(sub_samp),9)# select only nine cores
        kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
        richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
      } 
        else{
          sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
          if(dim(sub_samp)[1]<2)
          {
            richness[j]=NA 
          }
          else
          {
            sub_samp=rownames(sub_samp)# get the sample names
            kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)#based on this to extrapolate the richness
            data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
            ot=t(otu_table(data_sub))%>%data.frame()
            n=dim(ot)[2]# the number of "sites" for each plot
            ot1=rowSums(ot)
            out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 9, length.out=9)), se=FALSE)
            richness[j]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
          }
        }
      },
      error = function(e) 
      {
        # Handle the error and continue
        print(paste("Error on element", j, ":", e))
      },
      finally ={
        print(paste("Processed element", j))}
      )
    }
    
    richness_plot_neon[k,1]=plot_ID[k]
    richness_plot_neon[k,2]=richness%>%mean(na.rm = TRUE)
    richness_plot_neon[k,3]=richness%>%sd(na.rm = TRUE)},
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", k, ":", e))
    },
    finally ={
      print(paste("Processed element", k))}
    )
  }
  data[[m]]=richness_plot_neon# the first two was compuated on great lakes
}


richness_plot_neon%>%bind_rows(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],data[[6]],data[[7]],data[[8]])->data_neon_30m_scale

save(data_neon_30m_scale,file="data_neon_30m_scale.RData")



##  whole-community fungal diversity at the 40m by 40m scale

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
  else if(dim1>=2&dim1<=15)
  {
    sub_samp=rownames(sample_data(data_sub))
    kk=subset_samples(rare_all_assign,sample_names(rare_all_assign)%in%sub_samp)
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

richness_mean_subplot40_neon=matrix(ncol=2,nrow=length(a4))
for (i in 1:length(a4))
{
  
  richness_mean_subplot40_neon[i,1]=mean(richness[[i]])
  richness_mean_subplot40_neon[i,2]=sd(richness[[i]])
}


richness_subplot40_neon_all_extra=cbind(plotid=a4,mean_value=richness_mean_subplot40_neon[,1],sd_value=richness_mean_subplot40_neon[,2],rep(1600,length.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value","area")))

save(richness_subplot40_neon_all_extra,file="richness_subplot40_neon_all_extra.RData")


# for the eight individual guilds at the 40 by 40 m scale

guild_select=c("ectomycorrhizal","arbuscular_mycorrhizal","soil_saprotroph","litter_saprotroph","plant_pathogen","wood_saprotroph",
               "para","epiphyte")

## use the extrapolation approach to determine the 40 m scale richness for the neon site

data=list()
for (k in 1:8){
  #cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
  
  data_select=subset_taxa(rare_all_assign,guild==guild_select[k])
  
  set.seed=(5679)
  times=30
  a4=sample_data(data_select)
  a4=unique(a4$plotIDM)
  richness <- vector("list", length(a4))
  
  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, plotIDM==a4[i])
    dim1=dim(sample_data(data_sub))[1]
    
    tryCatch({if(dim1>=16)
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
      else if(dim1>=2&dim1<=15)
      {
        sub_samp=rownames(sample_data(data_sub))
        kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
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
    },
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", i, ":", e))
    },
    finally ={
      print(paste("Processed element", i))}
    )
  }
  
  richness_mean_subplot40_neon=matrix(ncol=2,nrow=length(a4))
  for (i in 1:length(a4))
  {
    
    richness_mean_subplot40_neon[i,1]=mean(richness[[i]])
    richness_mean_subplot40_neon[i,2]=sd(richness[[i]])
  }
  
  
  data[[k]]=cbind(plotid=a4,mean_value=richness_mean_subplot40_neon[,1],sd_value=richness_mean_subplot40_neon[,2],rep(1600,length.out=n),rep(guild_select[k],length.out=n))%>%data.frame()%>%rename_all(~paste0(c("plotid","mean_value","sd_value","area","guild")))
}
##save the data
bind_rows(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],data[[6]],data[[7]],data[[8]])->temp

# need to 
richness_subplot40_neon_all_extra%>%mutate(guild=rep("all",n()))%>%
  bind_rows(temp)->data_neon_40m_scale

save(data_neon_40m_scale,file="data_neon_40m_scale.RData")# for both whole-community and individual guilds



######################for the dob plots#######################

####dob at the 10m by 10 m plot
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

# get the richness at the 10 by 10 m

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



####dob at the 20m by 20 m plot

set.seed=(1203)
times=30
a4=sample_data(dob_permutation)
a4=unique(a4$subplotID20)
richness <- vector("list", length(a4))

for (i in 1:length(a4))# note that the i here is large
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob_permutation, subplotID20==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  tryCatch({if(dim1>=4)
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
    else if (dim1>=2&dim1<4)
    {
      sub_samp=rownames(sample_data(data_sub))
      kk=subset_samples(dob_permutation,sample_names(dob_permutation)%in%sub_samp)
      data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
      ot=t(otu_table(data_sub))%>%data.frame()
      n=dim(ot)[2]# the number of "sites" for each plot
      ot1=rowSums(ot)
      out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 4, length.out=4)), se=FALSE)
      richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    }
    else{
      
      richness[[i]]=NA
    }
  },
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))}
  )
}


richness_mean_subplot20_dob=matrix(ncol=2,nrow=length(a4))
for (i in 1:length(a4))
{
  richness_mean_subplot20_dob[i,1]=mean(richness[[i]])
  richness_mean_subplot20_dob[i,2]=sd(richness[[i]])
}

richness_mean_subplot20_dob%>%data.frame()%>%bind_cols(a4)%>%mutate(plotid=substr(a4,1,8))%>%
  group_by(plotid)%>%summarise(mean_value=mean(X1,na.rm=TRUE),sd_value=sd(X1,na.rm=TRUE))%>%
  mutate(area=rep(400,n()),guild=rep("all",n()))->data_dob_20_all

# for individual guilds

data=list()
for (m in 1:8)
{
data_select=subset_taxa(dob_permutation,guild==guild_select[m])
  
set.seed=(1203)
times=30
a4=sample_data(data_select)
a4=unique(a4$subplotID20)
richness <- vector("list", length(a4))

for (i in 1:length(a4))# note that the i here is large
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(data_select, subplotID20==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  tryCatch({if(dim1>=4)
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
    else if (dim1>=2&dim1<4)
    {
      sub_samp=rownames(sample_data(data_sub))
      kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
      data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
      ot=t(otu_table(data_sub))%>%data.frame()
      n=dim(ot)[2]# the number of "sites" for each plot
      ot1=rowSums(ot)
      out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 4, length.out=4)), se=FALSE)
      richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    }
    else{
      
      richness[[i]]=NA
    }
  },
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))}
  )
}

richness_mean_subplot20_dob=matrix(ncol=2,nrow=length(a4))
for (i in 1:length(a4))
{
  richness_mean_subplot20_dob[i,1]=mean(richness[[i]])
  richness_mean_subplot20_dob[i,2]=sd(richness[[i]])
}

richness_mean_subplot20_dob%>%data.frame()%>%bind_cols(a4)%>%mutate(plotid=substr(a4,1,8))%>%
  group_by(plotid)%>%summarise(mean_value=mean(X1,na.rm=TRUE),sd_value=sd(X1,na.rm=TRUE))%>%
  mutate(area=rep(400,n()),guild=rep("all",n()))->data[[m]]

}






####### dob plots at the 30m by 30 m scale
### at the 30 m by 30 m scale
location=sample_data(dob_permutation)%>%data.frame()%>%select(X1,X2)%>%rename_all(~paste0(c("gx","gy")))

rownames(location)=rownames(sample_data(dob_permutation))

location=sample_data(location)

dob_permutation=merge_phyloseq(dob_permutation,location)

ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

###

term=c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite")

guild_select=c("ectomycorrhizal","arbuscular_mycorrhizal","soil_saprotroph","litter_saprotroph","plant_pathogen","wood_saprotroph",
               "para","epiphyte")

taxa_df <- as.data.frame(tax_table(dob_permutation))

guild=tax_table(dob_permutation)%>%data.frame()%>%select(ta2)

guild%>%mutate(ta2=ifelse(ta2%in%term,"para",ta2))%>%rename(guild=ta2)->guild

taxa_df$guild=guild

new_tax_table <- as.matrix(taxa_df)


## at the 30m by 30 m scale for the dob site

plot_ID=sample_data(dob_permutation)%>%data.frame()
plot_ID=unique(plot_ID$plotIDM)
set.seed(520)
richness_plot_dob=matrix(ncol=3,nrow=length(plot_ID))
for(k in 1:length(plot_ID))
{
  cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
  
  tryCatch({core_number=numeric()
  for(i in 1:dim(ori)[1])
  {
    #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    
    d=sample_data(dob_permutation)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==plot_ID[k])
    
    core_number[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
  }
  core_number_8=core_number>8
  core_number_8[!core_number_8]=NA
  
  
  richness=numeric()
  for(j in 1:length(core_number))# only select the locations that over 9 cores can be sampled/81 cores
  {
    cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
    tryCatch({if (j%in%c(which(!is.na(core_number_8)))){
      # informs the processing
      sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
      sub_samp=sample(rownames(sub_samp),9)# select only nine cores
      kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
      richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
    } 
      else{
        sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
        if(dim(sub_samp)[1]<2)
        {
          richness[j]=NA 
        }
        else
        {
          sub_samp=rownames(sub_samp)# get the sample names
          kk=subset_samples(dob_permutation,sample_names(dob_permutation)%in%sub_samp)#based on this to extrapolate the richness
          data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
          ot=t(otu_table(data_sub))%>%data.frame()
          n=dim(ot)[2]# the number of "sites" for each plot
          ot1=rowSums(ot)
          out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 9, length.out=9)), se=FALSE)
          richness[j]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
          
        }
      }
    },
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", j, ":", e))
    },
    finally ={
      print(paste("Processed element", j))}
    )
  }
  
  richness_plot_dob[k,1]=plot_ID[k]
  richness_plot_dob[k,2]=richness%>%mean(na.rm = TRUE)
  richness_plot_dob[k,3]=richness%>%sd(na.rm = TRUE)},
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", k, ":", e))
  },
  finally ={
    print(paste("Processed element", k))}
  )
}


#for individual guilds for the dob sites

data=list()
for (m in 1:8)
{
  
  data_select=subset_taxa(dob_permutation,guild==guild_select[m])
  
  
  plot_ID=sample_data(data_select)%>%data.frame()
  plot_ID=unique(plot_ID$plotIDM)
  set.seed(520)
  richness_plot_dob=matrix(ncol=3,nrow=length(plot_ID))
  for(k in 1:length(plot_ID))
  {
    cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
    
    tryCatch({core_number=numeric()
    for(i in 1:dim(ori)[1])
    {
      #cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
      
      d=sample_data(data_select)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==plot_ID[k])
      
      core_number[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
    }
    core_number_8=core_number>8
    core_number_8[!core_number_8]=NA
    
    
    richness=numeric()
    for(j in 1:length(core_number))# only select the locations that over 9 cores can be sampled/81 cores
    {
      cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
      tryCatch({if (j%in%c(which(!is.na(core_number_8)))){
        # informs the processing
        sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
        sub_samp=sample(rownames(sub_samp),9)# select only nine cores
        kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
        richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
      } 
        else{
          sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
          if(dim(sub_samp)[1]<2)
          {
            richness[j]=NA 
          }
          else
          {
            sub_samp=rownames(sub_samp)# get the sample names
            kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)#based on this to extrapolate the richness
            data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
            ot=t(otu_table(data_sub))%>%data.frame()
            n=dim(ot)[2]# the number of "sites" for each plot
            ot1=rowSums(ot)
            out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 9, length.out=9)), se=FALSE)
            richness[j]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
            
          }
        }
      },
      error = function(e) 
      {
        # Handle the error and continue
        print(paste("Error on element", j, ":", e))
      },
      finally ={
        print(paste("Processed element", j))}
      )
    }
    
    richness_plot_dob[k,1]=plot_ID[k]
    richness_plot_dob[k,2]=richness%>%mean(na.rm = TRUE)
    richness_plot_dob[k,3]=richness%>%sd(na.rm = TRUE)},
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", k, ":", e))
    },
    finally ={
      print(paste("Processed element", k))}
    )
  }
  
  
  data[[m]]=cbind(richness_plot_dob,guild_select[m])# actually this is for the dob sites
  
  
}

richness_plot_dob%>%data.frame()%>%mutate(area=rep(900,n()),guild=rep("all",n()))%>%
  rename_all(~paste0(c("plotid","mean_value","sd_value","area","guild")))->temp


bind_rows(data[[1]]%>%data.frame(),
          data[[2]]%>%data.frame(),
          data[[3]]%>%data.frame(),
          data[[4]]%>%data.frame(),
          data[[5]]%>%data.frame(),
          data[[6]]%>%data.frame(),
          data[[7]]%>%data.frame(),
          data[[8]]%>%data.frame())%>%mutate(area=rep(900,n()))%>%select(X1,X2,X3,area,X4)%>%
  rename_all(~paste0(c("plotid","mean_value","sd_value","area","guild")))%>%
  bind_rows(temp)->data_dob_30m_scale

save(data_dob_30m_scale,file="data_dob_30m_scale.RData")




### for the dob sites at the 40 by 40 m scale
dob=dob_permutation
## when all the guilds were considered
set.seed=(5679)
times=30
a4=sample_data(dob_permutation)
a4=unique(a4$plotIDM)
richness <- vector("list", length(a4))

for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(dob_permutation, plotIDM==a4[i])
  dim1=dim(sample_data(data_sub))[1]
  
  tryCatch({if(dim1>=16)
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
    
    else if(dim1>=2&dim1<=15)
    {
      sub_samp=rownames(sample_data(data_sub))
      kk=subset_samples(dob_permutation,sample_names(dob_permutation)%in%sub_samp)
      data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
      ot=t(otu_table(data_sub))%>%data.frame()
      n=dim(ot)[2]# the number of "sites" for each plot
      ot1=rowSums(ot)# when the sum is o, it incurs errors
      out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 16, length.out=16)), se=FALSE)
      richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
    }
    
    else{
      richness[[i]]=NA
    }
  },
  error = function(e) 
  {
    # Handle the error and continue
    print(paste("Error on element", i, ":", e))
  },
  finally ={
    print(paste("Processed element", i))}
  )
}


richness_mean_subplot40_dob=matrix(ncol=2,nrow=length(a4))
for (i in 1:length(a4))
{
  
  richness_mean_subplot40_dob[i,1]=mean(richness[[i]])
  richness_mean_subplot40_dob[i,2]=sd(richness[[i]])
}

richness_mean_subplot40_dob%>%data.frame()%>%mutate(plotid=a4,area=rep(1600,n()),guild=rep("all",n()))%>%
  rename_all(~paste0(c("mean_value","sd_value","plotid","area","guild")))%>%
  select(plotid,mean_value,sd_value,area,guild)->data_dob_40m_all

save(data_dob_40m_all,file="data_dob_40m_all.RData")

###
taxa_df <- as.data.frame(tax_table(dob_permutation))

guild=tax_table(dob_permutation)%>%data.frame()%>%select(ta2)

guild%>%mutate(ta2=ifelse(ta2%in%term,"para",ta2))%>%rename(guild=ta2)->guild

taxa_df$guild=guild

new_tax_table <- as.matrix(taxa_df)

# Assign the modified taxonomy table back to the phyloseq object
tax_table(dob_permutation) <- new_tax_table


term=c("protistan_parasite","lichen_parasite","algal_parasite","mycoparasite","animal_parasite")

guild_select=c("ectomycorrhizal","arbuscular_mycorrhizal","soil_saprotroph","litter_saprotroph","plant_pathogen","wood_saprotroph",
               "para","epiphyte")

data=list()

for (k in 1:8){
  cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
  
  data_select=subset_taxa(dob_permutation,guild==guild_select[k])
  set.seed=(5679)
  times=30
  a4=sample_data(data_select)
  a4=unique(a4$plotIDM)
  richness <- vector("list", length(a4))
  
  for (i in 1:length(a4))
  {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    data_sub <- subset_samples(data_select, plotIDM==a4[i])
    dim1=dim(sample_data(data_sub))[1]
    
    tryCatch({if(dim1>=16)
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
      
      else if(dim1>=2&dim1<=15)
      {
        sub_samp=rownames(sample_data(data_sub))
        kk=subset_samples(data_select,sample_names(data_select)%in%sub_samp)
        data_sub =transform_sample_counts(kk, function(x) ifelse(x>0, 1, 0))
        ot=t(otu_table(data_sub))%>%data.frame()
        n=dim(ot)[2]# the number of "sites" for each plot
        ot1=rowSums(ot)# when the sum is o, it incurs errors
        out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(1, 16, length.out=16)), se=FALSE)
        richness[[i]]=out3$iNextEst$size_based%>%slice_tail()%>%pull(qD)
      }
      
      else{
        richness[[i]]=NA
      }
    },
    error = function(e) 
    {
      # Handle the error and continue
      print(paste("Error on element", i, ":", e))
    },
    finally ={
      print(paste("Processed element", i))}
    )
  }
  
  
  richness_mean_subplot40_dob=matrix(ncol=2,nrow=length(a4))
  for (i in 1:length(a4))
  {
    
    richness_mean_subplot40_dob[i,1]=mean(richness[[i]])
    richness_mean_subplot40_dob[i,2]=sd(richness[[i]])
  }
  
  
  data[[k]]=cbind(plotid=a4,mean_value=richness_mean_subplot40_dob[,1],sd_value=richness_mean_subplot40_dob[,2],rep(1600,length.out=n),rep(guild_select[k],length.out=n))%>%data.frame()%>%
    
    rename_all(~paste0(c("plotid","mean_value","sd_value","area","guild")))
  
}





