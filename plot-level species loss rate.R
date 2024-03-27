# for one site, we can estimated the sar by decreasing the soil core.
d=sample_data(rare_all)
table(d$Project)# the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:908]#a site correspondes to a plot
idneon=d1$plotID[909:6378]# an unique plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))
names(plotIDM)="plotIDM"# the plot id used for the SAR
row.names(plotIDM)=row.names(d)
plotIDM=sample_data(plotIDM)
d<- merge_phyloseq(rare_all, plotIDM)# merge the new plotid with the initial data 
# select an unique plot and build a SAR within the plot
a1= sample_data(d)# the unique plotID, we have 476 plots
a1=unique(a1$plotIDM)

set.seed(10105)
times=30
loss.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression with at least three cores
  {
    cl <- makeCluster(2)
    registerDoParallel(cl)
    loss.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
      
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) 
      { 
        # randomly sample j samples in the plot
        # the total species
        to=taxa_sums(data_sub)
        to=to>0
        to=sum(to,na.rm=TRUE)
        
        flag <- rep(TRUE, dim1[1])
        flag[sample(1:(dim1[1]-1), j-1)] <- FALSE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j-1] <- sum(otu_table(temp)["TRUE"] > 0)
        species[j]=to[1]
      }
      
      ex <- data.frame("A"=c(1:(dim1[1]-1),0),sp=species)
      
      temp<- summary(lm(sp~A,ex))[["coefficients"]]# extract the estimated c and z values
      return(c(temp[2,1], temp[1,1]))
    }
    stopCluster(cl)
  }
  
  else
  {
    loss.z[[i]]=matrix(10,nrow=2,ncol = dim1[1])#the number 10 is randomly selected, to create a matrix for the plots with < 3 cores to avoid NULL output
  }
}

all_z_loss=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for all the plots
for(i in 1:length(a1)){
  all_z_loss[i,]=loss.z[[i]][1,]
}

all_z_loss=data.frame(a1,all_z_loss)

save(all_z_loss,file="all_z_loss.ranall.RData")

all_c_loss=matrix(nrow=length(a1),ncol=30)# get the 30 simulated c values for the neon site
for(i in 1:length(a1)){
  all_c_loss[i,]=loss.z[[i]][2,]
}

all_c_loss=data.frame(a1,all_c_loss)

save(all_c_loss,file="all_c_loss.ranall.RData")


# to check the relationship between the z and the loss rate

all_z_loss=all_z_loss[,-2]
all_z_loss_mean=data.frame(plotID=a1,rate=apply(all_z_loss[,2:31],1,FUN=mean))


all_z_incre_mean=aggregate(z~plotID,data=model_data,FUN=mean)

comp_z_loss=merge(all_z_loss_mean,all_z_incre_mean,by="plotID")

core_number=df30mean[,c(1,6)]%>%unique()
names(core_number)[1]="plotID"

comp_z_loss=merge(comp_z_loss,core_number,by="plotID")

save(comp_z_loss,file="comp_z_loss.RData")

#### when the same number of cores were sampled within each plot

## a function for random sample

sample_ps <- function(ps, FUN = sample, ...){
  ids <- sample_names(ps)
  sampled_ids <- FUN(ids, ...)
  ps <- prune_samples(sampled_ids, ps)
  return(ps)
}

set.seed(10106)
times=30
loss.z10 <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 10)# we can construct a linear regression with at least three cores
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    loss.z10[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
     
      ka=sample_ps(data_sub,size=10)
       
      species <- vector(length = 10) # create a vector to save diversity
      for (j in 1:9) 
      { 
        to=taxa_sums(ka)
        to=to>0
        to=sum(to,na.rm=TRUE)
        

        flag <- rep(TRUE, 10)
        flag[sample(1:9, j)] <- FALSE
        temp <- merge_samples(ka, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
        species[10] <- to[1]
      }
      
      ex <- data.frame("A"=c(1:9,0),sp=species)
      
      temp<- summary(lm(sp~A,ex))[["coefficients"]]# extract the estimated c and z values
      return(c(temp[2,1], temp[1,1]))
    }
    stopCluster(cl)
  }
  
  else
  {
    loss.z10[[i]]=matrix(10,nrow=2,ncol = 10)#the number 10 is randomly selected, to create a matrix for the plots with < 3 cores to avoid NULL output
  }
}

all_z_loss.z10=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for all the plots
for(i in 1:length(a1)){
  all_z_loss.z10[i,]=loss.z10[[i]][1,]
}

save(all_z_loss.z10,file="all_z_loss.z10.RData")
all_c_loss.z10=matrix(nrow=length(a1),ncol=30)# get the 30 simulated c values for the neon site
for(i in 1:length(a1)){
  all_c_loss.z10[i,]=loss.z10[[i]][2,]
}
save(all_c_loss.z10,file="all_c_loss.z10.RData")
## get the mean
all_z_loss.z10_mean=data.frame(plotID=a1,rate10=apply(all_z_loss.z10,1,mean))

comp_z_loss= merge(all_z_loss.z10_mean,comp_z_loss,by="plotID")



