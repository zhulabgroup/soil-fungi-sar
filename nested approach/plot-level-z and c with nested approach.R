
# Nested method for the NEON site,which means that the prior sites were included when adding the new ones

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

# get the data of each of the permutation

set.seed(1015)
perdf <- list()
for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  
  data_sub <- subset_samples(d, plotIDM==a1[i])
  
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression with at least three cores
  {
    oot=otu_table(data_sub)
    oot=specaccum(oot,method="random")
    
    perdf[[i]]=oot$perm
    
  }
  
  else
    {
    perdf[[i]]=NA
  }
}


#

co=numeric()
for (i in 1:515)
  {
 co[i] =dim(as.data.frame(perdf[[i]]))[2]
}

# based on the permutation to get the mean c and the mean z for each plot



## for the first site
zzk=list()
  
for (i in 1:515)
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing 

  a5= data.frame(perdf[[i]])

  if (dim(a5)[1]>4) {
    
zz=matrix(ncol=2,nrow=100)

for (j in 1:100)
{
  ex<- as.data.frame(cbind("A"=c(1:dim(a5)[1]),species=a5[,j]))
  temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]
  zz[j,1]=temp[2,1]# the estimated z
  zz[j,2]=temp[1,1]# # the estimated C
}

zzk[[i]]=zz
  }  
  
  else if (dim(a5)[1]==4){
    
    zz=matrix(ncol=2,nrow=23)
    
    for (j in 1:23)
    {
      ex<- as.data.frame(cbind("A"=c(1:dim(a5)[1]),species=a5[,j]))
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]
      zz[j,1]=temp[2,1]# the estimated z
      zz[j,2]=temp[1,1]# # the estimated C
    }
    
    zzk[[i]]=zz
    
  }

  else if (dim(a5)[1]==3){
    
    zz=matrix(ncol=2,nrow=5)
    
    for (j in 1:5)
    {
      ex<- as.data.frame(cbind("A"=c(1:dim(a5)[1]),species=a5[,j]))
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]
      zz[j,1]=temp[2,1]# the estimated z
      zz[j,2]=temp[1,1]# # the estimated C
    }
    
    zzk[[i]]=zz
    
  }
  
}

## get the mean and sd of the estimated z and c value
op=matrix(ncol=2,nrow=515)
op_sd=matrix(ncol=2,nrow=515)
for (i in 1:515)
  {
  if(is.null(zzk[[i]]))
  {
    op[i,1]=NA
    op[i,2]=NA
    op_sd[i,1]=NA
    op_sd[i,2]=NA
    
  }
  else
    {
    a6=apply(zzk[[i]],2,FUN=mean)
    a7=apply(zzk[[i]],2,FUN=sd)
    op[i,1]=a6[1]
    op[i,2]=a6[2] 
    op_sd[i,1]=a7[1]
    op_sd[i,2]=a7[2] 
  }

}


### add the plot id
plot_level_zc_nest_mean=cbind(a1,op)
plot_level_zc_nest_sd=cbind(a1,op_sd)

save(plot_level_zc_nest_mean,file="plot_level_zc_nest_mean.RData")
save(plot_level_zc_nest_sd,file="plot_level_zc_nest_sd.RData")



