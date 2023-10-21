# estimate the z values based on functional guilds
ft=read.csv("FungalTraits_1.2_ver_16Dec_2020.csv",sep=",",header=T)
ft=ft[,c("GENUS","primary_lifestyle")]
# get the genus name of all the taxa in the full dataset
d=tax_table(rare_all)[,6]
d=data.frame(d)

names(ft)=c("genus","guild")
# assign guild types to each taxa
d=left_join(d,ft,by="genus")
# change it to the tax-table format
d=tax_table(d)
row.names(d)=row.names(tax_table(rare_all))
#adding the guild column to the existing phyloseq object
rare_all=merge_phyloseq(rare_all,d)

### choose a specific functional guild and to see how many plots it occurs


guild=tax_table(rare_all)
guild=data.frame(guild)
guild=unique(guild$ta2)
guild=guild[c(1,3:24,26:29)]# removed two categories of "" and 'NA'

b=vector("list",length(guild))
for(i in 1:length(guild))
  {
  a=subset_taxa(rare_all,ta2==guild[i])# must first select the guild and then check the sample and taxa sums
  a<- subset_samples(a , !is.na(lon) & !is.na(lat))
  a <- subset_taxa(a, taxa_sums(a) > 0)
  a <- subset_samples(a, sample_sums(a) > 0)
  b[[i]]=dim(otu_table(a))
}

c=matrix(nrow=27,ncol=2)
for( i in 1:27)
  {
  c[i,1]=b[[i]][1]
  c[i,2]=b[[i]][2]
}
cbind(guild,c)# look at the number of taxa for each guild and how many cores they occur

#some of the guilds are very rich in taxa, need to determine guild-specific z values
# i would not use Qin's criteria to group guilds as saprotroph on the litter, soil and wood my differ widely in their life strateges

#### for the "soil_saprotroph"
m=subset_taxa(rare_all,ta2=="soil_saprotroph")# select a guild
m<- subset_samples(rare_all, !is.na(lon) & !is.na(lat))
m<- subset_taxa(m, taxa_sums(m) > 0)
m<- subset_samples(m, sample_sums(m) > 0)


# add a plotID for the data

d=sample_data(m)
table(d$Project)# the first 908 rows are dob sites while the remaining 5470 rows are NEON sites
d1=data.frame(d[,c("geneticSampleID","Site")])
plotID=substr(d1$geneticSampleID,1,8)
d1=cbind(d1,plotID)
iddob=d1$Site[1:table(d$Project)[[1]]]#a site corresponds to a plot
idneon=d1$plotID[(table(d$Project)[[1]]+1):sum(table(d$Project))]# a plotID corresponds to a plot
plotIDM=data.frame(c(iddob,idneon))

names(plotIDM)="plotIDM"# the plotID used for the SAR
row.names(plotIDM)=row.names(sample_data(m))#an important step to merge the data
plotIDM=sample_data(plotIDM)
d<- merge_phyloseq(m, plotIDM)# merge the new plotid with the initial data 

a1= sample_data(d)
a1=unique(a1$plotIDM)

set.seed(10101)
times=30
power.z <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 3)# we can construct a linear regression model with at least three cores
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
      species <- vector(length = dim1[1]) # create a vector to save diversity
      for (j in 1:dim1[1]) 
      { 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], j)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
      }
      ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]# extract the estimated c and z values
      return(c(temp[2,1], temp[1,1]))
    }
    stopCluster(cl)
  }
  else
  {
    power.z[[i]]=matrix(10,nrow=2,ncol = dim1[1])#the number 10 is randomly selected to create a matrix for the plots with < 3 cores to avoid NULL output
  }
}


