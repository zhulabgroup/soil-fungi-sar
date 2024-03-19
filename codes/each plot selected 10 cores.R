# step 1 calculates the z and c values for both the neon and dob sites and the resulting outputs were saved locally
# this step firstly involves the rarefation of the data, which will be used in the downstream analyses when necessary.
library(doParallel)
library(phyloseq)
rarefy all the data
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="OH")# the data only include the O and M soil horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob<- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob<- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all=rarefy_even_depth(neon_dob, rngseed=10,sample.size = 3000, replace = F)#764 samples were removed

save(rare_all,file="rare_all.Rdata")# save the rarefied data

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

set.seed(10108)
times=30
power.z00 <- vector("list", length(a1))

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  dim1 <- dim(otu_table(data_sub)) # the number of samples in one site
  if (dim1[1] >= 10)# we can construct a linear regression with at least three cores
  {
    cl <- makeCluster(3)
    registerDoParallel(cl)
    power.z00[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar%{
      species <- vector(length = 10) # create a vector to save diversity
      for (j in 1:10) 
      { 
        cat('\r',paste(paste0(rep("*", round(j/ 1, 0)), collapse = ''), j, collapse = ''))# 
        # randomly sample j samples in the plot
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], j)] <- TRUE
        temp <- merge_samples(data_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
      }
      ex <- as.data.frame(cbind("A"=c(1:10),species))
      temp<- summary(lm(log(species)~log(A),ex))[["coefficients"]]# extract the estimated c and z values
      return(c(temp[2,1], temp[1,1]))
    }
    stopCluster(cl)
  }
  
  else
  {
    power.z00[[i]]=matrix(10,nrow=2,ncol = 10)#the number 10 is randomly selected, to create a matrix for the plots with < 3 cores to avoid NULL output
  }
}

all_z=matrix(nrow=length(a1),ncol=30)# get the 30 simulated z values for all the plots
for(i in 1:length(a1)){
  all_z[i,]=power.z00[[i]][1,]
}

all_z=data.frame(all_z)

all_z=data.frame(a1,all_z)

write.csv(all_z,"all_z.ran30_10cores.csv")

# get the means for each plot

all_z_eq10=apply(all_z[,2:31],1,data=all_z,FUN=mean)%>%data.frame()%>%cbind(a1)

names(all_z_eq10)=c("z","piddd")
# look at the z values determined with two different approaches

head(model_data)
all_z_ran=aggregate(z~piddd,data=model_data,FUN=mean)
com_z=merge(all_z_eq10,all_z_ran,by="piddd")
names(com_z)=c("piddd","z10","zall")
  
com_z=merge(com_z, df,by="piddd")
com_z=subset(com_z,z10<10)
# compare different approach

app1=melt(com_z[,c(1,2,3)])
app2=melt(com_z[,c(1,4,4)])
app3=cbind(app1,app2["value"])
names(app3)[4]="core"

p1=ggplot()+
geom_point(data=app3,pch=21,color="black",size=3,aes(x=core,y=value,fill=variable))+
  theme(legend.position = c(0.88, 0.182027558), 
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(face = "italic", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA),
        panel.background = element_rect(fill = "NA"))+
  geom_smooth(data=app3,aes(x=core,y=value,color=variable),
              method = "lm",se=FALSE,size=1)+
  xlab("Number of soil cores")+
  ylab("z")+
  scale_fill_manual("",breaks=c("z10","zall"),labels=c("z10","zall"),values=c("seagreen1","mediumpurple"))+
scale_color_manual("",breaks=c("z10","zall"),labels=c("z10","zall"),values=c("seagreen1","mediumpurple"))


p2=ggplot()+
  geom_point(data=com_z,pch=21,fill="gray",size=3,aes(x=z10,y=zall))+
  theme(legend.position = c(0.88, 0.182027558), 
        legend.text = element_text(size = 14), 
        text = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(face = "italic", size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", size = 1.5, fill = NA),
        panel.background = element_rect(fill = "NA"))+
  ylab(expression(italic(z)["All cores sampled"]))+
  xlab(expression(italic(z)["10 cores sampled"]))


