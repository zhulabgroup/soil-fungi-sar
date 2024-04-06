###
library(phyloseq)
library(dplyr)
### to devide the plots into four quardrants

d=sample_data(rare_all)
# select the NEON SITE

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

a1= subset_samples(d,Project=="NEON")# the unique plotID, we have 476 plots

data_sub <- sample_data(a1)%>%data.frame()

# get the subplot for each plot
subplot=strsplit(data_sub$geneticSampleID,"-")
subplot_matrix=matrix(ncol=6,nrow=dim(data_sub)[1])
for(i in 1:dim(data_sub)[1]){
  subplot_matrix[i,]=subplot[[i]]
}
subplot_matrix%>%data.frame()->subplot_ID

subplot_ID=subplot_ID[,3:4]

names(subplot_ID)=c("gx","gy")
subplot_ID$gx=as.numeric(subplot_ID$gx)

subplot_ID$gy=as.numeric(subplot_ID$gy)

# to assign each row to a subplot clockwisely 1,2,3,4 quardrant.


point_assign=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=0 && rowid$gx<=20 && rowid$gy>=0 && rowid$gy<=20)
  {
    point_assign[i]=1
  }
  else if(rowid$gx>=0 && rowid$gx<=20 && rowid$gy>20 && rowid$gy<40)
  {
    point_assign[i]=2
  }
  else if(rowid$gx>20 && rowid$gx<=40 && rowid$gy>20 && rowid$gy<40)
  {
    point_assign[i]=3
  }
  else
  {
    point_assign[i]=4
  }
}
subplot_ID=cbind(subplot_ID,point_assign)
# to see if it is possible to get a soil core that are at the 10 m level


# add the subplot ID to the phyloseq subject.
row.names(subplot_ID)=row.names(sample_data(a1))# need to be the same with the sample data
subplot_ID=sample_data(subplot_ID)
a1=merge_phyloseq(a1,subplot_ID)

# to see how many samples are there for each subplot for each plot
a2=sample_data(a1)

a2=unique(a2$plotIDM)

point_number=matrix(ncol=4,nrow=length(a2))

for (i in 1:length(a2))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_subb <- subset_samples(a1, plotIDM==a2[i])
  data_subb=sample_data(data_subb)%>%data.frame()
  
  ncores=table(data_subb$point_assign)%>%data.frame()
  ncores$Var1=as.factor(ncores$Var1)
  com_ple=data.frame(Var1=c(1:4),Freq=c(1,1,1,1))# to merge in case of NA
  
  com_ple$Var1=as.factor(com_ple$Var1)
t(left_join(com_ple,ncores,by="Var1")) %>%data.frame()->temp
  
  point_number[i,]=temp[3,]%>%as.numeric()
}

# the number of cores are quite similar among quandrants with about 3 points per/subplot
# create a id with a combination of both plotID and quardrant code
a3=sample_data(a1)
plotquaID=paste(a3$plotIDM,"-",a3$point_assign)%>%data.frame()
names(plotquaID)="plotquaID"

plotquaID$plotquaID=gsub(" ","",plotquaID$plotquaID)
# combind this unique ID with the whole sample data
row.names(plotquaID)=row.names(a3)
plotquaID=sample_data(plotquaID)

a1=merge_phyloseq(a1,plotquaID)
# get the total richness for each subplots

a4=sample_data(a1)
a4=unique(a4$plotquaID)

richness=numeric()
for (i in 1:length(a4))
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotquaID==a4[i])
  richness[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}

richness=data.frame(plotquaID=a4,richness=richness)# at the subplot level
richness$richness=as.numeric(richness$richness)
plotid=substr(richness$plotquaID,1,8)
richness=cbind(richness,plotid)

subplot_mean=aggregate(richness~plotid,data=richness,FUN=mean)
subplot_sd=aggregate(richness~plotid,data=richness,FUN=sd)

# we can exclude the subplot that have few cores to reduce the variation of subplot level richness

a4=sample_data(a1)
a4=unique(a4$plotquaID)

richness=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotquaID==a4[i])
  ncor=
  if (dim(sample_data(data_sub))[1]>2)
      {
        richness[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric() 
      }
else{
  richness[i]=NA
}
 
}

richness_morethan2=data.frame(plotquaID=a4,richness=richness)# more than 2 cores were selected

richness_morethan2$richness=as.numeric(richness_morethan2$richness)

plotid=substr(richness_morethan2$plotquaID,1,8)
richness_morethan2=cbind(richness_morethan2,plotid)

subplot_mean_morethan2=aggregate(richness~plotid,data=richness_morethan2,FUN=mean)

subplot_sd_morethan2=aggregate(richness~plotid,data=richness_morethan2,FUN=sd)


## the soil core level richness
# get the plot id 



## get the soil core level mean for the NEON sites

a4=sample_data(a1)
a4=unique(a4$plotIDM)

richness=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
  core_richness=estimate_richness(data_sub,measures = "Observed")
  richness[i]=mean(core_richness$Observed)
}

core_mean=data.frame(a4,richness)
### the plot level richness

a4=sample_data(a1)
a4=unique(a4$plotIDM)

richness=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
  
  total_richness=table(colSums(otu_table(data_sub))>0)["TRUE"]%>%as.numeric()
  
  richness[i]=total_richness
}

plot_total=data.frame(a4,richness)

##cbind all the three richness at three levels

names(subplot_mean)[1]="a4"

df=left_join(core_mean,subplot_mean,by="a4")

df=left_join(df,plot_total,by="a4")

names(df)=c("plotID","corerich","subplotrich","plotrich")

dk=melt(df)
area=rep(c(0.0036,400,1600),each=472)

dk=cbind(dk,area)

save(dk,file="plot_SAR.RData")

# combine the data with subplot withmore than 2 soil core
names(subplot_mean_morethan2)[1]="a4"

df_morethan2=left_join(core_mean,subplot_mean_morethan2,by="a4")

df_morethan2=left_join(df_morethan2,plot_total,by="a4")%>%setNames(c("plotID","corerich","subplotrich","plotrich"))

df_morethan2=melt(df_morethan2)
area=rep(c(0.0036,400,1600),each=472)

df_morethan2=cbind(df_morethan2,area)
# get the z value

df_morethan2_nona=df_morethan2%>%filter(!is.na(value))

zvalue=numeric()
pvalue=numeric()

a5=unique(df_morethan2_nona$plotID)

for (i in 1:length(a5))
{
  df1=subset(df_morethan2_nona,plotID==a5[i])
  
  ft=lm(log(value)~log(area),data=df1)%>%summary()
  
  zvalue[i]=ft$coefficients[2,1]
  pvalue[i]=ft$coefficients[1,4]
}




##3#
zvalue=numeric()
pvalue=numeric()
a5=unique(dk$plotID)
for (i in 1:472)
  {
  df1=subset(dk,plotID==a5[i])
  ft=lm(log(value)~log(area),data=df1)%>%summary()
  zvalue[i]=ft$coefficients[2,1]
  pvalue[i]=ft$coefficients[1,4]
}



####
for (i in 1:472)
{
  df1=subset(dk,plotID==a5[i])
  ft=plot(value~area,data=df1)
}

dd=core_number>13

for (i in dd ){
  # Generate random data for demonstration
  df1=subset(dk,plotID==a5[i])
  
  # Create a scatterplot
  plot <- ggplot(data = df1, aes(x = log(area), y = log(value))) +
    geom_point(size=4) +
    labs(title = paste("Plot", i), x = "are", y = "species")
  
  # Print the plot
  print(plot)
}
# to see the number of cores within each plot

a4=sample_data(a1)
a4=unique(a4$plotIDM)

core_number=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
      
  
  core_number[i]=dim(sample_data(data_sub))[1]
}
