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

# to get the richness at the scale of 30 by 30 m, which could be determined by four directions

point_assign1=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=0 && rowid$gx<=30 && rowid$gy>=0 && rowid$gy<=30)
  {
    point_assign1[i]=1
  }
  
  else
  {
    point_assign1[i]=4
  }
}
# the second approach
point_assign2=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=0 && rowid$gx<=30 && rowid$gy>=10 && rowid$gy<=40)
  {
    point_assign2[i]=1
  }
  
  else
  {
    point_assign2[i]=4
  }
}

# the third approach

point_assign3=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=10 && rowid$gx<=40 && rowid$gy>=10 && rowid$gy<=40)
  {
    point_assign3[i]=1
  }
  
  else
  {
    point_assign3[i]=4
  }
}

# the fourth approach

point_assign4=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=10 && rowid$gx<=40 && rowid$gy>=0 && rowid$gy<=30)
  {
    point_assign4[i]=1
  }
  
  else
  {
    point_assign4[i]=4
  }
}

subplot_ID_30=cbind(subplot_ID,point_assign1,point_assign2,point_assign3,point_assign4)

subplot_ID_30=subplot_ID_30[,-c(1:2)]

row.names(subplot_ID_30)=rownames(sample_data(a1))

subplot_ID_30=sample_data(subplot_ID_30)

a1=merge_phyloseq(a1,subplot_ID_30)# add to the phyloseq object
# need to combined the plotID
AA=data.frame(sample_data(a1)

AA=substr(AA$plotIDM,1,8)

subplot_ID_30_A=paste(AA,"-",point_assign1)
subplot_ID_30_B=paste(AA,"-",point_assign2)
subplot_ID_30_C=paste(AA,"-",point_assign3)
subplot_ID_30_D=paste(AA,"-",point_assign4)

AB=data.frame(subplot_ID_30_A=subplot_ID_30_A,subplot_ID_30_B=subplot_ID_30_B,subplot_ID_30_C=subplot_ID_30_C,subplot_ID_30_D=subplot_ID_30_D)
row.names(AB)=row.names(sample_data(a1))
AB=sample_data(AB)

a1=merge_phyloseq(a1,AB)


# to get the total richness within each subplots with a size of 30 by 30

a4=sample_data(a1)
a4=unique(a4$subplot_ID_30_A)

richness=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_A==a4[i])
  richness[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}
# for the second approach

SP_A=data.frame(cbind(a4,richness)

SP_B=data.frame(cbind(a4,richness2)

SP_C=data.frame(cbind(a4,richness3))

SP_D=data.frame(cbind(a4,richness4))             

# get the mean of the subplot

names(SP_A)[2]="richness"
names(SP_B)[2]="richness"
names(SP_C)[2]="richness"

names(SP_D)[2]="richness"

SP_30subplot=rbind(SP_A,SP_B,SP_C,SP_D)
SP_30subplot$richness=as.numeric(SP_30subplot$richness)

SP_30subplot_mean=aggregate(richness~a4,data=SP_30subplot,FUN=mean)

# get the rows that contain 1
plot_need=numeric()
for (i in 1:dim(SP_30subplot_mean)[1]){
  dm=SP_30subplot_mean[i,]
if (str_detect(dm$a4," - 1"))
    {
 plot_need[i]=1
}
else{
  plot_need[i]=NA
}
}

SP_30subplot_mean=cbind(SP_30subplot_mean,plot_need)

SP_30subplot_mean=SP_30subplot_mean%>%filter(!is.na(plot_need))


a4=sample_data(a1)
a4=unique(a4$subplot_ID_30_B)


 
richness2=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_B==a4[i])
  richness2[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}


# for the third apporach

a4=sample_data(a1)
a4=unique(a4$subplot_ID_30_C)

richness3=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_C==a4[i])
  richness3[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}

# for the fourth approach
a4=sample_data(a1)
a4=unique(a4$subplot_ID_30_D)

richness4=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, subplot_ID_30_D==a4[i])
  richness4[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}

# add the subplot ID to the phyloseq subject.
row.names(subplot_ID)=row.names(sample_data(a1))# need to be the same with the sample data
subplot_ID=sample_data(subplot_ID)
a1=merge_phyloseq(a1,subplot_ID)

# to see if it is possible to get a soil core that are at the 10 m level
# we will have 16 grids
# from bottom to the top

point_assign=numeric()
for(i in 1:dim(subplot_ID)[1])
{
  rowid=subplot_ID[i,]
  if(rowid$gx>=0 && rowid$gx<=10 && rowid$gy>=0 && rowid$gy<=10)
  {
    point_assign[i]=1
  }
  else if(rowid$gx>=0 && rowid$gx<=10 && rowid$gy>10 && rowid$gy<20)
  {
    point_assign[i]=2
  }
  else if(rowid$gx>0 && rowid$gx<=10 && rowid$gy>20 && rowid$gy<30)
  {
    point_assign[i]=3
  }
  
  else if(rowid$gx>0 && rowid$gx<=10 && rowid$gy>30 && rowid$gy<40)
  {
    point_assign[i]=4
  }
  
  else if(rowid$gx>10 && rowid$gx<=20 && rowid$gy>0 && rowid$gy<10)
  {
    point_assign[i]=5
  }
  else if(rowid$gx>10 && rowid$gx<=20 && rowid$gy>10 && rowid$gy<20)
  {
    point_assign[i]=6
  }
  
  else if(rowid$gx>10 && rowid$gx<=20 && rowid$gy>20 && rowid$gy<30)
  {
    point_assign[i]=7
  }
  else if(rowid$gx>10 && rowid$gx<=20 && rowid$gy>30 && rowid$gy<40)
  {
    point_assign[i]=8
  }
  else if(rowid$gx>20 && rowid$gx<=30 && rowid$gy>0 && rowid$gy<10)
  {
    point_assign[i]=9
  }
  
  else if(rowid$gx>20 && rowid$gx<=30 && rowid$gy>10 && rowid$gy<20)
  {
    point_assign[i]=10
  }
  else if(rowid$gx>20 && rowid$gx<=30 && rowid$gy>20 && rowid$gy<30)
  {
    point_assign[i]=11
  }
  else if(rowid$gx>20 && rowid$gx<=30 && rowid$gy>30 && rowid$gy<40)
  {
    point_assign[i]=12
  }
  
  else if(rowid$gx>30 && rowid$gx<=40 && rowid$gy>10 && rowid$gy<20)
  {
    point_assign[i]=13
  }
  else if(rowid$gx>30 && rowid$gx<=40 && rowid$gy>20 && rowid$gy<30)
  {
    point_assign[i]=14
  }
  else if(rowid$gx>30 && rowid$gx<=40 && rowid$gy>30 && rowid$gy<40)
  {
    point_assign[i]=15
  }
  
  else
  {
    point_assign[i]=16
  }
}


subplot_ID_5=data.frame(point_assign)#at the 5 m scale

names(subplot_ID_5)="point_assign5"

row.names(subplot_ID_5)=row.names(sample_data(a1))# need to be the same with the sample data

subplot_ID_5=sample_data(subplot_ID_5)

a1=merge_phyloseq(a1,subplot_ID_5)


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

## for the quardrat of 10 by 10 m

a3=sample_data(a1)
plotquaID5=paste(a3$plotIDM,"-",a3$point_assign5)%>%data.frame()

names(plotquaID5)="plotquaID5"

plotquaID5$plotquaID5=gsub(" ","",plotquaID5$plotquaID5)

# combind this unique ID with the whole sample data

row.names(plotquaID5)=row.names(sample_data(a1))

plotquaID5=sample_data(plotquaID5)

a1=merge_phyloseq(a1,plotquaID5)



# get the total richness for each subplots

a4=sample_data(a1)
a4=unique(a4$plotquaID5)

richness=numeric()
for (i in 1:length(a4))
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotquaID5==a4[i])
  richness[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}

richness=data.frame(plotquaID5=a4,richness=richness)# at the 10 by 10 subplot level

richness$richness=as.numeric(richness$richness)
plotid=substr(richness$plotquaID,1,8)
richness=cbind(richness,plotid)

subplot_mean=aggregate(richness~plotid,data=richness,FUN=mean)
subplot_sd=aggregate(richness~plotid,data=richness,FUN=sd)

# get the total richness for each 10 by 10 subplots

a4=sample_data(a1)
a4=unique(a4$plotquaID5)

richness=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotquaID5==a4[i])
  richness[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  
}

richness=data.frame(plotquaID=a4,richness=richness)# at the subplot level
richness$richness=as.numeric(richness$richness)
plotid=substr(richness$plotquaID,1,8)
richness=cbind(richness,plotid)

subplot_mean=aggregate(richness~plotid,data=richness,FUN=mean)
subplot_sd=aggregate(richness~plotid,data=richness,FUN=sd)
# to see how many cores are included in the 10 by 10 subplot

a4=sample_data(a1)
a4=unique(a4$plotquaID5)

richness_morethan2_10=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotquaID5==a4[i])
  if (dim(sample_data(data_sub))[1]>1)# selecte the subplot more than one core
  {
    richness_morethan2_10[i]=table(taxa_sums(data_sub)>0)["TRUE"]%>%as.numeric()
  }
  else{
    richness_morethan2_10[i]=NA
  } 
  
}

# get the mean value

richness_morethan2_10_mean=aggregate(richness~plotid,data=richness_morethan2_10,FUN=mean)

plotid10=substr(richness_morethan2_10$plotquaID5,1,8)

richness_morethan2_10=richness_morethan2_10[,-3]

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

# get the mean value of richness

richness_morethan2_10=cbind(plotquaID5=a4,richness=richness_morethan2_10)%>%data.frame()
richness_morethan2_10=richness_morethan2_10[,-2]

richness_morethan2_10$richness=as.numeric(richness_morethan2_10$richness)


plotid=substr(richness_morethan2_10$plotquaID5,1,8)

richness_morethan2_10=cbind(richness_morethan2_10,plotid)
# the mean richness within the 10 by 10 m subplots
subplot_mean_morethan2_10=aggregate(richness~plotid,data=richness_morethan2_10,FUN=mean)# we have 379 values

subplot_sd_morethan2_10=aggregate(richness~plotid,data=richness_morethan2_10,FUN=sd)# we have 379 values



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

for (i in 1:length(a2) )
  {
  data_subb <- subset_samples(a1, plotIDM==a2[i])
  data_subb=sample_data(data_subb)%>%data.frame()
  # Create a scatterplot
  plot <- ggplot(data = data_subb, aes(x = gx, y = gy) )+
    geom_point(size=4) +
    labs(title = paste("Plot", i), x = "are", y = "species")
  
  # Print the plot
  print(plot)
}
# to see the number of cores within each plot


#
a4=sample_data(a1)
a4=unique(a4$plotIDM)

core_number=numeric()
for (i in 1:length(a4))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(a1, plotIDM==a4[i])
      
  
  core_number[i]=dim(sample_data(data_sub))[1]
}


## for the model fitting
k=list()
p=unique(richness_scale_log$plotID)

for(i in 1:length(p))
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  d=subset(richness_scale,plotID==p[i])[,2:3]
  mf <- sar_average(data = d, grid_start = "none")
  k[[i]]=summary(mf)
}


kp=numeric()
p=unique(richness_scale$plotID)

for(i in 1:length(p))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  d=subset(richness_scale,plotID==p[i])[,2:3]
  kp[i] <- dim(d)[1]
  
}
