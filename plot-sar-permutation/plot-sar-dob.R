## determine the plot-level SAR for the dob sites

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

dob= subset_samples(d,Project=="DoB")# the unique plotID, we have 476 plots

data_sample <- sample_data(dob)%>%data.frame()

dob_plot=unique(data_sample$plotIDM)

d3=subset_samples(dob,plotIDM=="MT1")%>%sample_data()%>%data.frame()

1.# to assign a code to each soil core based on their locations within the 40 by 40 plots
# need to remove the site as some sites include codes such as 0 and 5, which would otherwise disable the function
# for the AB4 plot, no 40-m data was available
# for the CO1 plot, only data at the 0 and the 40 m were available

sample_ID=data_sample$geneticSampleID

sample_ID=substr(sample_ID,5,10)

assign_ID=numeric()
for (i in 1:908)
  {
  sub_sample=sample_ID[i]
  
  if(str_detect(sub_sample,"00"))
    {
    assign_ID[i]=0
  }
  else if(str_detect(sub_sample,"A05|B05|C05|A5|B5|C5"))
  {
    assign_ID[i]=5
  }
  else if(str_detect(sub_sample,"10"))
  {
    assign_ID[i]=10
  }
  else if(str_detect(sub_sample,"20"))
  {
    assign_ID[i]=20
  }
  else if(str_detect(sub_sample,"15"))
  {
    assign_ID[i]=15
  }
  else 
  {
    assign_ID[i]=40
  }
  }

# add the soil core ID to the phyloseq object

assign_ID=data.frame(assign_ID)
row.names(assign_ID)=row.names(sample_data(dob))
assign_ID=sample_data(assign_ID)

dob=subset_samples(dob,select=-assign_ID)#delete the initial one

dob=merge_phyloseq(dob,assign_ID)

2 # get the richness of different spatial scales of core, 5x5, 10x10,20x20 and 40 x40

richness_dob=matrix(nrow=length(dob_plot),ncol=5)
# the fifth core do not have data for the 0.00 code
for (i in c(1:4,6:38,40:43))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    d0=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID==0)# get the first soil core 
    d5=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5")))# get the first soil core 
    d10=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5","10")))# get the first soil core 
    d20=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5","10","20")))# get the first soil core 
    d40=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5","10","20","40")))# get the first soil core 
    richness_dob[i,1]=table(taxa_sums(d0)>0)["TRUE"]%>%as.numeric() 
    richness_dob[i,2]=table(taxa_sums(d5)>0)["TRUE"]%>%as.numeric() 
    richness_dob[i,3]=table(taxa_sums(d10)>0)["TRUE"]%>%as.numeric() 
    richness_dob[i,4]=table(taxa_sums(d20)>0)["TRUE"]%>%as.numeric() 
    richness_dob[i,5]=table(taxa_sums(d40)>0)["TRUE"]%>%as.numeric() 
  }
  for (i in c(5))
# when i =5, there is no data for the 0.00 core, so we get the mean value of all the cores for the core
{
  d0=subset_samples(dob,plotIDM==dob_plot[i])
  d0=estimate_richness(d0,measures = "Observed")
  
  d5=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID==5)# get the first soil core 
  
  d10=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("5","10")))# get the first soil core 
  d20=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("5","10","20")))# get the first soil core 
  d40=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("5","10","20","40")))# get the first soil core 
  
  richness_dob[i,1]=mean(d0$Observed)
  richness_dob[i,2]=table(taxa_sums(d5)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,3]=table(taxa_sums(d10)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,4]=table(taxa_sums(d20)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,5]=table(taxa_sums(d40)>0)["TRUE"]%>%as.numeric() 
  }

#excluded from the analysis as the maximal size is 20 by 20 m 
for (i in c(39))
  # when i =5, there is no data for the 0.00 core, so we get the mean value of all the cores for the core
{
  d0=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID==0)# get the first soil core 
  d5=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5")))# get the first soil core 
  d10=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5","10")))# get the first soil core 
  d20=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5","10","20")))# get the first soil core 
  d40=subset_samples(dob,plotIDM==dob_plot[i]&assign_ID%in%(c("0","5","10","20","40")))# get the first soil core 
  richness_dob[i,1]=table(taxa_sums(d0)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,2]=table(taxa_sums(d5)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,3]=table(taxa_sums(d10)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,4]=table(taxa_sums(d20)>0)["TRUE"]%>%as.numeric() 
  richness_dob[i,5]=table(taxa_sums(d40)>0)["TRUE"]%>%as.numeric() 
}

richness_dob=cbind(dob_plot,richness_dob)%>%data.frame()

3.# add the area to the data frame

richness_dob[,2:6]=lapply(richness_dob[,2:6], as.numeric)
richness_dob=reshape::melt(richness_dob)
area=rep(c(0.0036,25,100,400,1600),each=43)
richness_dob=cbind(richness_dob,area)

# to see the relatiship between area and richness


for (i in c(1:38,40:43))
{
  data_subb <- subset(richness_dob,dob_plot==dob_plot [i])
  
  fit <-sar_loga(data_subb[,4:3] )
  plot=plot(fit)
  
  # Print the plot
  print(plot)
  
}
#some point are abnormal 
#42
#41
#40
#27
#21
#15
#9
# get the z and p value for each site

zva=numeric()
pva=numeric()
for (i in 1:length(dob_plot))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_subb <- subset(richness_dob,dob_plot==dob_plot [i])
  
  data_subb$value=log(data_subb$value)
  data_subb$area=log(data_subb$area)
  
  fit <-lm(data_subb[2:dim(data_subb)[1],c(3,4)] )%>%summary()
  
  zva[i]=fit$coefficients[2,1]
  pva[i]=fit$coefficients[2,4]
}

op=cbind(dob_plot,zva)%>%data.frame()

op$zva=as.numeric(op$zva)

names(op)[1]="plotID"

plot_level_zc_nest_mean$meanz=as.numeric(plot_level_zc_nest_mean$meanz)
plot_level_zc_nest_mean$meanc=as.numeric(plot_level_zc_nest_mean$meanc)

dk=merge(op,plot_level_zc_nest_mean,by="plotID")
# to see how it is related to the random approach


df=model_data[,c(1,3)]

head(df)

df=aggregate(z~plotID,data=df,FUN=mean)

df=merge(op,df,by="plotID")

df=merge(plot_level_zc_nest_mean,df,by="plotID")

ggplot()+geom_point(data=dk,aes(x=lon,y=lat),size=1)+
  geom_text(data=dk,aes(x=lon,y=lat),size=3,label=substr(dk$geneticSampleID,1,7))
