## to determine the richness at the 30 by 30 square

##
## select the first point and to see how many cores are selected within

## the orig of the square
pp=sample_data(rare_all_assign)%>%data.frame()
pp=unique(pp$plotIDM)
ori=expand.grid(x=seq(0,10,0.5),y=seq(0,10,0.5))

# the maximum number of soil cores included in the 30 by 30 plot

max_number=numeric()
for(j in 1:length(pp))
  {
  cat("\r", paste(paste0(rep("*", round(j / 1, 0)), collapse = ""), j, collapse = "")) # informs the processing
  
  kk=numeric()
  for(i in 1:dim(ori)[1])
  {
    d=sample_data(rare_all_assign)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==pp[j])
    d$gx=as.numeric(d$gx)
    d$gy=as.numeric(d$gy)
    kk[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
  }
  
  max_number[j]=kk%>%max()
}

ggplot()+
  geom_point(data=d,aes(y=gy,x=gx))+
  xlim(0,40)+
  ylim(0,40)

  geom_vline(xintercept = 34.5)+
  geom_hline(yintercept =34.5 )+
  geom_vline(xintercept = 4.5)+
  geom_hline(yintercept =4.5 )

## for each plot, we can select the origin that lead to more number of cores included
  
  richness_plot=numeric()
  for(k in 1:length(pp))
  {
    cat("\r", paste(paste0(rep("*", round(k / 1, 0)), collapse = ""), k, collapse = "")) # informs the processing
    
    core_number=numeric()
    for(i in 1:dim(ori)[1])
    {
      d=sample_data(a1)%>%data.frame()%>%select(gx,gy,plotIDM)%>%filter(plotIDM==pp[k])
      
      core_number[i]=subset(d,gx>=ori$x[i]&gx<=ori$x[i]+30 &gy<=ori$y[i]+30&gy>=ori$y[i])%>%summarize(row_count = n())%>%as.numeric()
    }
    core_number_8=core_number>8
    core_number_8[!core_number_8]=NA
     set.seed(5204)
    richness=numeric()
  for(j in c(which(!is.na(core_number_8))))# only select the locations that over 9 cores can be sampled/81 cores
  {
   # informs the processing
     sub_samp=subset(d,gx>=ori$x[j]&gx<=ori$x[j]+30 &gy<=ori$y[j]+30&gy>=ori$y[j])
     sub_samp=sample(rownames(sub_samp),9)# select only nine cores
    kk=subset_samples(a1,sample_names(a1)%in%sub_samp)
    richness[j]=table(colSums(otu_table(kk))>0)["TRUE"]%>%as.numeric()
  } 
  richness_plot[k]=richness%>%mean(na.rm = TRUE)
  }
  
  
  ## get the richness
  
  richness_plot%>%data.frame()%>%mutate(plotid=pp)
  
  
  load("~/soil-sar/richness_subplot10_neon_standar.RData")
  load("~/soil-sar/richness_subplot20_neon_standar.RData")
  load("~/soil-sar/plot-sar-permutation/richness_subplot40_mean_neon_standar.RData") 
  
# cbind all the data sets
  
  richness_subplot10_neon_standar%>%mutate(area=rep(100,length.out=n))->richness_subplot10_neon_standar
  
  richness_subplot20_neon_standar%>%mutate(area=rep(400,length.out=n))->richness_subplot20_neon_standar
  
  richness_subplot30_neon_standar=richness_plot%>%data.frame()%>%mutate(plotid=pp)%>%mutate(area=rep(900,length.out=n))%>%rename_all(~paste0(c("mean_value","plotid","area")))%>%select(plotid,mean_value,area)
  
  richness_subplot40_neon_standar=richness_subplot40_mean_neon_standar%>%rename_all(~paste0(c("mean_value","plotid","area")))%>%select(plotid,mean_value,area)
  
  ## cbind the data across the scales
  
  richness_data_neon=bind_rows(richness_subplot10_neon_standar%>%select(plotid, mean_value, area),richness_subplot20_neon_standar%>%select(plotid, mean_value, area),richness_subplot30_neon_standar%>%select(plotid, mean_value, area),richness_subplot40_neon_standar%>%select(plotid, mean_value, area))%>%filter(!is.na(mean_value ))

  richness_data_neon%>%count(plotid)%>%filter(n>2)->richness_data_neon_sar# we have 156 plots
  
  save(richness_data_neon,file="richness_data_neon.RData")
  
  # for the estimation of the z and  c value
  
 
  a=unique(richness_data_neon_sar$plotid)
  para_estimate=data.frame(ncol=3,nrow=length(a))
  for (i in 1:length(a))
    {
    mod=lm(formula = log(mean_value) ~log(area), data =richness_data_neon%>%filter(plotid==a[i])) 
    para_estimate[i,1]=a[i]
    para_estimate[i,2]=coef(mod)[1]
    para_estimate[i,3]=coef(mod)[2]
    
  }
  para_estimate= para_estimate%>%rename_all(~paste0(c("plotID","logc","z")))
  
  richness_data_neon_sar%>%rename(plotID=plotid)%>%left_join(para_estimate,by="plotID")->tem
 
   ## good to thing that number of points does not affect the estimated z
  