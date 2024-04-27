####
d=sample_data(rare_all)
a1=unique(d$plotIDM)

# core-level richness for all the plots: one value for each plot
core_rich=data.frame(nrow=length(a1),ncol=2)
for (i in 1:length(a1))
  {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all, plotIDM==a1[i])
  
  rich=estimate_richness(data_sub,measures="Observed")%>%data.frame()
  core_rich[i,1]=a1[i]
  core_rich[i,2]=mean(rich$Observed)
}

# plot-level estimated richness


plot_rich=data.frame(nrow=length(a1),ncol=3)
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all, plotIDM==a1[i])
  data_sub =transform_sample_counts(data_sub, function(x) ifelse(x>0, 1, 0))
  ot=t(otu_table(data_sub))%>%data.frame()
 
  n=dim(ot)[2]# the number of "sites" for each plot
  
  if(n>2)
    {
    ot1=rowSums(ot)
    
    out3 <- iNEXT(c(n,ot1), q=0, datatype="incidence_freq", size=round(seq(10, 200, length.out=10)), se=FALSE)
    
    plot_rich[i,1]=a1[i]
    plot_rich[i,2]=out3$iNextEst$size_based[13,5]
    plot_rich[i,3]=out3$iNextEst$size_based[13,2]
    
  }
   
  else{

    plot_rich[i,1]=a1[i]
    plot_rich[i,2]=NA
    plot_rich[i,3]=NA
}
  
}




cor_number=numeric()
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all, plotIDM==a1[i])
  
  cor_number[i]=dim(sample_data(data_sub))[1]
  
}


### plot-level beta diversity

# for the neon data


# the location within the 40 by 40 plot
# for the dob data, the location of the points need to be transformed

loca <- data.frame(sample_data(rare_all))["geneticSampleID"] # need to be



assign_ID=data.frame(ncol=2,nrow=6378)

for (i in 1:6378)
  {
  
  sub_sample=loca[i,]

if (nchar(sub_sample)>15)
  {
  subb <- substr(sub_sample, 11, 20)
  subb <- data.frame(subb)
  d <- strsplit(subb$subb, "-")
  assign_ID[i,]=d[[1]][2:3]
}
  
else {
  if(str_detect(sub_sample,"00"))
  {
    assign_ID[i,]=c(0,0)
  }
  else if(str_detect(sub_sample,"A05|A5"))
  {
    assign_ID[i,]=c(0,5)
  }
  else if(str_detect(sub_sample,"A10"))
  {
    assign_ID[i,]=c(0,10)
  }
  else if(str_detect(sub_sample,"A20"))
  {
    assign_ID[i,]=c(0,20)
  }
  else if(str_detect(sub_sample,"A40"))
  {
    assign_ID[i,]=c(0,40)
  }
  
  else if(str_detect(sub_sample,"B05|B5"))
  {
    assign_ID[i,]=c(5,5)
  }
  else if(str_detect(sub_sample,"B10"))
  {
    assign_ID[i,]=c(10,10)
  }
  else if(str_detect(sub_sample,"B15"))
  {
    assign_ID[i,]=c(15,15)
  }
  else if(str_detect(sub_sample,"B20"))
  {
    assign_ID[i,]=c(20,20)
  }
  else if(str_detect(sub_sample,"B40"))
  {
    assign_ID[i,]=c(40,40)
  }
  else if(str_detect(sub_sample,"C05|C5"))
  {
    assign_ID[i,]=c(5,0)
  }
  else if(str_detect(sub_sample,"C10"))
  {
    assign_ID[i,]=c(10,0)
  }
  else if(str_detect(sub_sample,"C15"))
  {
    assign_ID[i,]=c(15,0)
  }
  else if(str_detect(sub_sample,"C20"))
  {
    assign_ID[i,]=c(20,0)
  }
  else 
  {
    assign_ID[i,]=c(40,0)
  }
}
  
}

 








# need to cbind the location data with the plotIDM


names(assign_ID) <- c("gx", "gy")
assign_ID$gx=as.numeric(assign_ID$gx)
assign_ID$gy=as.numeric(assign_ID$gy)
row.names(assign_ID) <- row.names(sample_data(rare_all))
assign_ID <- sample_data(assign_ID)
rare_all <- merge_phyloseq(rare_all, assign_ID) # adding the location data to the full data.

a1 <- sample_data(rare_all) # the unique plotID, we have 476 plots
a1 <- unique(a1$plotIDM) # 472 plotID

pair <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all, plotIDM == a1[i]) # all the samples in a given plotID
  location <- sample_data(data_sub)[, c("gx", "gy")]
  dis_tance <- dist(location, diag = TRUE, upper = TRUE)
  pair[[i]] <- matrix(dis_tance)
}

#
pair1=pair[[1]]
{
  for (i in 2:515)
    pair1=rbind(pair1,pair[[i]])  
}



# the beta diversity between anytwo cores within a 40 by 40 plot

beta <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all, plotIDM == a1[i]) # all the samples in a given plotID
  com <- otu_table(data_sub)
  ddis <- vegdist(com, method = "jaccard")
  beta[[i]] <- matrix(ddis) # the spatial distance between two locations
  
}
  
plot_id <- list()
for (i in 1:length(a1)) {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  data_sub <- subset_samples(rare_all, plotIDM == a1[i]) # all the samples in a given plotID
  com <- otu_table(data_sub)
  ddis <- vegdist(com, method = "jaccard")
  plot_id[[i]] <- rep(a1[i],times=length(matrix(ddis)) ) 
                    
}






# combine the community distance and beta diversity

decay <- list()
for (i in 1:length(a1))
{
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) #
  decay[[i]] <- cbind(pair[[i]], beta[[i]],plot_id[[i]])
}

decay_com <- decay[[1]]
for (i in 2:length(a1)) {
  decay_com <- rbind(decay_com, decay[[i]])
}

decay_com <- data.frame(decay_com)
names(decay_com) <- c("distance", "beta","plotID")

decay_com[,1:2]=sapply(decay_com[,1:2], as.numeric)
# look at how the relationship changes as increasing distance

a1=unique(decay_com$plotID)
r <- matrix(ncol=2,nrow=length(a1))
for (i in 1:length(a1))
{
  data = subset(decay_com, plotID==a1[i])
  if(dim(data)[1]>2)
    {
    cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
    d <- summary(lm(beta ~ distance, data))
    r[i,1] <- d$coefficients[2, 1]
    r[i,2] <- d$coefficients[2, 4] 
  }
  else{
    r[i,1] <- NA
    r[i,2] <- NA
  }
}

# 75% of the plot dones not show a distance decay 
# with only one point

mm=table(decay_com$plotID)%>%data.frame()
mm=subset(mm,Freq>2)


which(!is.na(mm))
sel_row=which(!is.na(mm))%>%as.numeric()# selected rows

# does core number affect

names(mm)=c("plotID","core")

a2=merge(decay_com,mm,by="plotID")

a2=aggregate(a2[,3:4],by=list(a2$plotID),FUN=mean)


ggplot()+
geom_point(data=subset(decay_com,plotID%in%mm$Var1),aes(y=log(beta),x=log(distance),color=plotID))+

guides(color="none")+
geom_smooth(data=subset(decay_com,plotID%in%mm$Var1),aes(y=log(beta),x=log(distance),color=plotID),method="lm",se=FALSE)





a1=unique(decay_com$plotID)

r <- numeric()

for (i in 1:length(a1))
{
 
  d <-plot(beta ~ distance, data = subset(decay_com, plotID==a1[i]))
  
  print(d)
}


