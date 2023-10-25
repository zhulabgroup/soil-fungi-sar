# look at the relationship between fungal diversity and environmental variables


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
plot_cor=sample_data(d)[,c("lon","lat","plotIDM")]#plot coordinates
# select an unique plot and build a SAR within the plot
a1= sample_data(d)# the unique plotID, we have 476 plots
a1=unique(a1$plotIDM)

# determine the plot level richness

rich=numeric()
for (i in 1:length(a1)){
  #cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(rare_all, plotIDM==a1[i])
  
  otu=otu_table(neon_sub)
  rich[i]=table(apply(otu,2,sum)>0)[[2]]# richness within a specific plot for a site
  
}





rich_plot=data.frame(a1,rich)
siteid=substr(rich_plot$a1,1,4)

rich_plot=cbind(rich_plot,siteid)
names(rich_plot)[1]="plotIDM"

plot_cor=data.frame(plot_cor)

rich_cor=unique(merge(rich_plot,plot_cor,by="plotIDM"))




