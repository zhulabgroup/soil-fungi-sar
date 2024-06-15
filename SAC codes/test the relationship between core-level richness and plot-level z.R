# to test how the richness per core affects the plot-level z values
# take the 'soil saprotrophic' guild for example
1.# for the "soil_saprotroph"

m=subset_taxa(rare_all,ta2=="soil_saprotroph")# select a guild
m<- subset_samples(m, !is.na(lon) & !is.na(lat))
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

rich<- numeric()# the mean richness per core for a site(alpha diversity)

for (i in 1:length(a1))
{
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  data_sub <- subset_samples(d, plotIDM==a1[i])
  d1=estimate_richness(data_sub,measures="Observed")
  rich[i]=mean(d1$Observed)
}

rich_ecm_core=data.frame(cbind(a1,rich))

ecm_z_ranall=data.frame(ecm_z_ranall)
k=apply(ecm_z_ranall[,3:32],1,mean)
k=cbind(ecm_z_ranall["a1"],k)
head(rich_ecm_core)
k=merge(rich_ecm_core,k,by="a1") # do not significantly related  
summary(lm(rich~k,data=subset(k,k<10)))
head(k)