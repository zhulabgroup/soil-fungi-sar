## the tif-formated file shoul be read in
devtools::install_github('seschaub/getSpei') 
require(getSpei)
require(ncdf4)
require(chron)
neon_dob <- readRDS("/.../.../phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="AH")
neon_dob <- subset_samples(neon_dob, get_variable(neon_dob, "horizon")!="OH")# the data only include the O and M soil horizon
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob<- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob<- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
# choose a 3000 reads as the fixed sampling depth
rare_all=rarefy_even_depth(neon_dob, rngseed=10,sample.size = 3000, replace = F)#764 samples were removed
save(rare_all,file="rare_all.Rdata")# save the rarefied data
##

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
plot_spei=sample_data(d)[,c("plotIDM","lon","lat")]
plot_spei=data.frame(plot_spei)
names(plot_spei)=c("site","longitude","latitude")
site_spei <- spec_spei(spei_files = c("spei01"), start_y = 2010, end_y = 2018,locations=plot_spei)# from 2010 to 2018, each year wiht 12-months observations.
plot_level_speci=site_spei[,c(1,8)]
names(plot_level_speci)=c("plotID","spei")
write.csv(plot_level_speci,"plot_level_speci.csv")
