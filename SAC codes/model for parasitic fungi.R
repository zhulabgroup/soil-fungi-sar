# model for parasitic fungi

# determine the mean richness within each core for parasitic fungi

1.# prepare the data
ft <- read.csv("FungalTraits_1.2_ver_16Dec_2020.csv", sep = ",", header = T)
ft <- ft[, c("GENUS", "primary_lifestyle")]
# get the genus name of all the taxa in the full data set
d <- tax_table(rare_all)[, 6]
d <- data.frame(d)

names(ft) <- c("genus", "guild")
# assign guild types to each taxa
d <- left_join(d, ft, by = "genus")
# change it to the tax-table format
d <- tax_table(d)
row.names(d) <- row.names(tax_table(rare_all))
# adding the guild column to the existing phyloseq object
rare_all <- merge_phyloseq(rare_all, d)

2.# the mean richness at the core level for parasitic fungi

pa <- data.frame(c("animal_parasite", "lichen_parasite", "protistan_parasite", "algal_parasite"))
names(pa) <- "pa"
a <- subset_taxa(rare_all, ta2%in%pa$pa) # must first select the guild and then check the sample and taxa sums
a <- subset_samples(a, !is.na(lon) & !is.na(lat))
a <- subset_taxa(a, taxa_sums(a) > 0)
a <- subset_samples(a, sample_sums(a) > 0)
a <- subset_samples(a, sample_sums(a) > 0)
rich_para=estimate_richness(a, measures="Observed")
k=data.frame(sample_data(a))["plotIDM"]
k=cbind(k,rich_para)
rich_papra_mean=aggregate(Observed~plotIDM,data=k,FUN=mean)
names(rich_papra_mean)[1]="plotID"

3.#get the z value for the parasitic fungi
# get the c value
a=list()
for (i in 1:dim(para_c_ranall)[1])
{
  d=data.frame(t(para_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(para_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(para_c_ranall$a1,each=30)
para_c_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values

# get the z value
a=list()
for (i in 1:dim(para_z_ranall)[1])
{
  d=data.frame(t(para_z_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(para_z_ranall)[1])
{
  b=rbind(b,a[[i]])
}

para_z_ranall_30=cbind(plotID,b)# each plot has 30 estimated z values
names(para_z_ranall_30)[2]="z"
para_model=cbind(para_c_ranall_30,para_z_ranall_30["z"])
names(para_model)[2]="logc"
para_model <- cbind(para_model, c = 2.71828^para_model$logc)
para_model=merge(para_model,model_var,by="plotID")

para_model=subset(para_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&rich>0)# only 104 plots from 33 sites

# the data that include core level richness
para_model_rich=subset(para_model,siteIDD!="GUAN"&z<10&richness>0)
para_model_rich=merge(para_model_rich,rich_papra_mean,by="plotID")

#standardized data
para_model[,c(2,5:27)]=apply(para_model[,c(2,5:27)],2,range01)
para_model_rich[,c(6:21,29)]=apply(para_model_rich[,c(6:21,29)],2,range01)

# colinearity
ggcorrplot(cor(para_model[,c(2,3,5:27)]), hc.order = TRUE, type = "lower", lab = TRUE)#

#bold was related with soil OC and was excluded,cec was related to soil OC and was removed
#coarse (removed) was related with fine and was removed
#cec(removed) and soil OC
#root cn and rootn(removed)
#bold (removed)and soil OC
#bio1 and bio4(removed)
3
# build a model for the para guild,with 104 plots included, fungal rich was most fluencial

mod <- lmer(z ~ Observed + funrich+organicCPercent + ph + nitrogen + sand +bio1+ bio2 + bio8 + bio12 + bio15 +bio18+ spei + richness  + (1 | siteIDD / plotID), data = para_model_rich)
effect_para=summary(mod)
effect_para=effect_para$coefficients
effect_para=data.frame(effect_para)[2:dim(effect_para)[1],]

p8=ggplot()+
  geom_point(data=effect_para,aes(x= Estimate,y=1:dim(effect_para)[1]),
             color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),size=3)+
  geom_segment(data=effect_para,size=0.8,color=rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple")),aes(x=effect_para$Estimate-1.96*effect_para$Std..Error,y=1:dim(effect_para)[1],xend=effect_para$Estimate+1.96*effect_para$Std..Error,yend=1:dim(effect_para)[1]))+
  geom_vline(xintercept = 0,color="red",linetype="dashed")+
  scale_y_continuous(breaks=1:14,labels = rev(c("Pla.rich", "Spei","Pre.WQ","Pre.seas.","MAP","MTWQ","MDR","MAT","Sand","SoilN","pH","SoilC","Plot.rich","Core.rich")))+
  theme(axis.text.y = element_text(colour = rev(c("seagreen1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","royalblue1","peru","peru","peru","peru","purple","purple"))))+
  theme(legend.key.size = unit(0.18, "inches"),   
        legend.position = c(0.4,0.85), 
        legend.text.align = 0, panel.background = element_blank(), 
        panel.border = element_rect(fill=NA,size=1,color="black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(size=15,color="black"),
        plot.title=element_text(hjust=0.5,face="bold",size=18),
        axis.text.x  = element_text(size=15))+
  xlim(-0.8,0.5)+
  ylab("")+
  annotate("text",x=0.35,y=14,label="**",size=8)+
  annotate("text",x=-0.5,y=13,label="**",size=8)+
  annotate("text",x=0.35,y=11,label="***",size=8)+
  annotate("text",x=0.35,y=10,label="",size=8)+
  annotate("text",x=-0.30,y=1,label="*",size=8)+
  ggtitle("Parasitic (N=409)")


