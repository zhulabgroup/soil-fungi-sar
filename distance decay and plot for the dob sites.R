# look at the distance decay pattern for the dob site

dob <- subset_samples(rare_all, Project == "DoB")
dob <- subset_samples(dob, horizon %in% c("O", "M"))
d <- sample_data(dob)

# extract the site and core info from sample name
site_names <- unique(d$Site)
sample_names <- matrix(unlist(strsplit(d$geneticSampleID, "[.]")), ncol = 3, byrow = T)
colnames(sample_names) <- c("Site1", "Location", "Horizon")
d <- cbind(d, sample_names)
sample_data(dob) <- d

dob <- subset_taxa(dob, taxa_sums(dob) > 0)
dob <- subset_taxa(dob, sample_sums(dob) > 0)
hist((sample_sums(otu_table(dob))))
quantile((sample_sums(otu_table(dob))), 0.1)

beta_diversity_mat <- vegdist(otu_table(dob), "jaccard")# beta diversity
d <- sample_data(dob)

dist_mat <- distm(d[, 1:2], fun = distHaversine)
dist_mat <- as.dist(dist_mat)

dist_mat=matrix(dist_mat)
beta_diversity_mat=matrix(beta_diversity_mat)

dob_dis=cbind(beta_diversity_mat,dist_mat)
dob_dis=data.frame(dob_dis)
names(dob_dis)=c("beta","distance")

# add a column that defines different scales of distances
# very tired way....

type=data.frame(dob_dis$distance)# quite weird that these codes do work as i expected.
dob_dis=cbind(dob_dis,type)

d1=subset(dob_dis,distance<=5)
d2=subset(dob_dis,distance>=57)
d3=subset(dob_dis,distance>=5&distance<57)

d1=cbind(d1,rep("fine",dim(d1)[1]))
d2=cbind(d2,rep("region",dim(d2)[1]))
d3=cbind(d3,rep("within",dim(d3)[1]))

names(d1)[4]="type"
names(d2)[4]="type"
names(d3)[4]="type"

ggplot(data=dob_dis,aes(x=distance,y=beta,color=type))+
  geom_point()+
  geom_smooth(method="lm")

# across different scales:within 5m, within 57 m and beyond 57m
# the scales were defined as "fine","within","region" scales

d5=rbind(d1,d1,d3,d1,d3,d2)
d6=rep(c("fine","within","region"),times=c(822,822+8617,822+8617+402339))
dob_dis_spa=cbind(d5,d6)

ggplot(data=dob_dis_spa,aes(x=log(distance+1),y=beta,color=d6))+geom_point(color="black",alpha=0.1)+
  geom_smooth(method="lm")+
  scale_color_manual("hehe",breaks=c("fine","region","within"),labels=c("<5 m","40x40 m","regional"),values=c("green","red","purple"))+
  theme(legend.position = "bottom", text = element_text(size=18), plot.title = element_text(size=15,hjust=0.5),axis.text.y = element_text(hjust = 0),axis.text.x = element_text(hjust = 1),axis.title.y = element_text(size=18),axis.title.x = element_text(size=18),axis.ticks.x = element_blank(),panel.background=element_rect(fill="NA"),panel.border = element_rect(color = "black", size = 1.5, fill = NA))
