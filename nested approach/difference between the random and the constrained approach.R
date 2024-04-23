##

# distance pattern



## the distance-decay pattern will use the full data generated in step 1
head(sample_data(rare_all))

d <- sample_data(rare_all)
table(d$Project) # the first 908 rows are dob sites with the remaining 5470 being NEON sites
d1 <- data.frame(d[, c("geneticSampleID", "Site")])
plotID <- substr(d1$geneticSampleID, 1, 8)
d1 <- cbind(d1, plotID)
iddob <- d1$Site[1:908] # a site correspondes to a plot
idneon <- d1$plotID[909:6378] # an unique plotID corresponds to a plot
plotIDM <- data.frame(c(iddob, idneon))
names(plotIDM) <- "plotIDM" # the plot id used for the SAR
row.names(plotIDM) <- row.names(d)
plotIDM <- sample_data(plotIDM)
d <- merge_phyloseq(rare_all, plotIDM) # merge the new plotid with the initial data
# select an unique plot and build a SAR within the plot
a1 <- sample_data(d) # the unique plotID, we have 476 plots
a1 <- unique(a1$plotIDM)
rare_all <- d


# for the neon data
sub_neon <- subset_samples(rare_all, get_variable(rare_all, "Project") == "NEON")
sub_neon <- subset_samples(sub_neon, get_variable(sub_neon, "horizon") != "AH")
sub_neon <- subset_samples(sub_neon, get_variable(sub_neon, "horizon") != "OH")

# the location within the 40 by 40 plot
loca <- data.frame(sample_data(sub_neon))["geneticSampleID"] # need to be
loca <- substr(loca$geneticSampleID, 11, 20)
loca <- data.frame(loca)
d <- strsplit(loca$loca, "-")

a <- matrix(nrow = dim(loca)[1], ncol = 2) # the location of each core within the 40 by 40 plot
for (i in 1:dim(loca)[1])
{
  a[i, ] <- as.numeric(d[[i]][2:3])
}

# need to cbind the location data with the plotIDM

a <- data.frame(a)
names(a) <- c("gx", "gy")
row.names(a) <- row.names(sample_data(sub_neon))
a <- sample_data(a)
sub_neon <- merge_phyloseq(sub_neon, a) # adding the location data to the full data.


## add a distance to the data

op=sample_data(sub_neon)%>%data.frame()

op=op[,c("gx","gy")]

# distance to the original point

op=mutate(op,distance=sqrt(gx^2+gy^2))



row.names(op)=row.names(sample_data(sub_neon))

op=sample_data(op)

sub_neon=merge_phyloseq(sub_neon,op)

# select the first plot

ot=list()
a1=sample_data(sub_neon)
a1=unique(a1$plotIDM)
for (i in 1:length(a1)){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing

data_sub <- subset_samples(sub_neon, plotIDM == a1[i]) # all the samples in a given plotID
dk=sample_data(data_sub)
ot[[i]]=otu_table(data_sub)[rownames(dk[order(dk$distance),]),]
}

# get the richness for each of the ten plot




set.seed(1005)

df_ran=list()

for (i in 1:10){
  dd1 <- specaccum(ot[[i]], method = 'random')

  df_ran[[i]]=dd1$richness

}

df_col=list()
for (i in 1:10){
  dd2 <- specaccum(ot[[i]], method = 'collector')
  
  df_col[[i]]=dd2$richness
}


## combine the data 
com_ran=list()
for (i in 1:10){
  
 two_app= cbind(df_ran[[i]],df_col[[i]])

 two_appp=two_app%>%melt()

com_ran[[i]]=mutate(two_appp,area=rep(c(1:(dim(two_appp)[1]/2)),times=2))

}

###creat the plots

p1=ggplot()+
  geom_point(data=com_ran[[1]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")


p2=ggplot()+
  geom_point(data=com_ran[[1]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[1]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")




p3=ggplot()+
  geom_point(data=com_ran[[2]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p4=ggplot()+
  geom_point(data=com_ran[[2]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[2]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p5=ggplot()+
  geom_point(data=com_ran[[3]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p6=ggplot()+
  geom_point(data=com_ran[[3]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[3]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p7=ggplot()+
  geom_point(data=com_ran[[4]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p8=ggplot()+
  geom_point(data=com_ran[[4]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[4]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p9=ggplot()+
  geom_point(data=com_ran[[5]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p10=ggplot()+
  geom_point(data=com_ran[[5]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[5]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p11=ggplot()+
  geom_point(data=com_ran[[6]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p12=ggplot()+
  geom_point(data=com_ran[[6]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[6]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p13=ggplot()+
  geom_point(data=com_ran[[7]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p14=ggplot()+
  geom_point(data=com_ran[[7]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[7]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p15=ggplot()+
  geom_point(data=com_ran[[8]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p16=ggplot()+
  geom_point(data=com_ran[[8]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[8]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")



p17=ggplot()+
  geom_point(data=com_ran[[1]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p18=ggplot()+
  geom_point(data=com_ran[[9]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[9]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")



p19=ggplot()+
  geom_point(data=com_ran[[10]],aes(x=area,y=value,color=as.factor(X2)))+
  
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")


p20=ggplot()+
  geom_point(data=com_ran[[10]],aes(x=log(area),y=log(value),color=as.factor(X2)))+
  geom_smooth(data=com_ran[[10]],aes(x=log(area),y=log(value),color=as.factor(X2)),method="lm")+
  scale_color_manual("",breaks=c(1,2),labels=c("Random","constrained"),values = c("blue","red"))+
  ylab("Log(Richness)")+
  xlab("log(Area) m^2")+
  guides(color="none")













df_ran_com=df_ran[[1]]
for (i in 2:10){
  
  df_ran_com=rbind(df_ran_com,df_ran[[i]])
}





df_ran=matrix(nrow=22,ncol=10)

for (i in 1:22){
  dd1 <- specaccum(ot[[i]], method = 'random')
  
  df_ran[i]=as.matrix(dd1$richness)
  
}






# get the richness of 

estimate_richness(ot)

sp <- specaccum(ot, method = 'random')

plot(sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'lightblue')

boxplot(sp, col = 'yellow', add = TRUE, pch = '+')

df=data.frame(sp$richness,1:22)
df1=data.frame(sp1$richness,1:22)
names(df1)=names(df)

df3=rbind(df,df1)

df3$X1.22=as.numeric(df3$X1.22)

ggplot()+
geom_point(data=df3,aes(x=X1.22 ,y=sp.richness,color=type))

ggplot()+
  geom_point(data=df3,aes(x=log(X1.22) ,y=log(sp.richness),color=type))+
  
  geom_smooth(data=df3,aes(x=log(X1.22) ,y=log(sp.richness),color=type),method="lm")+
  xlab("number of soil cores")+
  ylab("Species number")
  



df3=mutate(df3,type=rep(c("random","constrained"),each=22))


sp1 <- specaccum(ot, method = 'collector')

plot(sp1, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'lightblue')

boxplot(sp1, col = 'yellow', add = TRUE, pch = '+')

