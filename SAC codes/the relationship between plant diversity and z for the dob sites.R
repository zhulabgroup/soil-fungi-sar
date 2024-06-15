# the plant diversity for the dob sites
library(tidyverse)

pd=read.csv("veg.OTU.Table.all.plot.csv",sep=",",header=T)
pd[,2:70][pd[,2:70]>0]=1

pd=cbind(pd$plot,apply(pd[,2:70],1,sum))%>%data.frame()
names(pd)=c("plotID","rich")
head(z)
pd_dob=merge(z,pd,by="plotID")
pd_dob$rich=as.numeric(pd_dob$rich)

# create a plot
ggplot(data=pd_dob,aes(x=rich,y=z))+
  geom_point(data=pd_dob,aes(x=rich,y=z),size=2,alpha=0.5)+
  geom_smooth(method="lm",lty="dotted")+
  ylab(expression(italic(z)))+
  xlab("Plant richness")+
  theme(legend.key.size = unit(0.15, "inches"),legend.position = c(0.5,0.85), legend.text.align = 0, panel.background = element_blank(), panel.border = element_rect(fill=NA,size=1,color="black"),legend.text = element_text(size=10),legend.title = element_text(size=15),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.y  = element_text(size=15),plot.title=element_text(hjust=0.5,face="bold",size=18),axis.text.x  = element_text(size=15))
