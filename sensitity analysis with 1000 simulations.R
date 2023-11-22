library(ggcorrplot)
library(lme4)
library(lmerTest)
library(tidyverse)
library(cowplot)

# sensitive analysis for the climate and soil model
# testing how the results change when different number of repetition were included
## for the z values
1.# bind the 1000 simulated c and the z values
a=unique(all_z_ran1000$a1)

d=list()
for (i in 1:length(a)){
  d1=subset(all_z_ran1000,a1==a[i])
  d[[i]]=t(d1[,3:1002])
}

z_sen=d[[1]]
for(i in 2:length(a)){
  z_sen=rbind(z_sen,d[[i]])
}

## for the c values
d=list()
for (i in 1:515){
  
  d1=subset(all_c_ran1000,a1==a[i])
  d[[i]]=t(d1[,3:1002])
}

c_sen=d[[1]]
for(i in 2:515){
  c_sen=rbind(c_sen,d[[i]])
}

# sensitive analysis data
plotID=rep(a,each=1000)
data_sen=data.frame(cbind(plotID,c_sen,z_sen))
data_sen$V2=as.numeric(data_sen$V2)
data_sen$V3=as.numeric(data_sen$V3)
names(data_sen)=c("plotID","logc","z")

2.# to merge the data with the environmental variables


data_sen_mod=merge(data_sen,unique(model_data[,c(2:3,6:28)]),by="plotID")# better to use the unique climate variables
data_sen_mod=subset(data_sen_mod,siteIDD != "GUAN"&z<10)# excluded some plots that have less than 3 cores
# select the data
range01 <- function(x) ## to
{
  return((x - min(x)) / (max(x) - min(x)))
}

###
cc=2.71828^data_sen_mod$logc# first transformed the logc and then standardized the data
kk=cbind(data_sen_mod,cc)

data_sen_mod2=kk[,c(1,2,3,5:17,19,27,28)]# select only climate,fungal diversity and soil variables 

data_sen_mod2[,c(4:19)]=apply(data_sen_mod2[,c(4:19)],2,range01)%>% data.frame()
write.csv(data_sen_mod2,"data_sen_mod2.csv")

#. construct the model with plotID included as the only random effect model

mod <- lmer(z ~ logc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = data_sen_mod2)


pv=list()
ef=list()
for (j in 20:100)
  {
  cat('\r',paste(paste0(rep("*", round(j/ 1, 0)), collapse = ''), j, collapse = ''))# informs the processing
  a=unique(data_sen_mod2$plotID)
  ddf=list()
  for (i in 1:length(a))
  {
    a1=subset(data_sen_mod2,plotID==a[i])
    
    ddf[[i]]=sample_n(a1,j,replace = FALSE)#we selected 20 rows for one site 
  }
  
  ddf1=ddf[[1]]
  for(i in 2:length(a)){
    ddf1=rbind(ddf1,ddf[[i]])
  }
  
  mod <- lmer(z ~ cc + organicCPercent + ph + nitrogen + sand + bio2 + bio8 + bio18 + bio4 + bio12 + bio15 + spei + funrich + bio1 + (1 | plotID), data = ddf1)
  op=summary(mod)
  pv[[j]]=op$coefficients[1:15,5]# the p-value for each selected variables
  ef[[j]]=op$coefficients[1:15,1]# the effect size for each selected variables
}


pvv=pv[[20]]
for(j in 21:100){
  pvv=rbind(pvv,pv[[j]])
}
pvv=data.frame(pvv)
q=tidyr::gather(pvv)
boxplot(value~key,data=q)# to see the 

sen_pv=q
write.csv(sen_pv,"sen_pv.csv")

###
p1=ggboxplot(sen_pv, x = "key", y = "value", outlier.colour = "gray") +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  xlab("") +
  theme(legend.position = c(0.51, 0.85), legend.text = element_text(size = 14), text = element_text(size = 15), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90)) +
  scale_x_discrete(breaks=unique(sen_pv$key),label=c("Intercept","c","SoilC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","MAP","Pre.seas","SPEI","Fun.rich","MAT"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue")+
  ylab(expression(italic(P)*" value"))

###

eff=ef[[20]]
for(j in 21:100){
  eff=rbind(eff,ef[[j]])
}
eff=data.frame(eff)
q=tidyr::gather(eff)
boxplot(value~key,data=q)# to see the 

sen_eff=q
write.csv(sen_eff,"sen_eff.csv")

##

p2=ggboxplot(sen_eff, x = "key", y = "value", outlier.colour = "gray") +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  #geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  xlab("") +
  theme(legend.position = c(0.51, 0.85), legend.text = element_text(size = 14), text = element_text(size = 15), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90)) +
  scale_x_discrete(breaks=unique(sen_eff$key),label=c("Intercept","c","SoilC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","MAP","Pre.seas","SPEI","Fun.rich","MAT"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue")+
  ylab("Effect size")
# 
p1=ggplotGrob(p1)
p2=ggplotGrob(p2)
p1$widths=p2$widths
plot_grid(p2,p1,ncol=1)

#
ggboxplot(data_sen_mod2, x = "plotID", y = "z", outlier.colour = "mediumpurple",size=0.15,alpha=0.5) +
  #geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  xlab("plotID") +
  theme(legend.position = c(0.51, 0.85), legend.text = element_text(size = 14), text = element_text(size = 15), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.ticks.x = element_blank(),axis.text.x = element_text(angle=90)) +
  #scale_x_discrete(breaks=unique(sen_eff$key),label=c("Intercept","c","SoilC","pH","SoilN","Sand","Bio2","Bio8","Bio18","Tem.seas.","MAP","Pre.seas","SPEI","Fun.rich","MAT"))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "blue")+
  ylab(expression(italic(Z)*" value"))
  