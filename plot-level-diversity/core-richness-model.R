# core level richness



variables=c("grid","core_rich", "lon","lat", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18"  )


data_core_rich=subset(data_plot_rich,siteIDD!="GUAN")

mod_core_rich=aggregate(data_plot_rich[,variables],by=list(data_plot_rich$grid),FUN=mean,na.rm=TRUE)# aggregated at the 10 by 10 scale

# data preparation

powerTransform(mod_core_rich$core_rich)

mod_core_rich$richbc<-bcPower(mod_core_rich$core_rich, 0.9882336 )

# remove the NAs
#standardize the data

mod_core_rich[,6:15]=apply(mod_core_rich[,6:15],2,range01)

train_core_rich=mod_core_rich[,c("richbc","lon","lat")]

pre_core_rich=mod_core_rich[,c("richbc", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18")]


newdata=newdata[complete.cases(newdata),]
newloca=newdata[,c(1,2)]
newpre=newdata[,c(3:12)]
newpre=newpre[complete.cases(newpre),]
newpre=apply(newpre,2,range01)

RESPONSE<-train_core_rich$richbc
mod_predictors<-pre_core_rich[2:11]

coordinates(train_core_rich) = ~lon+lat
coordinates(newloca) = ~lon+lat


mc <- makeCluster(detectCores())
registerDoParallel(mc)

myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1808)

GLM_core<-train(mod_predictors,
                RESPONSE,
                method = "glm",
        
                trControl=myControl,
                preProc=c('center', 'scale'))
print(GLM_core)


train_core_rich$residuals.glm<-resid(GLM_core)



v.glm<-variogram(residuals.glm~ 1, data = train_core_rich,cutoff=35, width=35/15)

v.glm<-variogram(residuals.glm~ 1, data = train_core_rich)

plot(v.glm)

m.glm<-vgm(1200,"Exp",11,3000)# based on estimation of the semivariance

m.f.glm<-fit.variogram(v.glm, m.glm)
m.f.glm


plot(v.glm, pl=F, 
     model=m.f.glm,
     col="black", 
     cex=0.9, 
     lwd=0.5,
     lty=1,
     pch=19,
     main="Variogram and Fitted Model\n Residuals of GLM model",
     xlab="Distance (m)",
     ylab="Semivariance")

###



newloca$GLM_core <- predict(GLM_core, newpre)# the x and y disappear

SK.GLM_core<-krige(residuals.glm~ 1, 
                   loc=train_core_rich,        # Data frame
                   newdata=newloca,     # Prediction location
                   model = m.f.glm,     # fitted varigram model
                   beta = 0)    


newloca$SK.GLM_core<-SK.GLM_core$var1.pred
# Add RF predicted + SK predicted residuals

newloca$RK.GLM.bc<-(newloca$GLM_core+newloca$SK.GLM_core)

k1<-1/0.9882336  

newloca$RK.GLM <-((newloca$RK.GLM.bc *0.9882336 +1)^k1)

summary(newloca)


attributes_df <- data.frame(newloca$RK.GLM)

coordinates_df <- data.frame(coordinates(newloca))

pred_core=cbind(coordinates_df,attributes_df)
names(pred_core)[3]="rich"


p1=ggplot(pred_core) +
  geom_point(data=pred_core,pch=15,aes(x=lon, y=lat,color=rich), size=0.275)+
  scale_color_gradient(expression(Richness),low = "blue", high = "yellow")+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Predicted core-level richness")+
  ylab("")+
  xlim(-175,-42)+
  geom_point(data=plot_diversity,pch=21,aes(x=lon,y=lat,fill=project),alpha=1)+
  scale_fill_manual("Project",breaks=c("dob","neon"),label=c("DoB","NEON"),values=c("red","seagreen"))

# if we used poss-ion regression
# the effects of individual variables on the richness
k=summary(GLM_core)

effect_core=k$coefficients[2:11,]%>%data.frame()

p11=ggplot()+
  
  geom_segment(data=effect_core,aes(x=Estimate/10-1.96*Std..Error/10 ,y=1:10,yend=1:10,xend=Estimate/10+1.96*Std..Error/10),color="gray")+
  geom_point(data=effect_core,aes(x=Estimate/10,y=1:10),color="blue",size=4)+
  geom_vline(xintercept = 0,linetype="dashed")+
  theme(legend.position =c(0.2,0.35), 
        legend.key.size = unit(0.15, "inches"),
        guides(color = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_text(size=8),
        text = element_text(size = 18), 
        legend.text = element_text(size=8),
        plot.title = element_text(size = 15, hjust = 0.5), 
        axis.text.y = element_text(hjust = 0), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "NA"), 
        panel.border = element_rect(color = "black", size = 1.5, fill = NA))+
  xlab("Effect size")+
  ylab("")+
  scale_y_continuous(breaks = 1:10, label = c("SoilOC", "pH", "SoilN", "Sand", "MAT", "MDR", "Tem.seas.","MAP","Pre.seas.","PWQ")) +
  xlab("Effect size and 95% CI") +
  ggtitle("Core-level richness (33.3%)") 
  

