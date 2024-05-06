## beta diversity


variables=c("grid","beta_mean", "lon","lat", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18"  )


data_beta_rich=subset(data_plot_rich,siteIDD!="GUAN")

mod_beta_rich=aggregate(data_plot_rich[,variables],by=list(data_plot_rich$grid),FUN=mean,na.rm=TRUE)# aggregated at the 10 by 10 scale

# data preparation

powerTransform(mod_beta_rich$beta_mean)

mod_beta_rich$richbc<-bcPower(mod_beta_rich$beta_mean, 7.904799)

# remove the NAs
#standardize the data

mod_beta_rich[,6:15]=apply(mod_beta_rich[,6:15],2,range01)

train_beta_rich=mod_beta_rich[,c("richbc","lon","lat")]

pre_beta_rich=mod_beta_rich[,c("richbc", "organicCPercent", "ph","nitrogen", "sand","bio1", "bio2","bio4" , "bio12" ,"bio15" , "bio18")]


newdata=newdata[complete.cases(newdata),]
newloca=newdata[,c(1,2)]
newpre=newdata[,c(3:12)]
newpre=newpre[complete.cases(newpre),]
newpre=apply(newpre,2,range01)

RESPONSE<-train_beta_rich$richbc
mod_predictors<-pre_beta_rich[2:11]

coordinates(train_beta_rich) = ~lon+lat
coordinates(newloca) = ~lon+lat


mc <- makeCluster(detectCores())
registerDoParallel(mc)

myControl <- trainControl(method="repeatedcv", 
                          number=10, 
                          repeats=5,
                          allowParallel = TRUE)


set.seed(1809)

GLM_beta<-train(mod_predictors,
                RESPONSE,
                method = "glm",
                
                trControl=myControl,
                preProc=c('center', 'scale'))
print(GLM_beta)


train_beta_rich$residuals.glm<-resid(GLM_beta)



v.glm<-variogram(residuals.glm~ 1, data = train_beta_rich,cutoff=35, width=35/15)

v.glm<-variogram(residuals.glm~ 1, data = train_beta_rich)

plot(v.glm)

m.glm<-vgm(0.00014,"Sph",12,0.00016)# based on estimation of the semivariance

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



newloca$GLM_beta<- predict(GLM_beta, newpre)# the x and y disappear

SK.GLM_beta<-krige(residuals.glm~ 1, 
                   loc=train_beta_rich,        # Data frame
                   newdata=newloca,     # Prediction location
                   model = m.f.glm,     # fitted varigram model
                   beta = 0)    


newloca$SK.GLM_beta<-SK.GLM_beta$var1.pred
# Add RF predicted + SK preedicted residuals

newloca$RK.GLM.bc<-(newloca$GLM_beta+newloca$SK.GLM_beta)

k1<-1/7.904799 

newloca$RK.GLM <-((newloca$RK.GLM.bc *7.904799 +1)^k1)

summary(newloca)


attributes_df <- data.frame(newloca$RK.GLM)

coordinates_df <- data.frame(coordinates(newloca))

pred_beta=cbind(coordinates_df,attributes_df)
names(pred_beta)[3]="beta"


p2=ggplot(pred_beta) +
  scale_color_gradient(expression(Richness),low = "blue", high = "yellow")+
  geom_point(data=pred_beta,pch=15,aes(x=lon, y=lat,color=beta), size=0.275)+
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
  xlab("Predicted plot-level beta diversity")+
  ylab("")+
  xlim(-175,-42)+
  geom_point(data=plot_diversity,pch=21,aes(x=lon,y=lat,fill=project),alpha=1)+
  scale_fill_manual("Project",breaks=c("dob","neon"),label=c("DoB","NEON"),values=c("red","seagreen"))

###


k=summary(GLM_beta)

effect_beta=k$coefficients[2:11,]%>%data.frame()

p21=ggplot()+
  
  geom_segment(data=effect_beta,aes(x=Estimate-1.96*Std..Error ,y=1:10,yend=1:10,xend=Estimate+1.96*Std..Error),color="gray")+
  geom_point(data=effect_beta,aes(x=Estimate,y=1:10),color="blue",size=4)+
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
  ggtitle("Plot-level beta diversity (18.4%)") 

p1=ggplotGrob(p1)
p2=ggplotGrob(p2)
p3=ggplotGrob(p3)
p11=ggplotGrob(p11)
p21=ggplotGrob(p21)
p31=ggplotGrob(p31)
p1$widths=p11$widths
p2$widths=p21$widths
p3$widths=p31$widths

plot_grid(p11,p21,p31,p1,p2,p3,col=3)
