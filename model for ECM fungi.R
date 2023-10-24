# model for ECM
a=list()
for (i in 1:dim(ecm_c_ranall)[1])
{
  d=data.frame(t(ecm_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(ecm_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(ecm_c_ranall$a1,each=30)

ecm_c_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values

ecm_model=cbind(ecm_c_ranall_30,ecm_z_ranall_30["z"])
names(ecm_model)[2]="logc"
ecm_model=merge(ecm_model,model_var,by="plotID")
ecm_model=subset(ecm_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&rich>0)# only 104 plots from 33 sites
# head(ecm_model)
# consider the climate and soil data
ecm_model1=ecm_model[,c(1:27)]


#standardized data
ecm_model1[,c(2,5:27)]=apply(ecm_model1[,c(2,5:27)],2,range01)
# colinearity
ggcorrplot(cor(ecm_model1[,c(2,3,5:27)]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil OC and hence was excluded,cec was related to soil C and was removed
#bio4 was related with bio1 and was removed
#coarse (removed) and fine were related 
#cec(removed) and soil OC
#root cn and rootn(removed)
#bold (removed)and soil OC
#bio1 and bio4(removed)
3
# build a model for the ecm guild,with 104 plots included, soil pH and fundiversity

mod=lmer(z ~ logc+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ funrich +bio1 +rich+fine+ d15N +d13C +rootn +   rootc +rootcn +(1 |siteIDD/plotID),data=ecm_model1)
##
