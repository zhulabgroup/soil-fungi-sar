# model for epiphy
a=list()
for (i in 1:dim(epiphy_c_ranall)[1])
{
  d=data.frame(t(epiphy_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(epiphy_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(epiphy_c_ranall$a1,each=30)

epiphy_c_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values

epiphy_model=cbind(epiphy_c_ranall_30,epiphy_z_ranall_30["z"])
names(epiphy_model)[2]="logc"
epiphy_model=merge(epiphy_model,model_var,by="plotID")
epiphy_model=subset(epiphy_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&rich>0)# only 104 plots from 33 sites



#standardized data
epiphy_model[,c(2,5:27)]=apply(epiphy_model[,c(2,5:27)],2,range01)
# colinearity
ggcorrplot(cor(epiphy_model[,c(2,3,5:27)]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil OC and was excluded,cec was related to soil OC and was removed

#coarse (removed) was related with fine and was removed
#cec(removed) and soil OC
#root cn and rootn(removed)
#bold (removed)and soil OC
#bio1 and bio4(removed)
3
# build a model for the epiphy guild,with 104 plots included, fungal rich was most fluencial

mod=lmer(z ~ logc+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ rich+funrich +bio1 +fine+ d15N +d13C +rootc +rootcn +(1 |siteIDD/plotID),data=epiphy_model)
##
