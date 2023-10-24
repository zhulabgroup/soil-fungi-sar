# for specific functional guilds
# for the ACM functional guild
# get the explainatory variables from the full data based on the plotID
head(model_data)
model_var=model_data[,c(2,3,6:28)]
model_var=unique(model_var)

a=list()
for (i in 1:dim(acm_c_ranall)[1])
{
  d=data.frame(t(acm_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(acm_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(acm_c_ranall$a1,each=30)

acm_c_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values

acm_model=cbind(acm_c_ranall_30,acm_z_ranall_30["z"])
names(acm_model)[2]="logc"
acm_model=merge(acm_model,model_var,by="plotID")
acm_model=subset(acm_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&rich>0)# only 87 plots from 33 sites
# head(acm_model)
# consider the climate and soil data
acm_model1=acm_model[,c(1:17,19,27)]
acm_model1=cbind(acm_model1,c=2.71828^acm_model1$logc)
acm_model1=cbind(acm_model1,ct=log(acm_model1$c))

#standardized data
acm_model1[,5:20]=apply(acm_model1[,5:20],2,range01)
# colinearity
ggcorrplot(cor(acm_model1[,c(3,5:18)]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil c and hence was excluded
# build a model for the acm guild,355 plots

mod=lmer(z ~ organicCPercent + ph+ nitrogen+ cec+ sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +spei+ funrich +bio1 +(1 |siteIDD/plotID),data=acm_model1)
#

# consider the root model for the acm
a=list()
for (i in 1:dim(acm_c_ranall)[1])
{
  d=data.frame(t(acm_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(acm_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(acm_c_ranall$a1,each=30)

acm_c_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values

acm_model=cbind(acm_c_ranall_30,acm_z_ranall_30["z"])
names(acm_model)[2]="logc"
acm_model=merge(acm_model,model_var,by="plotID")
acm_model2=subset(acm_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&rich>0)# only 87 plots from 33 sites
# head(acm_model)
# consider the climate and soil data


#standardized data
acm_model2[,c(2,5:27)]=apply(acm_model2[,c(2,5:27)],2,range01)
# colinearity
ggcorrplot(cor(acm_model2[,c(2,5:27)]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil c and hence was excluded
#root n (excluded) and root cn were correlated
# build a model for the acm guild,355 plots

mod=lmer(z ~ logc+organicCPercent + ph+ nitrogen+ cec+ sand +bio2 +bio8+ bio18+  bio4 +bio12+ bio15  +spei+ rich+funrich +bio1 +fine +d15N + d13C+ rootc+ rootcn+ bio1+(1 |siteIDD/plotID),data=acm_model2)
#
# for the ecm guild

