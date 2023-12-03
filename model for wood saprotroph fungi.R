# model for woodsap
a=list()
for (i in 1:dim(woosap_c_ranall)[1])
{
  d=data.frame(t(woosap_c_ranall[i,3:32]))
  names(d)="id"
  a[[i]]=d
}

b=a[[1]]
for(i in 2:dim(woosap_c_ranall)[1])
{
  b=rbind(b,a[[i]])
}

plotID=rep(woosap_c_ranall$a1,each=30)

woosap_c_ranall_30=cbind(plotID,b)# for each plot,with 30 estimated z values

a <- list()
for (i in 1:dim(woosap_z_ranall)[1])
{
  d <- data.frame(t(woosap_z_ranall[i, 3:32]))
  names(d) <- "id"
  a[[i]] <- d
}

b <- a[[1]]
for (i in 2:dim(woosap_z_ranall)[1])
{
  b <- rbind(b, a[[i]])
}

woosap_z_ranall_30 <- b




woosap_model=cbind(woosap_c_ranall_30,b)
names(woosap_model)=c("plotID","logc","z")

woosap_model=merge(woosap_model,model_var,by="plotID")
woosap_model=subset(woosap_model,siteIDD!="GUAN"&z<10&fine>0&rootc>0&richness>0)# only 104 plots from 33 sites

woosap_model_rich=subset(woosap_model,siteIDD!="GUAN"&z<10&richness>0)# only 104 plots from 33 sites

woosap_model_rich <- cbind(woosap_model_rich, c = 2.71828^woosap_model_rich$logc)

#standardized data
woosap_model[,c(5:28)]=apply(woosap_model[,c(5:28)],2,range01)
woosap_model_rich[,c(5:28)]=apply(woosap_model_rich[,c(5:28)],2,range01)
# colinearity
ggcorrplot(cor(woosap_model[,c(2,3,5:27)]), hc.order = TRUE, type = "lower", lab = TRUE)#
#bold was related with soil OC and was excluded,cec was related to soil OC and was removed
#bio4 was related with bio1 and was removed
#coarse (removed) was related with fine and was removed
#cec(removed) and soil OC
#root cn and rootn(removed)
#bold (removed)and soil OC
#bio1 and bio4(removed)
3
# build a model for the woosap guild,with 104 plots included, rich and d13

mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1 +fine+ d15N +d13C +rootc +rootcn +(1 |siteIDD/plotID),data=woosap_model)
step(mod)

mod_woosap=lmer(z ~ c + bio15 + funrich + bio1 + rootcn + (1 | siteIDD/plotID),woosap_model)

p6=plot_model(mod_woosap,axis.labels = c(expression("Root"["cn"]),"MAT", "Fun.rich","Pre.seas."),colors=c("blue","red"),rm.terms = "c",title="Woo.sap. (N=103)")

plot_model(mod_woosap)
## excluding the root traits

mod=lmer(z ~ c+organicCPercent + ph+ nitrogen+ sand +bio2 +bio8+ bio18 +bio12+ bio15  +spei+ richness+funrich +bio1  +(1 |siteIDD/plotID),data=woosap_model_rich)
step(mod)

mod_woosap_rich=lmer(z ~ c + sand + bio2 + bio15 + spei + richness + funrich + bio1 + (1 | siteIDD/plotID),woosap_model_rich)
plot_model(mod_woosap_rich)

p60=plot_model(mod_woosap_rich,axis.labels = c("MAT", "Fun.rich","Pla.rich","Spei","MAP", "MDR","Sand"),colors=c("blue","red"),rm.terms = "c",title="Woo.sap. (N=427)")

