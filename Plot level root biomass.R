# root-biomass
# data describtion
#https://data.neonscience.org/data-products/DP1.10067.001

root.data <- loadByProduct(dpID="DP1.10067.001")#
root.data[[1]]
root.data[[2]]
root.data[[3]]# about sampling location and time
root.data[[4]]# root chemical traits
root.data[[5]]#dry mass with different root orders
root.data[[6]]# site description?
root.data[[7]]#description
root.data[[8]]# site description?
root.data[[9]]# site description?
root.data[[10]]# site description
head(root.data[[4]])
unique(root.data[[4]]["sampleType"])

root.mass=data.frame(root.data[[5]])# the biomass data
root.mass=root.mass[,c("siteID","plotID","sampleID","subsampleID","sizeCategory","rootStatus","dryMass","collectDate")]
# lump samples within a single core and plot
unique(root.mass$rootStatus)
# before 2016 roots were sorted based four orders and two root status but after 2016 may, roots were sorted into three orders
# so we need to use a criteria to lump the data
code=substr(root.mass$sampleID,12,14)
root.mass=cbind(root.mass,code)# samples of the same code are from the same sampling strip
face=substr(root.mass$sampleID,25,29)
root.mass=cbind(root.mass,face)# the face of the sample, samples of the same code and face are from the same core
head(root.mass)
subset(root.mass,plotID=="BART_074")
# to get the unique id of a soil core
core=paste(root.mass$plotID,"_",root.mass$code,"_",root.mass$face)
root.mass=cbind(root.mass,core)
coreID=unique(root.mass$core)
# it seems that some samples were remeasured so i took the remeasured samples
core_mass=list()
for (i in 1:length(coreID)){
  a=subset(root.mass,core==coreID[i])
  # get the mean of all the roots
  core_mass[[i]]=aggregate(dryMass~sizeCategory,data=a,FUN=mean) 
}

d=numeric()# the number of root grouop within each core(e.g. 0-1,2-10)
for (i in 1:length(coreID)){
  d[i]=dim(core_mass[[i]])[1]
}

core_mass1=core_mass[[1]]
for (i in 2:length(coreID)){
  core_mass1=rbind(core_mass1,core_mass[[i]])
}
## adding the coreID to the rbined data
coreID1=rep(coreID,times=d)

core_mass1=cbind(core_mass1,coreID1)# data.frame of root mass of each class for each soil core
# pool fine(<2mm) and coarse roots(2-10 mm) within each core

core_mass1[core_mass1==c("0-05","05-1","1-2","0-1")]="fine"
core_mass1[core_mass1=="2-10"]="coarse"

# aggregate samples by fine and coarse roots

core_mass_type=matrix(ncol=2,nrow=length(coreID))
for (i in 1:length(coreID)){
  a=subset(core_mass1,coreID1==coreID[i])
  # get the mean of all the roots
  core_mass_type[i,1]=aggregate(dryMass~sizeCategory,data=a,FUN=mean) [1,2]
  core_mass_type[i,2]=aggregate(dryMass~sizeCategory,data=a,FUN=mean) [2,2]
}
core_mass_type=cbind(core_mass_type,coreID)
core_mass_type=data.frame(core_mass_type)
names(core_mass_type)=c("coarse","fine","coreID")
core_mass_type$coarse=as.numeric(core_mass_type$coarse)
core_mass_type$fine=as.numeric(core_mass_type$fine)
head(core_mass_type)
# to look at the number of plot per site
plotID=substr(core_mass_type$coreID,1,8)
core_mass_type=cbind(core_mass_type,plotID)
# in some case, only one core was sampled, either in the S or in the N direction