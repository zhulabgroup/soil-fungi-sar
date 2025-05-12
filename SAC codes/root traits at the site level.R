# look the plot soil variables based on the neon sites

library(neonUtilities)

soil.data <- loadByProduct(dpID = "DP1.00096.001") #

soil.data <- loadByProduct(dpID = "DP1.10086.001") #what does the nitrogen and organic carbon percent mean?
###

#Standing stock of fine and coarse root biomass plus root carbon (C) and nitrogen (N) concentrations and stable isotopes from soil Megapits
# get the root trait data based on the megapit data
library(neonUtilities)
root.mega<- loadByProduct(dpID = "DP1.10066.001") 
names(root.mega)
root.mega[[1]]
root.mega[[2]]
root.mega[[3]]
root.mega[[4]]
root.mega[[5]]# root n, cn and isotopes, basd on the site
root.mega[[6]]#biomass, 
root.mega[[7]]# dry mass\
root.mega[[8]]
root.mega[[9]]
root.mega[[10]]# data ends here

