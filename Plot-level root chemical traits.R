# get the plot-level root traits#
library(neonUtilities)
root.data <- loadByProduct(dpID="DP1.10067.001")#
root.chemi=root.data[[4]]
head(root.chemi)
root.chemi=root.chemi[,c("siteID","plotID","d15N","d13C","nitrogenPercent","carbonPercent","CNratio")]
head(root.chemi)
table(root.chemi$plotID)
root.chemi.mean=aggregate(root.chemi[,3:7], by=list(root.chemi$plotID), mean)
names(root.chemi.mean)[1]="plotID"
write.csv(root.chemi.mean,"root.chemi.mean.csv")