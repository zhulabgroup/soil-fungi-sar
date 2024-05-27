# get the model
library(doParallel)
library(glmnet)
neon_dob_prevalent <- readRDS("~/soil-sar/SDM/neon_dob_prevalent_v4.1.Rds")

set.seed(10101)
train_set <- sample(c(TRUE, FALSE), size=length(sample_sums(neon_dob_prevalent)),
                    prob=c(0.70, 0.30), replace=TRUE)
test_set <- !train_set
neon_dob_prevalent_train <- prune_samples(train_set, neon_dob_prevalent)
neon_dob_prevalent_test <- prune_samples(test_set, neon_dob_prevalent)
n_taxa <- 8597
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent_train)), size=n_taxa))
neon_dob_prevalent_train_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_train)[spp_subset_ind],
                                              neon_dob_prevalent_train)
neon_dob_prevalent_test_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_test)[spp_subset_ind],
                                             neon_dob_prevalent_test)

# define the variables

# climate_only
terms1 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
            "mat_celsius_2", "map_mm_2")
# climate_and_soil
terms2 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "nitrogenPercent", "organicCPercent")
# important_vars (a posteriori naming)
terms3 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")# select the variables
# climate_and_ph
terms4 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH")

##

predictors_train1 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms1], "data.frame")
predictors_train_std1 <- scale(predictors_train1, center=TRUE, scale=TRUE)
complete_records1 <- which(complete.cases(predictors_train_std1))

predictors_train2 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms2], "data.frame")
predictors_train_std2 <- scale(predictors_train2, center=TRUE, scale=TRUE)
complete_records2 <- which(complete.cases(predictors_train_std2))

predictors_train3 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms3], "data.frame")
predictors_train_std3 <- scale(predictors_train3, center=TRUE, scale=TRUE)
complete_records3 <- which(complete.cases(predictors_train_std3))

predictors_train4 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms4], "data.frame")
predictors_train_std4 <- scale(predictors_train4, center=TRUE, scale=TRUE)
complete_records4 <- which(complete.cases(predictors_train_std4))

###
# do not remember to load the package otherwise NA would be returned

registerDoParallel(3)
lognets_train1 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std1)[complete_records1,],
                     y01[complete_records1], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

lognets_train2 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)# just a way for cross validation
  tryCatch(cv.glmnet(as.matrix(predictors_train_std2)[complete_records2,],
                     y01[complete_records2], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))# if an error occur re
}

lognets_train3 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std3)[complete_records3,],
                     y01[complete_records3], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

lognets_train4 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std4)[complete_records4,],
                     y01[complete_records4], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}
stopImplicitCluster()

lognet_thresholds_train1 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train1,
  terms = terms1)

lognet_thresholds_train2 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train2,
  terms = terms2)

lognet_thresholds_train3 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train3,
  terms = terms3)

lognet_thresholds_train4 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train4,
  terms = terms4)

###

newdata1 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms1]
newdata_std1 <- scaleToReference(newdata1, predictors_train1)

newdata2 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms2]
newdata_std2 <- scaleToReference(newdata2, predictors_train2)

newdata3 <- as(newdata, "data.frame")[,terms3]

newdata_std3 <- scaleToReference(newdata3, predictors_train3)

newdata_std3=newdata_std3[complete.cases(newdata_std3),]# get the completed cases

coord=newdata[which(complete.cases(newdata)),][,c(1,2)]

newdata4 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms4]
newdata_std4 <- scaleToReference(newdata4, predictors_train4)



registerDoParallel(3)
lognet_preds1 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train1[[i]], NA)) return(NA)
  prediction <- predict(lognets_train1[[i]],
                        newx = as.matrix(newdata_std1),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train1[[i]]
  as.numeric(classification)
}

lognet_preds2 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train2[[i]], NA)) return(NA)
  prediction <- predict(lognets_train2[[i]],
                        newx = as.matrix(newdata_std2),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train2[[i]]
  as.numeric(classification)
}

lognet_preds3 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train3[[i]], NA)) return(NA)
  prediction <- predict(lognets_train3[[i]],
                        newx = as.matrix(newdata_std3),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train3[[i]]
  as.numeric(classification)
}

lognet_preds4 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train4[[i]], NA)) return(NA)
  prediction <- predict(lognets_train4[[i]],
                        newx = as.matrix(newdata_std4),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train4[[i]]
  as.numeric(classification)
}


stopImplicitCluster()

save(lognet_preds3 ,file="lognet_preds3.RData")

# # get the species presence dat across all the sites

presence=matrix(nrow =n_taxa,ncol=dim(newdata_std3)[1])

for(i in 1:n_taxa){
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  presence[i,]=lognet_preds3[[i]]
}


presence=cbind(coord,t(presence))

richness=rowSums(presence[,3:8599],na.rm = TRUE)

presence_rich=cbind(coord,richness)

ggplot(presence_rich) +
  geom_tile(data=presence_rich,aes(x=lon, y=lat,color=richness),size=0.275)+
  scale_color_gradient(expression("Richness"),low = "blue", high = "yellow")+
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
  xlab("Predicted species richness")+
  ylab("")+
  xlim(-175,-42)
# get the coordinates for the new locations

save(presence_rich,file="presence_rich.RData")

###









validation_predictions <- list(
  lognet_preds1,
  lognet_preds2,
  lognet_preds3,
  lognet_preds4
)
names(validation_predictions) <- c("climate_only", "climate_and_soil", "important_vars", "climate_and_ph")

library(caret)
accuracy <- matrix(nrow=n_taxa, ncol=length(validation_predictions), dimnames=list(spp_subset_ind, names(validation_predictions)))
for(i in 1:n_taxa) {
  sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  y <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  for(mod in colnames(accuracy)) {
    accuracy[i,mod] <- mean(y == validation_predictions[[mod]][[i]], na.rm=TRUE)
  }
}


sensitivity <- specificity <- matrix(nrow=n_taxa, ncol=length(validation_predictions), dimnames=list(spp_subset_ind, names(validation_predictions)))
for(i in 1:n_taxa) {
  sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  y <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  if(all(c(0,1) %in% y)) {
    confusions <- lapply(
      validation_predictions,
      function(x) tryCatch(
        confusionMatrix(factor(x[[i]], levels=c(0,1)),
                        factor(y, levels=c(0,1)), positive="1"),
        error = function(e) return(NA))
    )
    sensitivity[i,] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Sensitivity"], error = function(e) return(NA)))
    specificity[i,] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Specificity"], error = function(e) return(NA)))
  } else {
    sensitivity[i,] <- NA
    specificity[i,] <- NA
  }
}


apply(accuracy, 2, mean, na.rm=TRUE)
apply(sensitivity, 2, mean, na.rm=TRUE)
apply(specificity, 2, mean, na.rm=TRUE)

apply(accuracy, 2, sd, na.rm=TRUE)
apply(sensitivity, 2, sd, na.rm=TRUE)
apply(specificity, 2, sd, na.rm=TRUE)

apply(sensitivity + specificity - 1, 2, mean, na.rm=TRUE) # True skill statistic
apply(sensitivity + specificity - 1, 2, sd, na.rm=TRUE) # True skill statistic

#

# based on the new locations to make predictions at the 10-min resolution

head(data_mean)

library(tidyr)
result <- separate(temp, coordsSite, into = c("lon", "lat"), sep = "_")# the coordinates used for the prediction of species 

# to get the soil variables of these locations


