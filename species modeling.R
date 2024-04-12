# devide the map into small grids
library(sf)
north_america <- ne_countries(continent="North America")

st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)
b <- as(extent(-170,-55,18,72), "SpatialPolygons")
crs(b) <- crs(north_america)
north_america_cropped  <- st_crop(north_america, b)
plot(north_america_cropped)


#
sf_map <- st_transform(north_america_cropped, crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

grid_10min <- st_make_grid(north_america_cropped, cellsize = 2)

grid_coords <- st_coordinates(grid_10min)%>%data.frame()->dd

# based on the variables

# a new way to define the study region
library(sp)
library(raster)
library(dplyr)
library(rnaturalearth)
library(sf)

us_extent <- extent(-170, -55, 18, 72)
raster_layer <- raster(ext = us_extent, res = 1/6)
grid_polygons <- as(raster_layer, "SpatialPolygonsDataFrame")
grid_coordinates <- coordinates(grid_polygons)%>%data.frame()
# extract the variables based on the coordinates


r_present <- getData("worldclim",var="bio",res=10)
r_present <- r_present[[c(1,4,12,15)]]
names(r_present) <- c("mat_celsius","temp_seasonality","map_mm","prec_seasonality")

# Run necessary transformations on wordclim-provided temperature data
r_present$mat_celsius <- r_present$mat_celsius/10
r_present$temp_seasonality <- r_present$temp_seasonality/1000

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland
north_america <- ne_countries(continent="North America")
st_is_valid(north_america, reason = TRUE)[!st_is_valid(north_america)]
sf_use_s2(use_s2 = FALSE)
plot(north_america)

b <- as(extent(-170,-55,18,72), "SpatialPolygons")# based on this to get the variables

crs(b) <- crs(north_america)
north_america_cropped  <- st_crop(north_america, b)






# construct the best models with 100 species for an example
# determine the important variables


set.seed(10101)
train_set <- sample(c(TRUE, FALSE), size=length(sample_sums(neon_dob_prevalent)),
                    prob=c(0.70, 0.30), replace=TRUE)
test_set <- !train_set

neon_dob_prevalent_train <- prune_samples(train_set, neon_dob_prevalent)
neon_dob_prevalent_test <- prune_samples(test_set, neon_dob_prevalent)
n_taxa <- 50
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent_train)), size=n_taxa))
neon_dob_prevalent_train_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_train)[spp_subset_ind],
                                              neon_dob_prevalent_train)
neon_dob_prevalent_test_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_test)[spp_subset_ind],
                                             neon_dob_prevalent_test)
##

terms1 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
            "mat_celsius_2", "map_mm_2")
# climate_and_soil
terms2 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "nitrogenPercent", "organicCPercent")
# important_vars (a posteriori naming)
terms3 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")
# climate_and_ph
terms4 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH")

predictors_train1 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms1], "data.frame")#
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


registerDoParallel(5)
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
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std2)[complete_records2,],
                     y01[complete_records2], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
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




# 

# some codes in the findthreshold funciton and take the first model as an example


terms=terms1

covars <- covarNamesFromLognets(lognets_train1)

pred_data <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms], "data.frame")
pred_data_std <- scale(pred_data, center=TRUE, scale=TRUE)
ind_cols <- match(covars, colnames(pred_data_std))
pred_data_std <- pred_data_std[,ind_cols]# to see which predictors are match

# to see how the thres was determined, let's take the first species as an examble 
thres <- foreach(i = seq_along(taxa_names(neon_dob_prevalent_train_sppsub)), .combine = c) %dopar% {
  if(identical(models[[i]], NA)) return(NA)
  sp <- taxa_names(neon_dob_prevalent_train_sppsub)[i]
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,sp]) > 0
  lognet_pred <- predict(lognets_train1[[i]], newx = pred_data_std, type="response", s = "lambda.min")# the response for each
  t_grid <- seq(min(lognet_pred, na.rm=TRUE), max(lognet_pred, na.rm=TRUE),
                by=(max(lognet_pred, na.rm=TRUE)-min(lognet_pred, na.rm=TRUE))/100)# there 100 values, something like sequentially sampled
  
  eval <- dismo::evaluate(p = lognet_pred[which(y01)], a = lognet_pred[which(!y01)], tr = t_grid)
  # plot(eval@TPR - eval@TNR)
  if(identical(opt, "eq")) {
    # Threshold to minimize diff b/w sensitivity and specificity
    t <- eval@t[which.min(abs(eval@TPR - eval@TNR))]
    return(t)
  }
  if(identical(opt, "maxtss")) {
    # threshold to maximize true skill statistic
    t <- eval@t[which.max(eval@TPR + eval@TNR)]
    return(t)
  }
}
stopImplicitCluster()
message("findThresholds: Ended at ", Sys.time())
return(thres)
}

#have not go through but it is related to the obseved presence and the predicted probality, based on which to get the threshold

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

newdata1 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms1]# variables in the test data
newdata_std1 <- scaleToReference(newdata1, predictors_train1)# variables in the train data/ f

newdata2 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms2]
newdata_std2 <- scaleToReference(newdata2, predictors_train2)

newdata3 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms3]
newdata_std3 <- scaleToReference(newdata3, predictors_train3)

newdata4 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms4]
newdata_std4 <- scaleToReference(newdata4, predictors_train4)# why this step
###3

registerDoParallel(5)
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

# validation

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

##

coef_matrix_lognet <- matrix(NA, nrow = length(lognets_train2), ncol=10)
rownames(coef_matrix_lognet) <- names(lognets_train2)
colnames(coef_matrix_lognet) <- rownames(coef(lognets_train2[[2]], s="lambda.min"))
system.time(
  for(i in 1:nrow(coef_matrix_lognet)) { # takes 1 minute
    if(!identical(lognets_train2[[i]], NA)) {
      coef_matrix_lognet[i,] <- coef(lognets_train2[[i]], s="lambda.min")[,1]
    }
  }
)
head(coef_matrix_lognet)

colnames(coef_matrix_lognet) <- c("Intercept", "MAT", "TSEA", "MAP", "PSEA",
                                  "MAT2", "MAP2", "pH", "%N", "%C")

coef_long_lognet <- coef_matrix_lognet[,-1] %>%
  as.data.frame() %>%
  tidyr::pivot_longer(MAT:`%C`)


mean_abs_coef <- coef_long_lognet %>%
  group_by(name) %>%
  summarise(mean = mean(abs(value), na.rm=TRUE)) %>%
  arrange(desc(mean)) %>%
  mutate(name = factor(name, levels=unique(name)))
mean_abs_coef

median_abs_coef <- coef_long_lognet %>%
  group_by(name) %>%
  summarise(median = median(abs(value), na.rm=TRUE)) %>%
  arrange(desc(median)) %>%
  mutate(name = factor(name, levels=unique(name)))
median_abs_coef


coef_long_lognet %>%
  mutate(name = factor(name, levels=unique(median_abs_coef$name)))


###

set.seed(10101)
n_taxa <- 500
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent)), size=n_taxa))
neon_dob_prevalent_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent)[spp_subset_ind],
                                        neon_dob_prevalent)


registerDoParallel(4)
lognets_train2 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  cat("\r", paste(paste0(rep("*", round(i / 1, 0)), collapse = ""), i, collapse = "")) # informs the processing
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std2)[complete_records2,],
                     y01[complete_records2], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}
stopImplicitCluster()





set.seed(10101)
n_taxa <- 100
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent)), size=n_taxa))
neon_dob_prevalent_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent)[spp_subset_ind],
                                        neon_dob_prevalent)

# Define variables
terms <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
           "mat_celsius_2", "map_mm_2", "soilInCaClpH", "nitrogenPercent", "organicCPercent")

neon_dob_permuted <- permuteClimate(neon_dob_prevalent_sppsub, randomseed=1010103)
predictors <- as(sample_data(neon_dob_permuted)[,terms], "data.frame")
predictors_std <- scale(predictors, center=TRUE, scale=TRUE)# models were built with none variables
complete_records <- which(complete.cases(predictors_std))

lognets_permuted <-  foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_permuted)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_std)[complete_records,],
                     y01[complete_records], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

###

thresholds_permuted <- findThresholds(neon_dob_permuted, lognets_permuted, terms = terms)# for the species to be presenting

metrics_2d_permuted <- getNicheMetrics2D(neon_dob_permuted, lognets_permuted, thresholds_permuted, 
                                         "mat_celsius", "map_mm",r_present_northam,n_bins_per_axis = 10, ncores=4,  
                                         pred_vars = terms,constant_values = NULL,x1_range=c(-12.2, 25.1),x2_range=c(116, 2556))


p_permuted <- plotNicheMetrics2D(neon_dob_permuted, metrics_2d_permuted, pred_vars = terms)

# to get the distribution probability of the 100 species, and this is just based on the 121 sites




getNicheMetrics2D2 <- function(physeq, models, thresholds, x1, x2,
                               x1_range=NULL, x2_range=NULL, constant_values=NULL,
                               r_climate, n_bins_per_axis=50,
                               hull=c("2D", "4D", "none"), inflate=0, ncores=1,
                               progressbar=TRUE, dev_version=FALSE,
                               pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")) {
  
  if(!x1 %in% pred_vars | !x2 %in% pred_vars) {
    stop("`x1` and `x2` must be among `pred_vars`")
  }
  hull <- match.arg(hull, c("2D", "4D", "none"))
  n_gradient_steps <- n_bins_per_axis + 1
  n_taxa <- length(taxa_names(physeq))
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]# the same as the intial one
  obs_env <- obs_env[complete.cases(obs_env),]
  
  env_medians <- apply(obs_env, 2, function(x) median(x, na.rm=TRUE))
  if(is.null(constant_values)) {
    env_constant <- env_medians
  } else {
    if(length(constant_values) != length(pred_vars)) {
      stop("`constant_values` must have same length as `pred_vars`.
           To use the median of observed values for a variable, enter NA for
           that variable.")
    }
    env_constant <- constant_values
    env_constant[which(is.na(constant_values))] <- env_medians[which(is.na(constant_values))]
    names(env_constant) <- pred_vars
  }
  message("Using as constant unless otherwise specified: ", paste0(pred_vars, " ", env_constant, ", "))
  
  if(is.null(x1_range)) {
    x1_range <- range(obs_env[,x1], na.rm=TRUE)
    message("Using as `x1_range`: ", paste0(x1_range, collapse=", "))
  }
  if(is.null(x2_range)) {
    x2_range <- range(obs_env[,x2], na.rm=TRUE)
    message("Using as `x2_range`: ", paste0(x2_range, collapse=", "))
  }
  
  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- tidyr::expand_grid(
    # quantize(obs_env[,x1], n_gradient_steps),
    # quantize(obs_env[,x2], n_gradient_steps)
    quantize(x1_range, n_gradient_steps),
    quantize(x2_range, n_gradient_steps)
  )
  names(surface_var) <- c(x1, x2)
  for(p in terms) {
    if(!p %in% colnames(surface_var)) {
      surface_var[,p] <- env_constant[p]
    }
  }
  
  covars <- covarNamesFromLognets(models)
  surface_var <- includeQuadraticTerms(surface_var)# all ^2 terms
  surface_var <- surface_var[,match(covars,colnames(surface_var))]# only select the variables that we interested
  predictors <- as(sample_data(neon_dob_permuted)[,covars], "data.frame")# why this is not keep constant
  
  predictors <- predictors[,match(covars, colnames(predictors))]
  surface_var_scaled <- scaleToReference(surface_var, predictors)# the former being the one with constants but the later is the observed
  
  # predictors=the data in the initial data 
  
  #surface_var= variables included in the initial model, with some kept constant
  # it suggests that 79 sites were inside the hull of the 
  
  #obs_env[,x1]=observed variables
  
  
  # If constraining by the convex hull around observations for the two climate axes of interest(species occurence data? be calculated?)
  if(identical(hull, "2D")) {
    x <- obs_env[,x1]
    y <- obs_env[,x2]
    ch_ind <- chull(x, y)
    infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))# observed values
    infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
    inhull <- pracma::inpolygon(surface_var[,x1,drop=TRUE], surface_var[,x2,drop=TRUE], infl_x, infl_y, boundary=TRUE)
  }# to check if the are within the polygon
  
  # Otherwise, inhull is NA
  if(identical(hull, "none")) {
    inhull <- NA
  }
  
  surface_var_scaled <- surface_var_scaled[inhull,]
  
  surface_suitability <- matrix(NA, nrow = nrow(surface_var), ncol = length(taxa_names(physeq)))
  
  registerDoParallel(ncores)
  temp <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  
  surface_suitability[inhull,] <- temp
  
  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))
  
  colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)
  
  
  
  
  
  
