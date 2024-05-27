# functions used in the excertise

findThresholds <- function(physeq, models, terms = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality"),
                           s=c("min", "1se"), opt=c("eq","maxtss"), ncores=1) {
  message("findThresholds: Started at ", Sys.time())
  s <- match.arg(s, c("min", "1se"))
  opt <- match.arg(opt, c("eq","maxtss"))
  s <- switch(s, "min" = "lambda.min", "1se" = "lambda.1se")
  
  covars <- covarNamesFromLognets(models)
  pred_data <- as(sample_data(physeq)[,terms], "data.frame")
  pred_data_std <- scale(pred_data, center=TRUE, scale=TRUE)
  ind_cols <- match(covars, colnames(pred_data_std))
  pred_data_std <- pred_data_std[,ind_cols]
  
  registerDoParallel(ncores)
  thres <- foreach(i = seq_along(taxa_names(physeq)), .combine = c) %dopar% {
    if(identical(models[[i]], NA)) return(NA)
    sp <- taxa_names(physeq)[i]
    y01 <- as.vector(otu_table(physeq)[,sp]) > 0
    lognet_pred <- predict(models[[i]], newx = pred_data_std, type="response", s = s)
    t_grid <- seq(min(lognet_pred, na.rm=TRUE), max(lognet_pred, na.rm=TRUE),
                  by=(max(lognet_pred, na.rm=TRUE)-min(lognet_pred, na.rm=TRUE))/100)
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

