
cv.glmmLasso <- function(fix, rnd, data, 
                         family = stats::gaussian(link = "identity"), 
                         kfold = 5, lambdas = NULL, nlambdas = 100, 
                         lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001), 
                         loss,
                         lambda.final=c('lambda.1se', 'lambda.min'),
                         ...)
{
  lambda.final <- match.arg(lambda.final)
  
  if(missing(loss))
  {
    # switch allows us to do take the family arg as assign the appropriate 
    # loss function 
    loss <- switch(family$family, 
                   'gaussian' = calc_mse,
                   'binomial' = calc_logloss,
                   'multinomial' = calc_multilogloss,
                   'poisson' = calc_deviance)
  }
  
  x <- useful::build.x(fix, data)
  nobs <- nrow(x)
  nvars <- ncol(x)
  
  # if lambda isn't specified by user, build the lambdas vector, this is 
  # static for all k folds
  if (is.null(lambdas))
  {
    # building the lambda vector
    lambdas <- buildLambdas(fix = fix,
                            rnd = rnd,
                            data = data, 
                            nlambdas = nlambdas, 
                            lambda.min.ratio= lambda.min.ratio)   
  }
  
  
  
  
  # building data frame to map a specific row to kth group
  # column 1 is the row, column 2 is a randomly assigned group
  # number of groups is determined by kfold value  
  rowDF <- tibble::tibble(
    row = seq(nobs),
    group = sample(rep(seq(kfold), length.out=nobs), replace = FALSE)
  )
  
  # sorting by group 
  rowDF <-  dplyr::arrange(rowDF, .data$group)
  
  
  #instantiating list to hold loss and models for each fold
  lossVecList <- vector(mode = 'list', length = kfold)
  modList_foldk <- vector(mode = 'list', length = kfold)
  
  for(k in 1:kfold)
  {
    testIndices <- dplyr::filter(rowDF, .data$group == k) %>% dplyr::pull(row)
    trainIndices <- rowDF$row[-testIndices]
    
    # fitting model
    # modList_foldk is a glmmLasso_MultLambdas object, which is a list of 
    # glmmLasso objects
    
    # for showing lambda at each iterations
    # message(sprintf('Round: %s\n ', k))
    modList_foldk[[k]] <- glmmLasso_MultLambdas(fix = fix,
                                                rnd = rnd,
                                                data = data %>% dplyr::slice(trainIndices),
                                                family = family,
                                                lambdas = lambdas,
                                                nlambdas = nlambdas,
                                                lambda.min.ratio = lambda.min.ratio,
                                                ...)
    
    
    
    # hacky way of getting the response variable out of the         
    response_var <- fix[[2]] %>% as.character()
    
    # pulling out actual data
    actualDataVector <- data %>% dplyr::slice(testIndices) %>% 
      dplyr::pull(response_var)
    
    # predicting values for each of the glmmLasso model (100 lambda) 
    # using matrix form for easier error calculation in loss()
    
    predictionMatrix <- predict.glmmLasso_MultLambdas(
      object = modList_foldk[[k]],
      newdata = data %>% dplyr::slice(testIndices)
    )
    
    # employing the loss function in form loss(actual,predicted)
    # using loss function, calculating a list of loss values for each vector 
    # of prediction
    # which comes from a glmmLasso model with a specific lambda 
    # storing loss values for each fold
    
    # TODO: think an error is thrown here 
    lossVecList[[k]] <- loss(actual = actualDataVector, predicted = predictionMatrix)
    # each element of this list should be 1 x nlambdas
  }
  
  #building matrix (k by nlambdas) to help calculate cross-validated mean error
  cvLossMatrix <- do.call(what = rbind, args = lossVecList)
  
  cvm = colMeans(cvLossMatrix)
  
  # calculating sd, cv, up, down
  cvsd <- apply(cvLossMatrix, 2, stats::sd, na.rm = TRUE)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd
  
  
  # finding the minimum cvm value in order pull out the lambda.min out of 
  # list of lambda
  minIndex <- which.min(cvm)    
  lambda.min <- lambdas[minIndex]
  
  # finding 1se index by doing vectorized comparision such that cvm <= cvup 
  # of minIndex
  my1seIndex <- min(which(cvm <= cvup[minIndex]))
  lambda.1se <- lambdas[my1seIndex]
  
  # chosing lambda.final to use by checking lambda.final option
  # note that first element lambda.final default value will return true for
  # lambda.1se 
  chosenLambda <- if(lambda.final == 'lambda.1se')
  {
    lambda.1se
  }else if(lambda.final == 'lambda.min')
  {
    lambda.min
  }
  
  
  
  glmmLasso.final <- glmmLasso::glmmLasso(fix = fix,
                                          rnd = rnd,
                                          data = data,
                                          family = family,
                                          lambda = chosenLambda)
  
  # add control list argument to this to make converge faster form one that 
  # create lambda.1se
  # TODO: (maybe) For final model fit, supply control list from the model that led to     either lambda.1se or lambda.min
  
  # mimicking cv.glmnet return objects
  return_List <- list(lambdas=lambdas,
                      cvm=cvm,
                      cvsd=cvsd,
                      cvup=cvup,
                      cvlo=cvlo,
                      glmmLasso.final=glmmLasso.final,
                      lambda.min=lambda.min,
                      lambda.1se=lambda.1se)
  
  
  class(return_List) <- 'cv.glmmLasso'
  
  
  return(return_List)
  
}



glmmLasso_MultLambdas <- function(fix, rnd, data, 
                                  family = stats::gaussian(link = "identity"), 
                                  lambdas = NULL,
                                  nlambdas = 100,
                                  lambda.min.ratio=ifelse(nobs < nvars, 0.01, 0.0001), 
                                  ...)
{
  
  # fitting first model to generate initial inputs for control parameter
  # here we use the first lambda (highest penalty) to start
  # based glmmLasso's author, glmmLasso is faster when final coefficient
  # estimates corresponding to a lambda is used as the starting value for
  # the next smaller lambda  
  
  # defining the number of observation
  nobs <- nrow(data)
  
  # defining the number of preditors based on the number of terms in fix formula
  nvars <- length(attr(stats::terms(fix), 'term.labels'))
  
  if (is.null(lambdas))
  {
    
    # building the lambda vector
    lambdas <- buildLambdas(fix = fix,
                            rnd = rnd,
                            data = data, 
                            nlambdas = nlambdas, 
                            lambda.min.ratio = lambda.min.ratio)    
  }
  
  
  
  # passing Q.start and Delta.start is modeled from glmmLasso demo file
  # from the "More Elegant section" 
  
  # Delta is matrix containing the estimates of fixed and random effects 
  # (columns) for each iteration (rows) of the main algorithm (i.e. before 
  # the final re-estimation step is performed, see details).
  # Passing the set of estimates from the last iteration as the 
  # 'start' parameter of the controlList
  
  # Q_long is a list containing the estimates of the random effects
  # variance-covariance parameters for each iteration of the main algorithm.
  # Passing the variance-covaiance matrix as the q_start parameter of
  # the controlList
  
  
  
  # initializing list of object to hold the model outputs 
  modList <- vector(mode = 'list', length = length(lambdas))
  
  
  # fit first lambda
  first_fit <- glmmLasso::glmmLasso(fix = fix,
                                    rnd = rnd,
                                    data = data,
                                    family = family,
                                    lambda = lambdas[1],
                                    ...)
  # builing the first Delta.start, transpose required to make dimension
  
  Delta.start <- first_fit$Deltamatrix[first_fit$conv.step, ] %>% t()
  Q.start <- list()
  Q.start[[1]] <- first_fit$Q_long[[first_fit$conv.step + 1]]
  if(nrow(Q.start[[1]])==1) Q.start[[1]] <- c(Q.start[[1]])
  
  for (l in seq_along(lambdas))
  {
    
    # for showing lambda at each iteration
    # message(sprintf('Lambda: %s\n ', lambdas[l]))
    
    fit <- glmmLasso::glmmLasso(fix = fix,
                                rnd = rnd,
                                data = data,
                                family = family,
                                lambda = lambdas[l],
                                control = list(start=Delta.start[l,],
                                               q_start=Q.start[[l]]),...)
    
    # storing model objects before storing to modList
    fit$lambda <- lambdas[l]
    fit$Delta.start <- Delta.start[l,]
    fit$Q.start <- Q.start[[l]]
    fit$data <- data
    fit$rnd <- rnd
    fit$fix <- fix
    fit$family <- family
    
    modList[[l]] <- fit
    Delta.start <- rbind(Delta.start, fit$Deltamatrix[fit$conv.step, ])
    Q.start[[l+1]] <- fit$Q_long[[fit$conv.step + 1]]
    if(nrow(Q.start[[l+1]])==1) Q.start[[l+1]] <- c(Q.start[[l+1]])
    
    
    
  }
  
  # the function returns a list of glmmLasso models 
  
  attr(modList, 'lambdas') <- lambdas
  
  class(modList) <- 'glmmLasso_MultLambdas'
  
  return(modList)
}



predict.glmmLasso_MultLambdas <- function(object, newdata, ...)
{
  # instantiating list to hold nlambdas number of n x 1 vectors 
  # pred_vec_list <- vector(mode = 'list', length = length(object))
  # storing returned vectors in a list 
  
  pred_vec_list <- purrr::map(.x = object, .f = stats::predict, 
                              newdata = newdata)
  
  pred_matrix <- do.call(what = cbind, args = pred_vec_list)
  
  return(pred_matrix)
}
####

calc_mse <- function(actual, predicted)
{
  return(colMeans((actual - predicted)^2)) 
}

#' 
#' @title calc_logloss
#' @description Functions for calculating logloss
#' @details Loss functions written for use in cv.glmmLasso 
#' @author Pirapong Jitngamplang, Jared Lander
#' @param actual actual data values 
#' @param predicted predicted data values
#' @return error between actual versus prediction
#'
calc_logloss <- function(actual, predicted)
{
  
  score <- -(actual * log(predicted) + (1 - actual) * log(1 -predicted))
  score[actual == predicted] <- 0
  score[is.nan(score)] <- Inf
  return(colMeans(score))
  
}

#' 
#' @title calc_multilogloss
#' @description Function for calculating multilogloss
#' @details loss functions written for use in cv.glmmLasso 
#' @author Pirapong Jitngamplang, Jared Lander
#' @param actual actual data values 
#' @param predicted predicted data values
#' @return error between actual versus prediction
#'

# modified from MultiLogLoss in MLMetrics package - credit to Yachen Yan

calc_multilogloss <- function(actual, predicted) 
{
  return(apply(predicted, 2, MLmetrics::MultiLogLoss, y_true = actual)) 
}


#' 
#' @title calc_deviance
#' @description Functions for calculating deviance
#' @details loss functions written for use in cv.glmmLasso 
#' @author Pirapong Jitngamplang, Jared Lander
#' @param actual actual data values 
#' @param predicted predicted data values
#' @param family default value is poisson
#' @param \dots can receive parameters accepted by dismo::calc.deviance
#' @return error between actual versus prediction
#'
calc_deviance <- function(actual, predicted, family = 'poisson',...)
{
  
  return(apply(predicted, 2, dismo::calc.deviance, obs = actual, family = family,
               ...))
}



####


computeLambdaMax <- function(fix, rnd, data, scale=TRUE)
{
  # converting formula into matrices to do lambdaMax calculation
  y <- useful::build.y(fix, data)
  x <- useful::build.x(fix, data)
  
  if(scale)
  {
    x <- scale(x)
  }
  
  # exp because of log scale
  # N*alpha*lambdaMax = max_l(<x_l, y>)
  lambdaMax <- exp(max(abs(colSums(x*y)), na.rm=TRUE) / nrow(data))
  
  # colSums(x*y) is same as crossprod(x,y)
  
  return(lambdaMax)
}

#' @title buildLambdas
#' @description generate lambda vector based on dataset given
#' @author Pirapong Jitngamplang, Jared Lander
#' @param fix A two-sided linear formula object describing the fixed-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. For categorical covariables use as.factor(.) in the formula. Note, that the corresponding dummies are treated as a group and are updated blockwise
#' @param rnd A two-sided linear formula object describing the random-effects part of the model, with the grouping factor on the left of a ~ operator and the random terms, separated by + operators, on the right; alternatively, the random effects design matrix can be given directly (with suitable column names). If set to NULL, no random effects are included.
#' @param data The data frame containing the variables named in formula.
#' @param nlambdas the number of lambdas values, default value is 100 if lambdas is not user-supplied
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size nobs relative to the number of variables nvars. If nobs > nvars, the default is 0.0001, close to zero. If nobs < nvars, the default is 0.01.
#' @return returns a vector of lambda
#'


buildLambdas <- function(fix, rnd, data, 
                         nlambdas = 100, 
                         lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001))
{
  # converting formula into matrices to do lambdaMax calculation
  x <- useful::build.x(fix, data)
  nobs <- nrow(x)
  nvars <- ncol(x)
  
  lambdaMax = computeLambdaMax(fix = fix, 
                               rnd = rnd,
                               data = data)
  
  lambda_vec <- seq(from = lambdaMax, 
                    to = lambdaMax * lambda.min.ratio, 
                    length.out = nlambdas) 
  # sorting such that first lambda is the largest
  lambda_vec <- sort(lambda_vec, decreasing = TRUE)
  
  return(lambda_vec)
}