# This code was built using the following packages and versions
# R version 3.6.1 (2019-07-05)
# splines (version: 3.6.1)
# nlme (version: 3.1.140)
# foreign (version: 0.8.71)

#################
# Load Packages #
#################

library(nlme)
library(splines)
library(foreign)

#################
# Load Data set #
#################


data <- suppressWarnings(read.spss(file.path("Long format met NIHSS sens en neglect.sav"), 
                                   to.data.frame = TRUE, use.value.labels = TRUE))

data <- data[order(data$Number, data$Days),]

data$FE <- factor(data$FE, levels = c(0,1,2), labels = c("none","partial","full"))
data$PREFERRED_HAND <- factor(data$PREFERRED_HAND, levels = c("right", "left", "no clear preference"), labels = c("right", "left", "no clear preference"))
data$RTPA[data$RTPA == 4] <- NA
data$RTPA <- droplevels(data$RTPA)

data$Sens[data$Sens == 10] <- NA
data$Sens[data$Sens == 3] <- NA
data$Sens <- droplevels(data$Sens)

data <- data[complete.cases(data), ]

## check patient 123 is duplicated
drops <- "Index1"
dataN <- data[ , !(names(data) %in% drops)] 
id_m <- which(duplicated(dataN) | duplicated(dataN[nrow(data):1, ])[nrow(dataN):1])

data <- data[-c(id_m),]

data.id <- data[!duplicated(data$Number), ]
data.idTAIL <- data[tapply(row.names(data), data$Number, tail, 1), ]

##### exclude patients with 1 row
##### remove 3 patients (21 observations) with PREFERRED_HAND == no clear preference

table(table(data$Number))

dat1mes <- lapply(split(data, data$Number), function(x) as.data.frame(x))
iDs <- lapply(dat1mes, function(x) if (dim(x)[1] == 1) print(x$Number))
vecIDS <- do.call(c, iDs)

data2 <-  data[!data$Number %in% vecIDS, ]
data2.id <- data2[!duplicated(data2$Number), ]

data2 <- data2[data2$PREFERRED_HAND != "no clear preference", ]
data2.id <- data2.id[data2.id$PREFERRED_HAND != "no clear preference", ]

data2$PREFERRED_HAND <- droplevels(data2$PREFERRED_HAND)
data2.id$PREFERRED_HAND <- droplevels(data2.id$PREFERRED_HAND)

#############################
# Cross Validation Function #
#############################

CV_MSE <- function(i){ 
  
  ### Packages
  library(nlme)
  library(splines)
  
  ### Functions
  IndvPred_lme <- function(lmeObject, newdata, timeVar, times = NULL, M = 200L,
                           interval = c("confidence", "prediction"),
                           all_times = FALSE,
                           level = 0.95, return_data = FALSE, seed = 1L) {
    if (!inherits(lmeObject, "lme") && !inherits(lmeObject, "lmeComponents"))
      stop("Use only with 'lme' or 'lmeComponents' objects.\n")
    interval <- match.arg(interval)
    if (inherits(lmeObject, "lme")) {
      data <- lmeObject$data
      formYx <- formula(lmeObject)
      mfX <- model.frame(terms(formYx), data = data)
      TermsX <- attr(mfX, "terms")
      formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
      mfZ <- model.frame(terms(formYz), data = data)
      TermsZ <- attr(mfZ, "terms")
      idVar <- names(lmeObject$modelStruct$reStruct)
      betas <- fixef(lmeObject)
      sigma <- lmeObject$sigma
      D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", sigma^2)[[1]]
      V <- vcov(lmeObject)
      times_orig <- data[[timeVar]]
      times_orig <- times_orig[!is.na(times_orig)]
    } else {
      formYx <- lmeObject$formYx
      TermsX <- lmeObject$TermsX
      formYz <- lmeObject$formYz
      TermsZ <- lmeObject$TermsZ 
      idVar <- lmeObject$idVar
      betas <- lmeObject$betas
      sigma <- lmeObject$sigma
      D <- lmeObject$D
      V <- lmeObject$V
      times_orig <- lmeObject$times_orig
    }
    # drop missing values from newdata
    all_vars <- unique(c(all.vars(TermsX), all.vars(TermsZ)))
    newdata_nomiss <- newdata[complete.cases(newdata[all_vars]), ]
    mfX_new <- model.frame(TermsX, data = newdata_nomiss)
    X_new <- model.matrix(formYx, mfX_new)
    mfZ_new <- model.frame(TermsZ, data = newdata_nomiss)
    Z_new <- model.matrix(formYz, mfZ_new)
    na_ind <- attr(mfX_new, "na.action")
    y_new <- model.response(mfX_new, "numeric")
    if (length(idVar) > 1)
      stop("the current version of the function only works with a single grouping variable.\n")
    if (is.null(newdata[[idVar]]))
      stop("subject id variable not in newdata.")
    id_nomiss <- match(newdata_nomiss[[idVar]], unique(newdata_nomiss[[idVar]]))
    n <- length(unique(id_nomiss))
    ######################################################################################
    modes <- matrix(0.0, n, ncol(Z_new))
    post_vars <- DZtVinv <- vector("list", n)
    for (i in seq_len(n)) {
      id_i <- id_nomiss == i
      X_new_id <- X_new[id_i, , drop = FALSE]
      Z_new_id <- Z_new[id_i, , drop = FALSE]
      Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + sigma^2 * diag(sum(id_i)))
      DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
      modes[i, ] <- c(DZtVinv[[i]] %*% (y_new[id_i] - X_new_id %*% betas))
      t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
      t2 <- DZtVinv[[i]] %*% X_new_id %*% V %*% 
        crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
      post_vars[[i]] <- D - t1 + t2
    }
    fitted_y <- c(X_new %*% betas) + rowSums(Z_new * modes[id_nomiss, , drop = FALSE])
    ######################################################################################
    if (is.null(times) || !is.numeric(times)) {
      times <- seq(min(times_orig), max(times_orig), length.out = 100)
    }
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    last_time <- tapply(newdata[[timeVar]], id, max)
    times_to_pred <- lapply(last_time, function (t) 
      if (all_times) times else times[times > t])
    id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
    #newdata_pred <- newdata_pred[id_pred, ]
    newdata_pred <- right_rows(newdata, newdata[[timeVar]], id, times_to_pred)
    newdata_pred[[timeVar]] <- unlist(times_to_pred)
    mfX_new_pred <- model.frame(TermsX, data = newdata_pred, na.action = NULL)
    X_new_pred <- model.matrix(formYx, mfX_new_pred)
    mfZ_new_pred <- model.frame(TermsZ, data = newdata_pred, na.action = NULL)
    Z_new_pred <- model.matrix(formYz, mfZ_new_pred)
    predicted_y <- c(X_new_pred %*% betas) + 
      rowSums(Z_new_pred * modes[id_pred, , drop = FALSE])
    set.seed(seed)
    betas_M <- MASS::mvrnorm(M, betas, V)
    modes_fun <- function (betas) {
      t(mapply("%*%", DZtVinv, split(y_new - X_new %*% betas, id_nomiss)))
    }
    modes_M <- lapply(split(betas_M, row(betas_M)), modes_fun)
    matrix_row <- function (m, i) m[i, , drop = FALSE]
    modes_M <- lapply(seq_len(n), function (i) t(sapply(modes_M, matrix_row, i = i)))
    b_M <- modes_M
    for (i in seq_len(n)) {
      b_M[[i]] <- t(apply(modes_M[[i]], 1, MASS::mvrnorm, n = 1, Sigma = post_vars[[i]]))
    }
    n_pred <- length(predicted_y)
    sampled_y <- matrix(0.0, n_pred, M)
    for (m in seq_len(M)) {
      betas_m <- betas_M[m, ]
      b_m <- t(sapply(b_M, function (x) x[m, ]))
      mean_m <- c(X_new_pred %*% betas_m) + 
        rowSums(Z_new_pred * b_m[id_pred, , drop = FALSE])
      sampled_y[, m] <- if (interval == "confidence") mean_m 
      else rnorm(n_pred, mean_m, lmeObject$sigma)
    }
    low <- apply(sampled_y, 1, quantile, probs = (1 - level) / 2)
    upp <- apply(sampled_y, 1, quantile, probs = 1 - (1 - level) / 2)
    rm(list = ".Random.seed", envir = globalenv())
    if (!return_data) {
      list(times_to_pred = times_to_pred, predicted_y = predicted_y, 
           low = low, upp = upp)
    } else {
      out_data <- rbind(newdata, newdata_pred)
      out_data$pred <- c(fitted_y, predicted_y)
      out_data$low <- c(rep(NA, length(fitted_y)), low)
      out_data$upp <- c(rep(NA, length(fitted_y)), upp)
      out_data[order(out_data[[idVar]], out_data[[timeVar]]), ]
    }
  }
  
  right_rows <- function(data, times, ids, Q_points) {
    fids <- factor(ids, levels = unique(ids))
    if (!is.list(Q_points))
      Q_points <- split(Q_points, row(Q_points))
    ind <- mapply(findInterval, Q_points, split(times, fids))
    ind[ind < 1] <- 1
    rownams_id <- split(row.names(data), fids)
    ind <- matrix(ind, ncol = 1)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data[c(ind), ]
  }
  
  ### Split data sets
  trainingData <- data2[!data2$Number %in% i, ]
  testingData <- data2[data2$Number %in% i, ]
  
  ### Models
  model_noCov <- lme(ARAT ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)), 
                     data = trainingData,
                     random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), 
                     na.action = na.exclude,
                     control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  
  model_sigAllMain_sigAllInter <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                        SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                        MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                        MILEG + 
                                        FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                        FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                        Neglect, data = trainingData,
                                      random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                                      control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  
  model_sigCov_minRTPA <- lme(ARAT ~  #RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                MILEG + 
                                FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                Neglect, data = trainingData,
                              random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                              control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_sigCov_minFE <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              MILEG + 
                              #FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              Neglect, data = trainingData,
                            random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                            control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_sigCov_minSA <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              #SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              MILEG + 
                              FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                              Neglect, data = trainingData,
                            random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                            control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_sigCov_minMIARM <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 #MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 MILEG + 
                                 FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 Neglect, data = trainingData,
                               random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_sigCov_minFMARM <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 MILEG + 
                                 FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 #FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 Neglect, data = trainingData,
                               random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_sigCov_minMILEG <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 #MILEG + 
                                 FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 Neglect, data = trainingData,
                               random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))  
  
  model_sigCov_minNeglect <- lme(ARAT ~  RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                   SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                   MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                   MILEG + 
                                   FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                   FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)), data = trainingData,
                                 random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                                 control = lmeControl(maxIter = 1e8, msMaxIter = 1e8)) 
  
  model_SA_FE_withInter <- lme(ARAT ~  SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                 FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)), data = trainingData,
                               random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_SA_FE_noInter <- lme(ARAT ~   ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                               FE + SA, data = trainingData,
                             random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                             control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_SA_FE_Neglect_MILEG_withInter <- lme(ARAT ~  SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                               FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) +
                                               Neglect + MILEG, data = trainingData,
                                             random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude,
                                             control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  model_allMain_allInter <- lme(ARAT ~ AGE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + GENDER * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) +
                                  BAMFORD * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) +
                                  AFFECTED_BODYSIDE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + PREFERRED_HAND * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  +
                                  CIRSTOT * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  +
                                  SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  +
                                  MILEG * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) +
                                  FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + NIHSS * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  +
                                  Sens * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + Neglect * ns(Days, knots = c(6, 11, 19, 34, 50, 91)), data = trainingData,
                                random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude)
  
  
  model_allMain_sigAllInter <- lme(ARAT ~ AGE + GENDER + 
                                     BAMFORD + RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                     AFFECTED_BODYSIDE + PREFERRED_HAND  + 
                                     CIRSTOT  + 
                                     SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                     MILEG + FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                                     FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + NIHSS  + 
                                     Sens + Neglect, data = trainingData,
                                   random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude)
  
  
  
  ### Predictions
  ##### function to calculate predictions
  
  pred_fun <- function(model, method){
    k = 1
    dat_pred1 <- list()
    for (o in unique(testingData$Number)) {
      newdata <- testingData[testingData$Number == o,]
      for (j in 1:(length(testingData[["Days"]][testingData$Number == o])-1)) {
        dat_pred1[[k]] <- newdata[-(1:j),]
        k <- k + 1
      }
    }
    
    k = 1
    dat_pred1USED <- list()
    for (o in unique(testingData$Number)) {
      newdataUSED <- testingData[testingData$Number == o,]
      for (j in 1:(length(testingData[["Days"]][testingData$Number == o])-1)) {
        dat_pred1USED[[k]] <- newdataUSED[(1:j),]
        k <- k + 1
      }
    }
    
    
    times1 <- lapply(dat_pred1, FUN = function(x) x[["Days"]])
    try(pred1 <- mapply(FUN = function(x, y) {
      IndvPred_lme(lmeObject = model, newdata = x, timeVar = "Days", times = y,
                   M = 200L,
                   interval = c("confidence"),
                   all_times = TRUE,
                   level = 0.95, return_data = FALSE, seed = 1L)$predicted_y
    }, dat_pred1USED, times1, SIMPLIFY = FALSE))
    
    dataRealTime1 <- lapply(dat_pred1, function(x) x$Days)
    dataRealARAT1 <- lapply(dat_pred1, function(x) x$ARAT)  
    dataID1 <- lapply(dat_pred1, function(x) x$Number)
    
    usedMeas1 <- lapply(dat_pred1USED, function(x) nrow(x))
    dataLastTime1 <- lapply(dat_pred1USED, function(x) tail(x$Days, n = 1))
    
    predictionMat1 <- do.call(rbind, mapply(cbind, dataID1, dataRealTime1, pred1, dataLastTime1, dataRealARAT1, usedMeas1, method))
    predictionMat1 <- as.data.frame(predictionMat1, stringsAsFactors = FALSE)
    colnames(predictionMat1) <- c("ID", "Pred_time", "Pred_ARAT", "last_time", "real_ARAT", "Measurements_used", "Method")
    predictionMat1
  }   
  
  # model_noCov
  predictionMat1 <- pred_fun(model = model_noCov, method = 1)
  
  # model_sigAllMain_sigAllInter
  predictionMat2 <- pred_fun(model = model_sigAllMain_sigAllInter, method = 2)
  
  # model_sigCov_minRTPA
  predictionMat3 <- pred_fun(model = model_sigCov_minRTPA, method = 3)
  
  # model_sigCov_minFE
  predictionMat4 <- pred_fun(model = model_sigCov_minFE, method = 4)
  
  # model_sigCov_minSA
  predictionMat5 <- pred_fun(model = model_sigCov_minSA, method = 5)
  
  # model_sigCov_minMIARM
  predictionMat6 <- pred_fun(model = model_sigCov_minMIARM, method = 6)
  
  # model_sigCov_minFMARM
  predictionMat7 <- pred_fun(model = model_sigCov_minFMARM, method = 7)
  
  # model_sigCov_minMILEG
  predictionMat8 <- pred_fun(model = model_sigCov_minMILEG, method = 8)
  
  # model_sigCov_minNeglect
  predictionMat9 <- pred_fun(model = model_sigCov_minNeglect, method = 9)
  
  # model_SA_FE_withInter
  predictionMat10 <- pred_fun(model = model_SA_FE_withInter, method = 10)
  
  # model_SA_FE_Neglect_MILEG_withInter
  predictionMat11 <- pred_fun(model = model_SA_FE_Neglect_MILEG_withInter, method = 11)
  
  # model_allMain_allInter
  predictionMat12 <- pred_fun(model = model_allMain_allInter, method = 12)
  
  # model_allMain_noAllInter
  predictionMat13 <- pred_fun(model = model_allMain_sigAllInter, method = 13)
  
  # model_SA_FE_noInter
  predictionMat14 <- pred_fun(model = model_SA_FE_noInter, method = 14)
  
  ### Combine results
  rbind(predictionMat1, predictionMat2, predictionMat3, predictionMat4, predictionMat5, 
        predictionMat6, predictionMat7, predictionMat8, predictionMat9, predictionMat10, 
        predictionMat11, predictionMat12, predictionMat13, predictionMat14)
}




################################
# Run multiple (50) iterations #
################################
# In order to use multiple cores in your computer change makeCluster(1)


M = 10
Res <- list()
set.seed(123)
for (p in 1:M) {
  V <- 5
  n <- nrow(data2.id)
  splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
  
  cl <- makeCluster(4)
  clusterExport(cl, c("data2"))
  Res[[p]] <- try(parLapply(cl, splits, CV_MSE))
  stopCluster(cl)
  print(p)
}



####################
# Save the results #
####################
# The cross validation function returns a list

save(Res, file = "CV_res_new4.RData")


#####################################
# Merge results from all iterations #
#####################################
# The cross validation function returns a list - the following code with combine the lists

load("CV_res_new4.RData")

M = 10
CV_res <- list()
for (k in c(1:M)){
  CV_res[[k]] <- do.call("rbind", Res[[k]])
}


CV_res_final <- do.call(rbind, CV_res)



CV_res_final <- do.call(rbind, CV_res)


CV_res_final$Method <- factor(CV_res_final$Method, levels = c(0:4), 
                              labels = c("Sc0", "Sc1", "Sc2", "Sc3", "Sc4"))


# Set all predictions above 16 set to 16
# Correct for ceiling and floor effect

CV_res_final$Pred_ARAT[CV_res_final$Pred_ARAT < 0] <- 0
CV_res_final$Pred_ARAT[CV_res_final$Pred_ARAT > 57] <- 57



####################
# Plot the results #
####################

# Ontain bias
CV_res_final$bias <- abs(unlist(CV_res_final$"Pred_BI") - unlist(CV_res_final$real_BI))CV_res_final$AE <- abs(CV_res_final$Pred_ARAT - CV_res_final$real_ARAT)

# Plot
boxplot(CV_res_final$AE ~ CV_res_final$Method, ylim = c(0,35), ylab = "Absolute Error")