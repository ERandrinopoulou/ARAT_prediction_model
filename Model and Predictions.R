# This code was built using the following packages and versions
# R version 3.6.1 (2019-07-05)
# splines (version: 3.6.1)
# nlme (version: 3.1.140)
# foreign (version: 0.8.71)
# ggplot2 (version: 3.2.1)

#################
# Load Packages #
#################

library(nlme)
library(splines)
library(foreign)
library(ggplot2)

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

##################################################################################
############################# Investigating linearity ############################
##################################################################################


model_nsCHECK1 <- lme(ARAT ~ ns(Days, df = 2), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, df = 2))), na.action = na.exclude)
model_nsCHECK2 <- lme(ARAT ~ ns(Days, df = 3), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, df = 3))), na.action = na.exclude)
model_nsCHECK3 <- lme(ARAT ~ ns(Days, df = 4), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, df = 4))), na.action = na.exclude)
model_nsCHECK4 <- lme(ARAT ~ ns(Days, knots = c(6, 11, 19, 34, 50, 91)), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude)

set.seed(123)
set.seed(1234)
ind <- sample(unique(data$Number), 49, replace = F)

data$fitted_marg <- fitted(model_nsCHECK1, level = 0)
data$fitted_subj <- fitted(model_nsCHECK1, level = 1)
p1 <- xyplot(ARAT + fitted_marg + fitted_subj ~ Days | Number, data = data,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 3)
               y.mat <- matrix(y, ncol = 3)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "blue")
             }, subset = Number %in% ind, layout = c(7, 7), as.table = TRUE,
             xlab = "Time", ylab = "ARAT")

data$fitted_marg <- fitted(model_nsCHECK2, level = 0)
data$fitted_subj <- fitted(model_nsCHECK2, level = 1)
p2 <- xyplot(ARAT + fitted_marg + fitted_subj ~ Days | Number, data = data,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 3)
               y.mat <- matrix(y, ncol = 3)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "red")
             }, subset = Number %in% ind, layout = c(7, 7), as.table = TRUE,
             xlab = "Time", ylab = "ARAT")


data$fitted_marg <- fitted(model_nsCHECK3, level = 0)
data$fitted_subj <- fitted(model_nsCHECK3, level = 1)
p3 <- xyplot(ARAT + fitted_marg + fitted_subj ~ Days | Number, data = data,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 3)
               y.mat <- matrix(y, ncol = 3)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "green")
             }, subset = Number %in% ind, layout = c(7, 7), as.table = TRUE,
             xlab = "Time", ylab = "ARAT")


data$fitted_marg <- fitted(model_nsCHECK4, level = 0)
data$fitted_subj <- fitted(model_nsCHECK4, level = 1)
p4 <- xyplot(ARAT + fitted_marg + fitted_subj ~ Days | Number, data = data,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 3)
               y.mat <- matrix(y, ncol = 3)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "black")
             }, subset = Number %in% ind, layout = c(7, 7), as.table = TRUE,
             xlab = "Time", ylab = "ARAT")

AIC(model_nsCHECK1, model_nsCHECK2, model_nsCHECK3, model_nsCHECK4)
p1 + as.layer(p2) + as.layer(p3) + as.layer(p4)


##################################################################################
############## Final model investigating covariates and interactions #############
##################################################################################


full.model <- lme(ARAT ~ AGE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + GENDER * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                    BAMFORD * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + RTPA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                    AFFECTED_BODYSIDE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + PREFERRED_HAND * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + 
                    CIRSTOT * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + 
                    SA * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + MIARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + 
                    MILEG * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + FE * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + 
                    FMARM * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + NIHSS * ns(Days, knots = c(6, 11, 19, 34, 50, 91))  + 
                    Sens * ns(Days, knots = c(6, 11, 19, 34, 50, 91)) + Neglect * ns(Days, knots = c(6, 11, 19, 34, 50, 91)), data = data, 
                  random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude)
anova(full.model, type = "marginal")

####

full.model <- lme(ARAT ~ AGE * ns(Days, df = 3) + GENDER * ns(Days, df = 3) + 
                    BAMFORD * ns(Days, df = 3)  + RTPA * ns(Days, df = 3) + 
                    AFFECTED_BODYSIDE * ns(Days, df = 3) + PREFERRED_HAND * ns(Days, df = 3)  + 
                    CIRSTOT * ns(Days, df = 3)  + 
                    SA * ns(Days, df = 3) + MIARM * ns(Days, df = 3)  + 
                    MILEG * ns(Days, df = 3)  + FE * ns(Days, df = 3) + 
                    FMARM * ns(Days, df = 3) + NIHSS * ns(Days, df = 3)  + 
                    Sens * ns(Days, df = 3) + Neglect * ns(Days, df = 3), data = data, 
                  random = list(Number = pdDiag(form = ~  ns(Days, df = 3))), na.action = na.exclude)
anova(full.model, type = "marginal")

opt.model1 <- lme(ARAT ~ AGE + GENDER + 
                    BAMFORD + RTPA * ns(Days, df = 3) + 
                    AFFECTED_BODYSIDE + PREFERRED_HAND  + 
                    CIRSTOT  + 
                    SA * ns(Days, df = 3) + MIARM * ns(Days, df = 3) + 
                    MILEG + FE * ns(Days, df = 3) + 
                    FMARM * ns(Days, df = 3) + NIHSS  + 
                    Sens + Neglect, data = data, 
                  random = list(Number = pdDiag(form = ~  ns(Days, df = 3))), na.action = na.exclude)
anova(opt.model1, type = "marginal")

opt.model2 <- lme(ARAT ~  RTPA * ns(Days, df = 3) + 
                    SA * ns(Days, df = 3) + MIARM * ns(Days, df = 3) + 
                    MILEG + FE * ns(Days, df = 3) + 
                    FMARM * ns(Days, df = 3) + 
                    Neglect, data = data, 
                  random = list(Number = pdDiag(form = ~  ns(Days, df = 3))), na.action = na.exclude)
anova(opt.model2, type = "marginal")


anova(update(full.model,method="ML"),update(full.model,~.-BAMFORD * ns(Days, df = 3),method="ML"))


anova(full.model, type = "marginal")
anova(update(full.model, method="ML"),update(full.model,~.-RTPA * ns(Days, df = 3), method="ML"))


anova(full.model, opt.model1, opt.model2)


#############################
# Fit a Mixed-effects Model #
#############################

final.model1 <- lme(ARAT ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)) * SA + 
                      ns(Days, knots = c(6, 11, 19, 34, 50, 91)) * FE, 
                    data = data[!data$Number %in% indW, ],
                    random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))),
                    na.action = na.exclude,
                    control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))


FinalRes <- summary(final.model1)$tTable

ind1 <- grep(":",rownames(FinalRes))
ind2 <- grep("ns",rownames(FinalRes))
ind3 <- grep("(Intercept)",rownames(FinalRes))

varKeepF <- rownames(FinalRes)[-c(ind1, ind2, ind3)]


tableF <- summary(final.model1)$tTable
vecKeepF <-  tableF[c(2), ]
for (p in 2:length(varKeepF)) {
  ii <- rownames(FinalRes)[rownames(FinalRes) %in% varKeepF[p]]  #grep(covarKeep[p], rownames(tableF))
  vecKeepF <- rbind(vecKeepF, FinalRes[c(ii), ])
}

final_modelF <- data.frame(names = varKeepF, vecKeepF)
final_modelF <- final_modelF[,c(1,2,3,6)]


#############################
# Functions for Predictions #
#############################

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

DynPlots <- function(model.output = model.output, newdata, timeVar, 
                        main_title = "Dynamic predictions", legend = TRUE,
                        nameX = "time since stroke (days)",
                        nameY = "ARAT") {
  
  
  # Generating individual prediction ------------------------------------
  
  data <- model.output$data
  formYx <- formula(model.output)
  y <- formYx[[2]]
  
  IndvPrediction <- IndvPred_lme(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                 interval = "prediction", return_data = TRUE)
  
  pred <- IndvPrediction[which(!is.na(IndvPrediction$low)),]
  
  # Generating plot -----------------------------------------------------
  plot <- ggplot() + theme_bw()+
    geom_point(aes(x = newdata[[timeVar]], y = newdata[[y]], colour = "Obsv"), size = 1.5) +
    
    geom_line(aes(x = pred[[timeVar]], y = pred[["pred"]], colour = "Pred"), size = 1.2) +
    geom_line(aes(x = pred[[timeVar]], y = pred[["low"]], colour = "CI"), linetype = 3, size = 1.2) +
    geom_line(aes(x = pred[[timeVar]], y = pred[["upp"]], colour = "CI"), linetype = 3, size = 1.2) +
    geom_vline(xintercept = tail(newdata[[timeVar]], n = 1), linetype = "longdash", colour = c("black")) +
    
    scale_x_continuous(name = nameX, limits = c(0, 210), expand = c(0,0)) +
    scale_y_continuous(name = nameY, limits = c(0, 60), expand = c(0,0)) +
    scale_colour_manual(name = "", values = c("black", "black", "black")) + 
    
    ggtitle(main_title) +
    theme(
      plot.title = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.position = "none"
    )  
  plot
}

###################################################
# Dynamic Prediction Illustration for Patient 120 #
###################################################

newPatient <- data[data$Number == 120, ]

newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:4, ]

DynPlots(model.output = final.model1, newdata = newPatient1, 
                   timeVar = "Days",
                   main_title = "",
                   nameX = "", nameY = "ARAT")
DynPlots(model.output = final.model1, newdata = newPatient2, 
                   timeVar = "Days",
                   main_title = "")