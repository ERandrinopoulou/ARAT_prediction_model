#################
# Load packaged #
#################

library(shiny)
library(lattice)
library(latticeExtra)
library(ggplot2)
library(gridExtra)
library(grid)
library(JMbayes)
library(splines)
library(lattice)

#########################
# Load data and results #
#########################

load("CV_res_new4.RData")
load("data.RData")

M = 10
CV_res <- list()
for (k in c(1:M)){
  CV_res[[k]] <- do.call("rbind", Res[[k]])
}


CV_res_final <- do.call(rbind, CV_res)



CV_res_final$AE <-  abs(CV_res_final$Pred_ARAT - CV_res_final$real_ARAT)

CV_res_final$Pred_ARAT[CV_res_final$Pred_ARAT < 0] <- 0
CV_res_final$AE[CV_res_final$AE > 57] <- 57


CV_res_final$Method <- factor(CV_res_final$Method, levels = 1:14, 
                              labels = c("No Cov", "model_sigAllMain_sigAllInter", "model_sigCov_minRTPA", "model_sigCov_minFE", 
                                         "model_sigCov_minSA", "model_sigCov_minMIARM", "model_sigCov_minFMARM", "model_sigCov_minMILEG", 
                                         "model_sigCov_minNeglect",
                                         "model_SA_FE_withInter", 
                                         "model_SA_FE_Neglect_MILEG_withInter", "model_allMain_allInter", "model_allMain_sigAllInter",
                                         "model_SA_FE_noInter"))

##################################
# Checks for decreasing patients #
##################################

indW <- numeric()


for (i in unique(data$Number)) {
  for (j in 1:(length(data$Number[data$Number == i]) - 1)){
    try(if (((data$ARAT[data$Number == i][j])) - ( data$ARAT[data$Number == i][j + 1]) >= 7 ) {
      indW[i] <- i
    })
    
  }
}


indW <- indW[!is.na(indW)]

data$Number <- as.factor(data$Number)
CHECK_pat1 <- xyplot(ARAT ~ Days | Number, data = data[data$Number %in% indW,], ylim = c(-1, 60), 
                     xlim = c(-1,300), type = c("l"), col = "blue",  grid = TRUE, 
                     main = "Patients with decreased evolution.")


data.id <- data[!duplicated(data$Number), ]

############
# Figure 1 #
############
#500 x 500

tiff("Figure1.tiff", width = 5, height = 5, units = 'in', res = 300)


ProfPatARAT1 <- xyplot(ARAT ~ Days, group = Number, data = data, type = c("l"), cex = 0.6,col = "grey", 
                       strip = FALSE, pch = 16, grid = FALSE, lwd = 2, ylim = c(-1,57.4), xlim = c(0,300), 
                       xlab = list("Days since stroke"), ylab = list("ARAT"),  
                       par.settings = list(fontsize = list(text = 13, points = 18)))

subdata1 <- data[data$Number == 3,]
IndProfPatARAT1 <- xyplot(ARAT ~ Days, groups = Number, type = c("p", "l"), pch = 16, cex = 0.6, col = "blue4", lwd = 4, 
                          data = subdata1, ylim = c(-1,57.4), xlim = c(0,300), 
                          par.settings = list(fontsize = list(text = 13, points = 18)))

subdata2 <- data[data$Number == 266,]
IndProfPatARAT2 <- xyplot(ARAT ~ Days, groups = Number, type = c("p", "l"), pch = 16, cex = 0.6, col = "blue", lwd = 4, 
                          data = subdata2, ylim = c(-1,57.4), xlim = c(0,300), 
                          par.settings = list(fontsize = list(text = 13, points = 18)))

#subdata3 <- data[data$Number == 32,]
#subdata3 <- data[data$Number == 39,]
subdata3 <- data[data$Number == 51,]
IndProfPatARAT3 <- xyplot(ARAT ~ Days, groups = Number, type = c("p", "l"), pch = 16, cex = 0.6, col = "coral4", lwd = 4,
                          data = subdata3, ylim = c(-1,57.4), xlim = c(0,300), 
                          par.settings = list(fontsize = list(text = 13, points = 18)))

subdata4 <- data[data$Number == 351,]
IndProfPatARAT4 <- xyplot(ARAT ~ Days, groups = Number, type = c("p", "l"), pch = 16, cex = 0.6, col = "darkcyan", lwd = 4, 
                          data = subdata4, ylim = c(-1,57.4), xlim = c(0,300), 
                          par.settings = list(fontsize = list(text = 13, points = 18)))


ProfPatARAT1 + as.layer(IndProfPatARAT1) + as.layer(IndProfPatARAT2) + as.layer(IndProfPatARAT3) + 
  as.layer(IndProfPatARAT4) 


dev.off()

############
# Figure 4 #
############
model_SA_FE <- lme(ARAT ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)) * SA + 
                           ns(Days, knots = c(6, 11, 19, 34, 50, 91)) * FE, 
                   data = data[!data$Number %in% indW, ],
                   random = list(Number = pdDiag(form = ~  ns(Days, knots = c(6, 11, 19, 34, 50, 91)))),
                   na.action = na.exclude,
                   control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))


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
  ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
  data[c(ind), ]
}


# DynPlotsOLD_noCI <- function(model.output = model.output, newdata, timeVar, 
#                         main_title = "Dynamic predictions", legend = TRUE,
#                         nameX = "time since stroke (days)",
#                         nameY = "ARAT") {
#   
#   
#   # Generating individual prediction ------------------------------------
#   
#   data <- model.output$data
#   formYx <- formula(model.output)
#   y <- formYx[[2]]
#   
#   IndvPrediction <- IndvPred_lme(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
#                                  interval = "prediction", return_data = TRUE)
#   
#   pred <- IndvPrediction[which(!is.na(IndvPrediction$low)),]
#   
#   # Generating plot -----------------------------------------------------
#   plot <- ggplot() +theme_bw()+
#     geom_point(aes(x = newdata[[timeVar]], y = newdata[[y]], colour = "Obsv"), size = 1.5) +
#     
#     geom_line(aes(x = pred[[timeVar]], y = pred[["pred"]], colour = "Pred"), size = 1.2) +
#     #geom_line(aes(x = pred[[timeVar]], y = pred[["low"]], colour = "CI"), linetype = 3, size = 1.2) +
#     #geom_line(aes(x = pred[[timeVar]], y = pred[["upp"]], colour = "CI"), linetype = 3, size = 1.2) +
#     geom_vline(xintercept = tail(newdata[[timeVar]], n = 1), linetype = "longdash", colour = c("black")) +
#     
#     scale_x_continuous(name = nameX, limits = c(0, 180), expand = c(0,0)) +
#     scale_y_continuous(name = nameY, limits = c(0, 60), expand = c(0,0)) +
#     scale_colour_manual(name = "", values = c("black", "black")) + 
#     
#     ggtitle(main_title) +
#     theme(
#       plot.title = element_text(size = 12),
#       axis.title = element_text(size = 15),
#       axis.text = element_text(size = 12),
#       legend.position = "none"
#     )  
#   plot
# }



DynPlots <- function(model.output = model.output, newdata, timeVar, 
                     main_title = "Dynamic predictions"){
  
  # Load individual prediction ------------------------------------
  data <- model.output$data
  formYx <- formula(model.output)
  yOutcome <- formYx[[2]]
  
  IndvPrediction95 <- IndvPred_lme(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                   interval = "prediction", return_data = TRUE)
  
  IndvPrediction68 <- IndvPred_lme(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                   interval = "prediction", return_data = TRUE, level = 0.68)
  
  pred95 <- IndvPrediction95[which(!is.na(IndvPrediction95$low)),]
  pred68 <- IndvPrediction68[which(!is.na(IndvPrediction68$low)),]
  
  nopred <- IndvPrediction95[which(is.na(IndvPrediction95$low)),]
  
  timeVariable <- pred95[[timeVar]]
  
  # Generating plot -----------------------------------------------------
  xyplot(pred ~ timeVariable , main = main_title, data = pred95,
         type = "l", col = rgb(0.6769,0.4447,0.7114, alpha = 1), lty = c(1, 2, 2), lwd = 3,
         ylim = c(0,65), xlim = c(0,200), ylab = list(yOutcome, cex = 1.5), xlab = list("Days since stroke", cex = 1.5),
         scales = list(x = list(cex = 1.3) , y = list(cex = 1.3)),
         panel = function(x, y,  ...) {
           panel.xyplot(x, y, ...)
           panel.polygon(c(pred95[,"Days"], rev(pred95[,"Days"])), 
                         c(pred95[,"upp"], rev(pred95[,"low"])),
                         border = NA,
                         col = rgb(0.6769,0.4447,0.7114, alpha = 0.2))
           panel.polygon(c(pred68[,"Days"], rev(pred68[,"Days"])), 
                         c(pred68[,"upp"], rev(pred68[,"low"])),
                         border = NA,
                         col =rgb(0.6769,0.4447,0.7114, alpha = 0.4))
           panel.points(x = nopred[[timeVar]], y = nopred[[yOutcome]], cex = 1.2, pch = 16, col = "black");
           panel.lines(x = rep(tail(nopred[[timeVar]], n = 1), 200), y = seq(-100, 100, length = 200), col = "grey", lty = 3, lwd = 2)
         })
  
  
}

##################
newPatient <- data[data$Number == 120, ]
newPatient <- data[data$Number == 134, ]
newPatient <- data[data$Number == 3, ]

newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:3, ]

p11 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
            timeVar = "Days",
            main_title = "")#,
            #nameX = "", nameY = "ARAT")
p12 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
            timeVar = "Days",
            main_title = "")#, 
            #nameX = "", nameY = "")
#grid.arrange(p11, p12, nrow = 2)

##################
newPatient <- data[data$Number == 266, ]

newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:4, ]

p21 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                  timeVar = "Days",
                  main_title = "")#, 
                  #nameX = "time since stroke (days)", nameY = "ARAT")
p22 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                  timeVar = "Days",
                  main_title = "")#, 
                  #nameX = "time since stroke (days)", nameY = "")
#grid.arrange(p21, p22, nrow = 2)

##################
newPatient <- data[data$Number == 25, ]
newPatient <- data[data$Number == 51, ]

newPatient1 <- newPatient[1:4, ]
newPatient2 <- newPatient[1:5, ]

p31 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                  timeVar = "Days",
                  main_title = "")#, 
                  #nameX = "", nameY = "ARAT")
p32 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                  timeVar = "Days",
                  main_title = "")#, 
                  #nameX = "", nameY = "")

#grid.arrange(p31, p32, nrow = 2)

##################
newPatient <- data[data$Number == 351, ]


newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:4, ]

p41 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                  timeVar = "Days",
                  main_title = "")#, 
                  #nameX = "time since stroke (days)", nameY = "ARAT")
p42 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                  timeVar = "Days",
                  main_title = "")#, 
                  #nameX = "time since stroke (days)", nameY = "")

#grid.arrange(p41, p42, nrow = 2)

# grid.arrange(p11, p12, p21, p22, p31, p32, p41, p42, nrow = 4, widths=c(1,1))
# grid.text("Patient 1", x = unit(0.5, "npc"), 
#           y = unit(0.99, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
# grid.text("Patient 2", x = unit(0.5, "npc"), 
#           y = unit(0.74, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
# grid.text("Patient 3", x = unit(0.5, "npc"), 
#           y = unit(0.49, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
# grid.text("Patient 4", x = unit(0.5, "npc"), 
#           y = unit(0.24, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
# 
# 
# grid.text("Time point 1", x = unit(0.28, "npc"), 
#           y = unit(0.975, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
# grid.text("Time point 2", x = unit(0.8, "npc"), 
#           y = unit(0.975, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))

#700 x 500

tiff("Figure4.tiff", width = 7, height = 5, units = 'in', res = 300)

grid.arrange(p11, p12, p21, p22, nrow = 2, widths=c(1,1))
grid.text("Patient 1", x = unit(0.5, "npc"),
          y = unit(0.982, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("Patient 2", x = unit(0.5, "npc"),
          y = unit(0.5, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))

grid.text("Time point 1", x = unit(0.28, "npc"),
          y = unit(0.97, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("Time point 2", x = unit(0.8, "npc"),
          y = unit(0.97, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))

dev.off()

######################
# Obtain mean Errors #
######################

Median_measUsed <- function(input.model, input.meas){ 
  CV_res_final_used <- CV_res_final[CV_res_final$Measurements_used >= input.meas & CV_res_final$Pred_time >= 160 &
                                      CV_res_final$Method == input.model, ]
  median(CV_res_final_used$AE[CV_res_final_used$Measurements_used == input.meas])
  #quantile(CV_res_final_used$AE[CV_res_final_used$Measurements_used == input.meas])[c(2,4)]
}

Quar_measUsed <- function(input.model, input.meas){ 
  CV_res_final_used <- CV_res_final[CV_res_final$Measurements_used >= input.meas & CV_res_final$Pred_time >= 160 &
                                      CV_res_final$Method == input.model, ]
  #median(CV_res_final_used$AE[CV_res_final_used$Measurements_used == input.meas])
  quantile(CV_res_final_used$AE[CV_res_final_used$Measurements_used == input.meas])[c(2,4)]
}

medians_allmeth1 <- medians_allmeth2 <- medians_allmeth3 <- medians_allmeth4 <- medians_allmeth5 <- medians_allmeth6 <- medians_allmeth7 <- numeric()

for (i in levels(CV_res_final$Method)){
  medians_allmeth1[i] <- Median_measUsed(i, 1)
}
for (i in levels(CV_res_final$Method)){
  medians_allmeth2[i] <- Median_measUsed(i, 2)
}
for (i in levels(CV_res_final$Method)){
  medians_allmeth3[i] <- Median_measUsed(i, 3)
}
for (i in levels(CV_res_final$Method)){
  medians_allmeth4[i] <- Median_measUsed(i, 4)
}
for (i in levels(CV_res_final$Method)){
  medians_allmeth5[i] <- Median_measUsed(i, 5)
}
for (i in levels(CV_res_final$Method)){
  medians_allmeth6[i] <- Median_measUsed(i, 6)
}
for (i in levels(CV_res_final$Method)){
  medians_allmeth7[i] <- Median_measUsed(i, 7)
}
res_medians <- data.frame(MAE_1meas = medians_allmeth1, MAE_2meas = medians_allmeth2,
           MAE_3meas = medians_allmeth3, MAE_4meas = medians_allmeth4, MAE_5meas = medians_allmeth5,
           MAE_6meas = medians_allmeth6, MAE_7meas = medians_allmeth7)
round(res_medians, digits = 2)
#

quar_allmeth1 <- quar_allmeth2 <- quar_allmeth3 <- quar_allmeth4 <- quar_allmeth5 <-quar_allmeth6 <- quar_allmeth7 <- matrix(NA, length(levels(CV_res_final$Method)), 2)

for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth1[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 1)
}
for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth2[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 2)
}
for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth3[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 3)
}
for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth4[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 4)
}
for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth5[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 5)
}
for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth6[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 6)
}
for (i in 1:length(levels(CV_res_final$Method))){
  quar_allmeth7[i,] <- Quar_measUsed(levels(CV_res_final$Method)[i], 7)
}
res_Quar <- data.frame(MAE_1meas = quar_allmeth1, MAE_2meas = quar_allmeth2,
                          MAE_3meas = quar_allmeth3, MAE_4meas = quar_allmeth4, 
                          MAE_5meas = quar_allmeth5,
                          MAE_6meas = quar_allmeth6, MAE_7meas = quar_allmeth7)
round(res_Quar, digits = 2)
rownames(res_Quar) <- levels(CV_res_final$Method)


############
# Figure 2 #
############
#500 x 500

tiff("Figure2.tiff", width = 5, height = 5, units = 'in', res = 300)

## measuremenents used
Plot_measUsed <- function(input.model){ 
  CV_res_final_used <- CV_res_final[CV_res_final$Measurements_used >= 1 & CV_res_final$Pred_time >= 160 &
                                      CV_res_final$Method == input.model, ]
  #par(mar = c(4, 5, 2, 4))
  boxplot(CV_res_final_used$AE ~ CV_res_final_used$Measurements_used, horizontal=FALSE, 
          cex.axis = 0.9, ylab = "Absolute Error",
          cex.lab = 1.1, xlab = "Number of serial measurements used", main = "Mixed-Effects Model", 
          cex.main = 1.3, ylim = c(0, 60),
          outline = FALSE)
}

#Plot_measUsed("model_SA_FE_withInter")
Plot_measUsed("model_SA_FE_noInter")

dev.off()


############
# Figure 3 #
############
#900 x 500

Plot_lasttime_perARATbase <- function(modelselection, baselineARAT, title_main){
  
  CV_res_final_used <- CV_res_final[CV_res_final$Measurements_used >= 1 & CV_res_final$Pred_time >= 160, ]
  
  data.id$indPatient <- as.numeric(cut(data.id$ARAT, c(0, 22, 47, 70), right = TRUE, include.lowest = TRUE) )
  data.id$indPatient <- factor(data.id$indPatient, levels = c(1,2,3),
                               labels = c("[0-22]", "(22-47]", "(47-57]"))
  indPatF <- data.id$Number[data.id$indPatient == baselineARAT]
  

  bb <- CV_res_final_used[CV_res_final_used$Method == modelselection,]
  bb$group[bb$last_time >=0 & bb$last_time <5] <- 1
  bb$group[bb$last_time >=5 & bb$last_time <11] <- 2
  bb$group[bb$last_time >=11 & bb$last_time <22] <- 3
  bb$group[bb$last_time >=22 & bb$last_time <60] <- 4
  bb$group[bb$last_time >=60 & bb$last_time <160] <- 5
  bb$group <- factor(bb$group, labels = c("[0-5)", "[5-11)", "[11-22)", "[22-60)", "[60-159)"))
  
  boxplot(bb$AE[bb$ID %in% indPatF] ~ bb$group[bb$ID %in% indPatF], horizontal=FALSE, ylab = "", 
          cex.lab = 2.5, xlab = "", main = title_main, 
          cex.main = 2, ylim = c(0, 60), cex.axis = 1.2, outline = FALSE, las = 2)
}

tiff("Figure3.tiff", width = 9, height = 5, units = 'in', res = 300)

par(mgp=c(7.5,0.9,0))
par(mar=c(7.5,4.8,2,2))
par(mfrow = c(1,3))
Plot_lasttime_perARATbase("model_SA_FE_withInter", "(47-57]", title_main = "Baseline ARAT 48-57")
mtext(side = 1, text="", line = 2.2)
mtext(side = 2, text="Absolute error", line = 3, cex = 1.3)

Plot_lasttime_perARATbase("model_SA_FE_withInter", "(22-47]", title_main = "Baseline ARAT 23-47")
mtext(side = 1, text="Last observed outcome (Days since stroke)", line = 5.5, cex = 1.3)
mtext(side = 2, text="Absolute error", line = 3, cex = 1.3)

Plot_lasttime_perARATbase("model_SA_FE_withInter", "[0-22]", title_main = "Baseline ARAT 0-22")
mtext(side = 1, text="", line = 2.2)
mtext(side = 2, text="Absolute error", line = 3, cex = 1.3)

dev.off()


###########################################################################
###########################################################################

##################
# Figure 1 Suppl #
##################


model_nsCHECK1 <- lme(ARAT ~ ns(Days, df = 2), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, df = 2))), na.action = na.exclude)
model_nsCHECK2 <- lme(ARAT ~ ns(Days, df = 3), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, df = 3))), na.action = na.exclude)
model_nsCHECK3 <- lme(ARAT ~ ns(Days, df = 4), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, df = 4))), na.action = na.exclude)
model_nsCHECK4 <- lme(ARAT ~ ns(Days, knots = c(6, 11, 19, 34, 50, 91)), data = data, 
                      random = list(Number = pdDiag(form = ~ ns(Days, knots = c(6, 11, 19, 34, 50, 91)))), na.action = na.exclude)

set.seed(123)
set.seed(2019)
ind <- sample(unique(data$Number), 49, replace = F)
ind <- sample(unique(data$Number), 12, replace = F)

data$fitted_marg <- fitted(model_nsCHECK1, level = 0)
data$fitted_subj <- fitted(model_nsCHECK1, level = 1)
p1 <- xyplot(ARAT  + fitted_subj ~ Days | Number, data = data, strip = FALSE,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 2)
               y.mat <- matrix(y, ncol = 2)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black", pch = 16)
               #panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "blue")
             }, subset = Number %in% ind, layout = c(4,3), as.table = TRUE,
             xlab = list("time since stroke (days)", cex = 1.2), ylab = list("ARAT", cex = 1.2),scales=list(x=list(cex=1.1), y=list(cex=1.1)), 
             key=list(columns=4, 
                      text=list(lab=c("spline 1 knot", "spline 2 knots", "spline 3 knots", "spline 6 knots")),
                      col=c("blue", "red", "green", "black"), 
                      lines=list(lty=c(1), lwd=2), cex = 1.2))

data$fitted_marg <- fitted(model_nsCHECK2, level = 0)
data$fitted_subj <- fitted(model_nsCHECK2, level = 1)
p2 <- xyplot(ARAT  + fitted_subj ~ Days | Number, data = data, strip = FALSE,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 2)
               y.mat <- matrix(y, ncol = 2)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black", pch = 16)
               #panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "red")
             }, subset = Number %in% ind, layout = c(4, 3), as.table = TRUE,
             xlab = list("time since stroke (days)", cex = 1.2), ylab = list("ARAT", cex = 1.2),scales=list(x=list(cex=1.1), y=list(cex=1.1)),
             key=list(columns=4, 
                      text=list(lab=c("spline 1 knot", "spline 2 knots", "spline 3 knots", "spline 6 knots")),
                      col=c("blue", "red", "green", "black"), 
                      lines=list(lty=c(1), lwd=2), cex = 1.2))


data$fitted_marg <- fitted(model_nsCHECK3, level = 0)
data$fitted_subj <- fitted(model_nsCHECK3, level = 1)
p3 <- xyplot(ARAT + fitted_subj ~ Days | Number, data = data, strip = FALSE,
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 2)
               y.mat <- matrix(y, ncol = 2)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black", pch = 16)
               #panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "green")
             }, subset = Number %in% ind, layout = c(4, 3), as.table = TRUE,
             xlab = list("time since stroke (days)", cex = 1.2), ylab = list("ARAT", cex = 1.2),scales=list(x=list(cex=1.1), y=list(cex=1.1)),
             key=list(columns=4, 
                      text=list(lab=c("spline 1 knot", "spline 2 knots", "spline 3 knots", "spline 6 knots")),
                      col=c("blue", "red", "green", "black"), 
                      lines=list(lty=c(1), lwd=2), cex = 1.2))


data$fitted_marg <- fitted(model_nsCHECK4, level = 0)
data$fitted_subj <- fitted(model_nsCHECK4, level = 1)
p4 <- xyplot(ARAT + fitted_subj ~ Days | Number, data = data, strip = FALSE, 
             panel = function (x, y, ...) {
               x.mat <- matrix(x, ncol = 2)
               y.mat <- matrix(y, ncol = 2)
               panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black", pch = 16)
               #panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "light grey")
               panel.xyplot(x.mat[, 2], y.mat[, 2], type = "l", lwd = 2, col = "black")
             }, subset = Number %in% ind, layout = c(4, 3), as.table = TRUE,
             xlab = list("time since stroke (days)", cex = 1.2), ylab = list("ARAT", cex = 1.2),scales=list(x=list(cex=1.1), y=list(cex=1.1)),
             key=list(columns=4, 
                      text=list(lab=c("spline 1 knot", "spline 2 knots", "spline 3 knots", "spline 6 knots")),
                      col=c("blue", "red", "green", "black"), 
                      lines=list(lty=c(1), lwd=2), cex = 1.2))


p1 + as.layer(p2) + as.layer(p3) + as.layer(p4)



##################
# Figure 2 Suppl #
##################



Plot_models <- function(input.meas){ 
  CV_res_final_used <- CV_res_final[CV_res_final$Measurements_used >= 1 & CV_res_final$Pred_time >= 160 &
                                      CV_res_final$Measurements_used == input.meas, ]
  CV_res_final_used <- CV_res_final_used[CV_res_final_used$Method %in% 
                                           c("model_allMain_allInter", "model_allMain_sigAllInter",
                                           "model_sigAllMain_sigAllInter", 
                                           "model_SA_FE_withInter",
                                           "No Cov"),]
  
  CV_res_final_used$Method <- droplevels(CV_res_final_used$Method)
  
  CV_res_final_used$Method <- ordered(CV_res_final_used$Method,  c("model_allMain_allInter",
                                                                   "model_allMain_sigAllInter", 
                                       "model_sigAllMain_sigAllInter", "model_SA_FE_withInter",
                                       "No Cov"))
 
  levels(CV_res_final_used$Method) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
  
  par(mar = c(5, 5, 2, 1))
  boxplot(CV_res_final_used$AE ~ CV_res_final_used$Method, horizontal=FALSE, xlab = "", 
          cex.axis = 1.2,
          cex.lab = 1.4, ylab = "Measurements used", main = paste0(input.meas, " measurements used"), 
          cex.main = 1.4, ylim = c(0, 60),
          outline = FALSE, las = 2)
}

par(mfrow = c(2,3), oma = c(8,2,3,1), mgp = c(4,1,0), mar = c(7,3,2,1))

Plot_models(1)
Plot_models(2)
Plot_models(3)
Plot_models(4)
Plot_models(5)
Plot_models(6)





##################
# Figure 3 Suppl #
################## 

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
  ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
  data[c(ind), ]
}



DynPlots <- function(model.output = model.output, newdata, timeVar, 
                     main_title = ""){
  
  # Load individual prediction ------------------------------------
  data <- model.output$data
  formYx <- formula(model.output)
  yOutcome <- formYx[[2]]
  
  IndvPrediction95 <- IndvPred_lme(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                   interval = "prediction", return_data = TRUE)
  
  IndvPrediction68 <- IndvPred_lme(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                   interval = "prediction", return_data = TRUE, level = 0.68)
  
  pred95 <- IndvPrediction95[which(!is.na(IndvPrediction95$low)),]
  pred68 <- IndvPrediction68[which(!is.na(IndvPrediction68$low)),]
  
  nopred <- IndvPrediction95[which(is.na(IndvPrediction95$low)),]
  
  timeVariable <- pred95[[timeVar]]
  
  # Generating plot -----------------------------------------------------
  xyplot(pred ~ timeVariable , main = main_title, data = pred95,
         type = "l", col = rgb(0.6769,0.4447,0.7114, alpha = 1), lty = c(1, 2, 2), lwd = 3,
         ylim = c(0,65), xlim = c(0,200), ylab = list(yOutcome, cex = 1.2), xlab = list("Days since stroke", cex = 1.2),
         scales = list(x = list(cex = 1.3) , y = list(cex = 1.3)),
         panel = function(x, y,  ...) {
           panel.xyplot(x, y, ...)
           panel.polygon(c(pred95[,"Days"], rev(pred95[,"Days"])), 
                         c(pred95[,"upp"], rev(pred95[,"low"])),
                         border = NA,
                         col = rgb(0.6769,0.4447,0.7114, alpha = 0.2))
           panel.polygon(c(pred68[,"Days"], rev(pred68[,"Days"])), 
                         c(pred68[,"upp"], rev(pred68[,"low"])),
                         border = NA,
                         col =rgb(0.6769,0.4447,0.7114, alpha = 0.4))
           panel.points(x = nopred[[timeVar]], y = nopred[[yOutcome]], cex = 1.2, pch = 16, col = "black");
           panel.lines(x = rep(tail(nopred[[timeVar]], n = 1), 200), y = seq(-100, 100, length = 200), col = "grey", lty = 3, lwd = 2)
         })
  
  
}


newPatient <- data[data$Number == 134, ]
newPatient <- data[data$Number == 3, ]

newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:3, ]

p11 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                   timeVar = "Days",
                   main_title = "")#,
                   #nameX = "", nameY = "ARAT")
p12 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                   timeVar = "Days")#,
                   #nameX = "", main_title = "")
#grid.arrange(p11, p12, nrow = 2)

##################
newPatient <- data[data$Number == 266, ]

newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:4, ]

p21 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                   timeVar = "Days",
                   main_title = "")#,
                   #nameX = "", nameY = "ARAT")
p22 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                   timeVar = "Days")#,
                   #nameX = "", main_title = "")
#grid.arrange(p21, p22, nrow = 2)

##################
newPatient <- data[data$Number == 25, ]
newPatient <- data[data$Number == 51, ]

newPatient1 <- newPatient[1:4, ]
newPatient2 <- newPatient[1:5, ]

p31 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                   timeVar = "Days",
                   main_title = "")#,
                   #nameX = "", nameY = "ARAT")
p32 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                   timeVar = "Days")#,
                   #nameX = "", main_title = "")

#grid.arrange(p31, p32, nrow = 2)

##################
newPatient <- data[data$Number == 351, ]

newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:4, ]

p41 <- DynPlots(model.output = model_SA_FE, newdata = newPatient1, 
                   timeVar = "Days",
                   main_title = "")#,
                   #nameX = "time since stroke (days)", nameY = "ARAT")
p42 <- DynPlots(model.output = model_SA_FE, newdata = newPatient2, 
                   timeVar = "Days",
                   main_title = "")#,
                   #nameX = "time since stroke (days)", nameY = "")

#grid.arrange(p41, p42, nrow = 2)
### 700 X 787
grid.arrange(p11, p12, p21, p22, p31, p32, p41, p42, nrow = 4, widths=c(1,1))
grid.text("Patient 1", x = unit(0.5, "npc"), 
          y = unit(0.99, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("Patient 2", x = unit(0.5, "npc"), 
          y = unit(0.74, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("Patient 3", x = unit(0.5, "npc"), 
          y = unit(0.49, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("Patient 4", x = unit(0.5, "npc"), 
          y = unit(0.24, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))


grid.text("Time point 1", x = unit(0.28, "npc"), 
          y = unit(0.98, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("Time point 2", x = unit(0.8, "npc"), 
          y = unit(0.98, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
