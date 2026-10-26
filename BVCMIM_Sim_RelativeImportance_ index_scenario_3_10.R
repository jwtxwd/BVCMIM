############################################################
# Relative importance at baseline (h1) and trajectory (h2) #
############################################################
Sys.setenv(HOME="C:/Users/Wei Jia")
# set.seed(02032024)
set.seed(02025024)
library(cmdstanr)
# library(posterior)
# library(bayesplot)
# color_scheme_set("brightblue")

setwd('C:/Users/Wei Jia/OneDrive - University of Cincinnati/Mixture/')
source("BVCMIM/BVCMIM_Sim_Pred.R")
mat_array = array(0, dim = c(10, 2, 5))
mat_sd_array = array(0, dim = c(10, 2, 5))
mat_truth_array = array(0, dim = c(10, 2, 5))
h1_result = array(0, dim = c(10, 4))
h2_result = array(0, dim = c(10, 4))
for (i in 1:10){
  # simulated_data = readRDS('BVCMIM/Simulated_data/scenario_3_5_rho_pos_neg.RDS')
  # simulated_data = readRDS(paste('BVCMIM/simulations/scenario_3/Simulated_data/5_rho_pos_neg/scenario_3_5_rho_pos_neg_', i, '.RDS', sep = ''))
  simulated_data = readRDS(paste('BVCMIM/simulations/scenario_3/Simulated_data/5_rho/scenario_3_5_rho_', i, '.RDS', sep = ''))
  # fit = readRDS('BVCMIM/results/scenario_3_5_rho_pos_neg.RDS')
  # fit = readRDS(paste('BVCMIM/simulations/scenario_3/results/5_rho_pos_neg/scenario_3_5_rho_pos_neg_', i, '.RDS', sep=''))
  fit = readRDS(paste('BVCMIM/simulations/scenario_3/results/5_rho/scenario_3_5_rho_', i, '.RDS', sep=''))
  q = simulated_data$q
  M = simulated_data$M 
  Z = simulated_data$Z
  Y = simulated_data$Y
  X = simulated_data$X
  U = simulated_data$U
  W = simulated_data$W
  n = simulated_data$n
  par(mfrow=c(1,2))
  
  mat 		= matrix(NA, nrow=q, ncol=5)
  mat_sd 		= matrix(NA, nrow=q, ncol=5)
  true_mat	= matrix(rep(0,q*5), nrow=q, ncol=5)
  mat_truth = matrix(NA, nrow=q, ncol=5)
  
  # hfun_1 <- function(z) z
  # hfun_2 <- function(z) z
  # hfun_3 <- function(z) 1.5/(1+exp(-15*z))
  hfun_1 <- function(z) z #dnorm(z-0.2,0,0.5) #(0.85*z)^2-0.5 # h1 <- 3/(1+exp(-2*wx1))
  hfun_2 <- function(z) z
  hfun_3 <- function(z) 1.5/(1+exp(-15*z)) #1.5/(1+exp(-15*z)) -0.25*z
  
  weights_1_1_mean = as.matrix(apply(fit$draws('weights_1_1'), 3, mean))
  weights_1_2_mean = as.matrix(apply(fit$draws('weights_1_2'), 3, mean))
  weights_1_3_mean = as.matrix(apply(fit$draws('weights_1_3'), 3, mean))
  weights_1_4_mean = as.matrix(apply(fit$draws('weights_1_4'), 3, mean))
  weights_1_5_mean = as.matrix(apply(fit$draws('weights_1_5'), 3, mean))
  weights_2_1_mean = as.matrix(apply(fit$draws('weights_2_1'), 3, mean))
  weights_2_2_mean = as.matrix(apply(fit$draws('weights_2_2'), 3, mean))
  weights_2_3_mean = as.matrix(apply(fit$draws('weights_2_3'), 3, mean))
  weights_2_4_mean = as.matrix(apply(fit$draws('weights_2_4'), 3, mean))
  weights_2_5_mean = as.matrix(apply(fit$draws('weights_2_5'), 3, mean))
  index_1_1 = Z[,1:8] %*% weights_1_1_mean
  index_1_1 = scale(index_1_1)
  index_1_2 = Z[,9:24] %*% weights_1_2_mean
  index_1_2 = scale(index_1_2)
  index_1_3 = Z[,25:26] %*% weights_1_3_mean
  index_1_3 = scale(index_1_3)
  index_1_4 = Z[,27:32] %*% weights_1_4_mean
  index_1_4 = scale(index_1_4)
  index_1_5 = Z[,33:36] %*% weights_1_5_mean
  index_1_5 = scale(index_1_5)
  h_1 = hfun_1(index_1_1) + hfun_2(index_1_3) + hfun_3(index_1_5) + 1*(hfun_1(index_1_1))*(hfun_3(index_1_5))
  index_1 = cbind(index_1_1, index_1_2, index_1_3, index_1_4, index_1_5)
  # weights_2_1 = c(8, 7, 6, 5, 4, 3 , 2, 1)
  # weights_2_1 <- weights_2_1/sqrt(sum(weights_2_1^2))
  # weights_2_2 = c(8, 7, 6, 5, 4, 3, 2, 1, 8, 7, 6, 5, 4, 3 , 2, 1)
  # weights_2_2 <- weights_2_2/sqrt(sum(weights_2_2^2))
  # weights_2_3 = c(0,0)
  # weights_2_4 = c(6, 5, 4, 3, 2, 1)
  # weights_2_4 <- weights_2_4/sqrt(sum(weights_2_4^2))
  # weights_2_5 = c(0,0,0,0)
  index_2_1 = Z[,1:8] %*% weights_2_1_mean
  index_2_1 = scale(index_2_1)
  index_2_2 = Z[,9:24] %*% weights_2_2_mean
  index_2_2 = scale(index_2_2)
  index_2_3 = Z[,25:26] %*% weights_2_3_mean
  index_2_3 = scale(index_2_3)
  index_2_4 = Z[,27:32] %*% weights_2_4_mean
  index_2_4 = scale(index_2_4)
  index_2_5 = Z[,33:36] %*% weights_2_5_mean
  index_2_5 = scale(index_2_5)
  h_2 = hfun_1(index_2_1) + hfun_2(index_2_2) + hfun_3(index_2_4) + 1*(hfun_1(index_2_1))*(hfun_3(index_2_4))
  index_2 = cbind(index_2_1, index_2_2, index_2_3, index_2_4, index_2_5)
  
  weights_1_1 = c(0.5, 0.6, 0.5, 0, 0, 0 , 0, 0)
  # weights_1_1 <- weights_1_1/sqrt(sum(weights_1_1^2))
  weights_1_2 = c(0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0)
  # weights_1_2 <- weights_1_2/sqrt(sum(weights_1_2^2))
  weights_1_3 = c(0,0)
  # weights_1_3 <- weights_1_3/sqrt(sum(weights_1_3^2))
  weights_1_4 = c(0, 0, 0, 0, 0, 0)
  # weights_1_4 <- weights_1_4/sqrt(sum(weights_1_4^2))
  weights_1_5 = c(0, 0, 0.5, 0.5)
  # weights_1_5 <- weights_1_5/sqrt(sum(weights_1_5^2))
  index_1_1_r = Z[,1:8] %*% weights_1_1
  index_1_1_r = scale(index_1_1_r)
  index_1_2_r = Z[,9:24] %*% weights_1_2
  # index_1_2_r = scale(index_1_2_r)
  index_1_3_r = Z[,25:26] %*% weights_1_3
  # index_1_3_r = scale(index_1_3_r)
  index_1_4_r = Z[,27:32] %*% weights_1_4
  # index_1_4_r = scale(index_1_4_r)
  index_1_5_r = Z[,33:36] %*% weights_1_5
  index_1_5_r = scale(index_1_5_r)
  index_1_r = cbind(index_1_1_r, rep(0, length(index_1_1_r)), index_1_3_r, rep(0, length(index_1_1_r)), index_1_5_r)
  h_1_r = hfun_1(index_1_1_r) + hfun_2(index_1_3_r) + hfun_3(index_1_5_r) + 2*(hfun_1(index_1_1_r))*(hfun_3(index_1_5_r))
  
  weights_2_1 = c(0, 0, 0, 0, 0, 0 , 0, 0)
  # weights_2_1 <- weights_2_1/sqrt(sum(weights_2_1^2))
  weights_2_2 = c(0.6, 0.4, 0.5, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0)
  # weights_2_2 <- weights_2_2/sqrt(sum(weights_2_2^2))
  weights_2_3 = c(0,0)
  # weights_2_3 <- weights_2_3/sqrt(sum(weights_2_3^2))
  weights_2_4 = c(0, 0, 0, 0, 0.4, 0.2)
  # weights_2_4 <- weights_2_4/sqrt(sum(weights_2_4^2))
  weights_2_5 = c(0, 0, 0, 0)
  # weights_2_5 <- weights_2_5/sqrt(sum(weights_2_5^2))
  index_2_1_r = Z[,1:8] %*% weights_2_1
  # index_2_1_r = scale(index_2_1_r)
  index_2_2_r = Z[,9:24] %*% weights_2_2
  index_2_2_r = scale(index_2_2_r)
  index_2_3_r = Z[,25:26] %*% weights_2_3
  # index_2_3_r = scale(index_2_3_r)
  index_2_4_r = Z[,27:32] %*% weights_2_4
  index_2_4_r = scale(index_2_4_r)
  index_2_5_r = Z[,33:36] %*% weights_2_5
  # index_2_5_r = scale(index_2_5_r)
  index_2_r = cbind(index_2_1_r, index_2_2_r, rep(0, length(index_1_1_r)), index_2_4_r, rep(0, length(index_1_1_r)))
  h_2_r = hfun_1(index_2_2_r) + hfun_2(index_2_3_r) + hfun_3(index_2_4_r) + 2*(hfun_1(index_2_2_r))*(hfun_3(index_2_4_r))
  
  true_h1 = h_1_r
  true_h2 = h_2_r
  
  hhat_1	= apply(fit$draws('h1'), 3, mean)
  hhat_2	= apply(fit$draws('h2'), 3, mean)
  
  model_1 = lm(hhat_1~true_h1)
  model_2 = lm(hhat_2~true_h2)
   
  h1_result[i,] = c(summary(model_1)$coef[1,1], summary(model_1)$coef[2,1], summary(model_1)$r.sq, sqrt( sum( (hhat_1 - true_h1)^2) / n))
  h2_result[i,] = c(summary(model_2)$coef[1,1], summary(model_2)$coef[2,1], summary(model_2)$r.sq, sqrt( sum( (hhat_2 - true_h2)^2) / n))
  
  qs = c(0.25, 0.75, 0.50)
  # compute_truth = function(z){
  #   true_h1 = hfun_1(z[1,1]) + hfun_3(z[1,5]) + 2*hfun_1(z[1,1])*hfun_3(z[1,5]) - (hfun_1(z[2,1]) + hfun_3(z[2,5]) + 2*hfun_1(z[2,1])*hfun_3(z[2,5]))
  #   true_h2 = hfun_1(z[3,2]) + hfun_3(z[3,4]) + 2*hfun_1(z[3,2])*hfun_3(z[3,4]) - (hfun_1(z[4,2]) + hfun_3(z[4,4]) + 2*hfun_1(z[4,2])*hfun_3(z[4,4]))
  #   list(true_h1, true_h2)
  # }
  compute_truth = function(Z){
    true_h1 = (hfun_1(Z[1,1:8] %*% weights_1_1) + hfun_3(Z[1,33:36] %*% weights_1_5) + 2*(hfun_1(Z[1,1:8] %*% weights_1_1))*hfun_3(Z[1,33:36] %*% weights_1_5))- (hfun_1(Z[2,1:8] %*% weights_1_1) + hfun_3(Z[2,33:36] %*% weights_1_5) + 2*(hfun_1(Z[2,1:8] %*% weights_1_1))*hfun_3(Z[2,33:36] %*% weights_1_5))
    true_h2 = (hfun_1(Z[1,9:24] %*% weights_2_2) + hfun_3(Z[1,27:32] %*% weights_2_4) + 2*(hfun_1(Z[1,9:24] %*% weights_2_2))*hfun_3(Z[1,27:32] %*% weights_2_4))- (hfun_1(Z[2,9:24] %*% weights_2_2) + hfun_3(Z[2,27:32] %*% weights_2_4) + 2*(hfun_1(Z[2,9:24] %*% weights_2_2))*hfun_3(Z[2,27:32] %*% weights_2_4))
    list(true_h1, true_h2)
  }



  group_1_ind = c(1:8)
  group_2_ind = c(9:24)
  group_3_ind = c(25,26)
  group_4_ind = c(27:32)
  group_5_ind = c(33:36)
  ind_ls = list(group_1_ind, group_2_ind, group_3_ind, group_4_ind, group_5_ind)
  for (j in 1:5) {
    cross_sec 		= rbind(apply(Z, 2, median), apply(Z, 2, median))
    ind = ind_ls[[j]]
    for (idx in ind){
      cross_sec[,idx] 	= c(quantile(Z[,idx], qs[2]), quantile(Z[,idx], qs[1]))
    }
    cross_sec_true 		= rbind(apply(index_1_r, 2, median), apply(index_1_r, 2, median), apply(index_2_r, 2, median), apply(index_2_r, 2, median))
    cross_sec_true[,j] 	= c(quantile(index_1_r[,j], qs[2]), quantile(index_1_r[,j], qs[1]), quantile(index_2_r[,j], qs[2]), quantile(index_2_r[,j], qs[1]))
    hgrid_cross_sec = newh.postmean(Z = Z, Znew = cross_sec, fit = fit) #returns 75% h2, 25% h2, 75% h1, 25% h1
    mat[1,j] 	= hgrid_cross_sec$postmean[3] - hgrid_cross_sec$postmean[4]
    mat[2,j] 	= hgrid_cross_sec$postmean[1] - hgrid_cross_sec$postmean[2]
    mat_sd[1,j] 	= sqrt(hgrid_cross_sec$postvar[3,3] + hgrid_cross_sec$postvar[4,4] - 2*hgrid_cross_sec$postvar[3,4])
    mat_sd[2,j] 	= sqrt(hgrid_cross_sec$postvar[1,1] + hgrid_cross_sec$postvar[2,2] - 2*hgrid_cross_sec$postvar[1,2])
    cross_sec_truth = compute_truth(cross_sec)
    mat_truth[1,j] = cross_sec_truth[[1]]
    mat_truth[2,j] = cross_sec_truth[[2]]
  }
  
  # for (j in 1:5) {
  #   cross_sec_1 		= rbind(apply(index_1, 2, median), apply(index_1, 2, median))
  #   cross_sec_2     = rbind(apply(index_2, 2, median), apply(index_2, 2, median))
  #   
  #   cross_sec_1[,j] 	= c(quantile(index_1[,j], qs[2]), quantile(index_1[,j], qs[1]))
  #   cross_sec_2[,j] 	= c(quantile(index_2[,j], qs[2]), quantile(index_2[,j], qs[1]))
  #   
  #   cross_sec_true 		= rbind(apply(index_1_r, 2, median), apply(index_1_r, 2, median), apply(index_2_r, 2, median), apply(index_2_r, 2, median))
  #   cross_sec_true[,j] 	= c(quantile(index_1_r[,j], qs[2]), quantile(index_1_r[,j], qs[1]), quantile(index_2_r[,j], qs[2]), quantile(index_2_r[,j], qs[1]))
  #   # hgrid_cross_sec = newh.postmean(Z = Z, Znew = cross_sec, fit = fit) #returns 75% h2, 25% h2, 75% h1, 25% h1
  #   hgrid_cross_sec = newh_postmean_index(index_1, index_2, cross_sec_1, cross_sec_2, fit)
  #   mat[1,j] 	= hgrid_cross_sec$postmean[3] - hgrid_cross_sec$postmean[4]
  #   mat[2,j] 	= hgrid_cross_sec$postmean[1] - hgrid_cross_sec$postmean[2]
  #   mat_sd[1,j] 	= sqrt(hgrid_cross_sec$postvar[3,3] + hgrid_cross_sec$postvar[4,4] - 2*hgrid_cross_sec$postvar[3,4])
  #   mat_sd[2,j] 	= sqrt(hgrid_cross_sec$postvar[1,1] + hgrid_cross_sec$postvar[2,2] - 2*hgrid_cross_sec$postvar[1,2])
  #   cross_sec_truth = compute_truth(rbind(cross_sec_1, cross_sec_2))
  #   mat_truth[1,j] = cross_sec_truth[[1]]
  #   mat_truth[2,j] = cross_sec_truth[[2]]
  # }
  mat_array[i,,] = mat
  mat_sd_array[i,,] = mat_sd
  mat_truth_array[i,,] = mat_truth
}

mat = apply(mat_array, c(2,3), mean)
mat_sd = apply(mat_sd_array, c(2,3), mean)
mat_truth = apply(mat_truth_array, c(2,3), mean)

Summary_h 	= matrix(NA, nrow=q, ncol=4)
Summary_h[1,] = apply(h1_result, 2, mean)
Summary_h[2,] = apply(h2_result, 2, mean)
Summary_h

par(mfrow=c(1,2))
#Plot
ylim 		= c(-1, 7)
plot(1:5, mat[1,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Exposure index", ylab="Main effect of each exposure index",
    main="Baseline (h1)"
)
points(1:5, mat_truth[1,], col='red', pch=16)
# hack: we draw arrows but with very special "arrowheads"
arrows(1:5, mat[1,] - 1.96*mat_sd[1,], 1:5, mat[1,] + 1.96*mat_sd[1,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c(expression(z^{'*'}~1), expression(z^{'*'}~2), expression(z^{'*'}~3), expression(z^{'*'}~4), expression(z^{'*'}~5)), cex.axis=1.2)
abline(h=0)

#Plot
ylim 		= c(-1, 7)
plot(1:5, mat[2,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Exposure index", ylab="Main effect of each exposure index",
    main="Trajectory (h2)"
)
points(1:5, mat_truth[2,], col='red', pch=16)
# hack: we draw arrows but with very special "arrowheads"
arrows(1:5, mat[2,] - 1.96*mat_sd[2,], 1:5, mat[2,] + 1.96*mat_sd[2,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c(expression(z^{'*'}~1), expression(z^{'*'}~2), expression(z^{'*'}~3), expression(z^{'*'}~4), expression(z^{'*'}~5)), cex.axis=1.2)
abline(h=0)

