############################################################
# Relative importance at baseline (h1) and trajectory (h2) #
############################################################
# Sys.setenv(HOME="C:/Users/Wei Jia")
# set.seed(02032024)
set.seed(02025024)
library(cmdstanr)
# library(posterior)
# library(bayesplot)
library(ggplot2)
# color_scheme_set("brightblue")

setwd('C:/Users/Wei Jia/OneDrive - University of Cincinnati/Mixture')
source("BVCMIM/BVCMIM_Sim_Pred.R")
simulated_data = readRDS('BVCMIM/chem/data/chem_data.RDS')
fit = readRDS('BVCMIM/chem/results/chem_menmdi_012.RDS')
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
hfun_3 <- function(z) 4/(1+exp(-12*z)) #1.5/(1+exp(-15*z)) -0.25*z

group_1_ind = c(1:8)
group_2_ind = c(9:24)
group_3_ind = c(25,26)
group_4_ind = c(27:32)
group_5_ind = c(33:36)
ind_ls = list(group_1_ind, group_2_ind, group_3_ind, group_4_ind, group_5_ind)

# hfun_3 = function(z) z
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
h_1 = hfun_1(index_1_1) + hfun_2(index_1_3) + hfun_3(index_1_5) + 0.5*(hfun_1(index_1_1))*(hfun_3(index_1_5))
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
h_2 = hfun_1(index_2_1) + hfun_2(index_2_2) + hfun_3(index_2_4) + 0.5*(hfun_1(index_2_1))*(hfun_3(index_2_4))
index_2 = cbind(index_2_1, index_2_2, index_2_3, index_2_4, index_2_5)

#   xlab("expos1")
qs = c(0.25, 0.75, 0.50)

index_1_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_1_z = index_1_z[rep(seq_len(nrow(index_1_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_1_ind){
  index_1_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_1_z_75 = index_1_z
for (idx in group_5_ind){
  index_1_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_1_z_25 = index_1_z
for (idx in group_5_ind){
  index_1_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_1_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_1_z_25, fit = fit)$postmean[101:200]
weights_1_1_mean = as.matrix(apply(fit$draws('weights_1_1'), 3, mean))
weights_1_5_mean = as.matrix(apply(fit$draws('weights_1_5'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_1_z[,group_1_ind] %*% weights_1_1_mean, index_1_z[,group_1_ind] %*% weights_1_1_mean), q = c(rep(75, 100), rep(25, 100)))
fig_1_5 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F) + 
  # facet_grid(variable2 ~ variable1) +
  # ggtitle("h(index 1 | quantiles of index 5)") +
  xlab("Phth & BPA") +
  ylab("h1") +
  labs(color = "Quantile of PFAS")
fig_1_5
