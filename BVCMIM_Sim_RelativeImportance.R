############################################################
# Relative importance at baseline (h1) and trajectory (h2) #
############################################################
Sys.setenv(HOME="C:/Users/Wei Jia")
# set.seed(02032024)
set.seed(02025024)
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

setwd("C:/Mixture")
source("BVCMIM/BVCMIM_Sim_Pred.R")
simulated_data = readRDS('BVCMIM/Simulated_data/scenario_3_1_rho.RDS')
fit = readRDS('BVCMIM/results/scenario_3_1_rho.RDS')
q = simulated_data$q
M = simulated_data$M 
Z = simulated_data$Z
Y = simulated_data$Y
X = simulated_data$X
U = simulated_data$U
W = simulated_data$W
n = simulated_data$n
par(mfrow=c(1,2))

mat 		= matrix(NA, nrow=q, ncol=M)
mat_sd 		= matrix(NA, nrow=q, ncol=M)
true_mat	= matrix(rep(0,q*M), nrow=q, ncol=M)
mat_truth = matrix(NA, nrow=q, ncol=M)

hfun_1 <- function(z) z
hfun_2 <- function(z) z
hfun_3 <- function(z) z
#scenario 1 single index
# weights_1_1 = c(8, 7, 6, 5, 4, 3 , 2, 1)
# weights_1_1 <- weights_1_1/sqrt(sum(weights_1_1^2))
# weights_1_3 = c(2,3)
# weights_1_3 <- weights_1_3/sqrt(sum(weights_1_3^2))
# weights_1_5 = c(4, 3, 2, 1)
# weights_1_5 <- weights_1_5/sqrt(sum(weights_1_5^2))
# index_1_1 = Z[,1:8] %*% weights_1_1
# index_1_1 = scale(index_1_1)
# index_1_3 = Z[,25:26] %*% weights_1_3
# index_1_3 = scale(index_1_3)
# index_1_5 = Z[,33:36] %*% weights_1_5
# index_1_5 = scale(index_1_5)
# h_1 = hfun_1(index_1_1) + hfun_2(index_1_3) + hfun_3(index_1_5) + 0.5*(hfun_1(index_1_1))*(hfun_3(index_1_5))
# 
# weights_2_1 = c(8, 7, 6, 5, 4, 3 , 2, 1)
# weights_2_1 <- weights_2_1/sqrt(sum(weights_2_1^2))
# weights_2_2 = c(8, 7, 6, 5, 4, 3, 2, 1, 8, 7, 6, 5, 4, 3 , 2, 1)
# weights_2_2 <- weights_2_2/sqrt(sum(weights_2_2^2))
# weights_2_4 = c(6, 5, 4, 3, 2, 1)
# weights_2_4 <- weights_2_4/sqrt(sum(weights_2_4^2))
# index_2_1 = Z[,1:8] %*% weights_2_1
# index_2_1 = scale(index_2_1)
# index_2_2 = Z[,9:24] %*% weights_2_2
# index_2_2 = scale(index_2_2)
# index_2_4 = Z[,27:32] %*% weights_2_4
# index_2_4 = scale(index_2_4)
# h_2 = hfun_1(index_2_1) + hfun_2(index_2_2) + hfun_3(index_2_4) + 0.5*(hfun_1(index_2_1))*(hfun_3(index_2_4))

# weights_1_1 = c(8, 7, 0, 0, 0, 0 , 0, 0)
# weights_1_1 <- weights_1_1/sqrt(sum(weights_1_1^2))
# weights_1_3 = c(2,0)
# weights_1_3 <- weights_1_3/sqrt(sum(weights_1_3^2))
# weights_1_5 = c(4, 3, 0, 0)
# weights_1_5 <- weights_1_5/sqrt(sum(weights_1_5^2))
# index_1_1 = Z[,1:8] %*% weights_1_1
# index_1_1 = scale(index_1_1)
# index_1_3 = Z[,25:26] %*% weights_1_3
# index_1_3 = scale(index_1_3)
# index_1_5 = Z[,33:36] %*% weights_1_5
# index_1_5 = scale(index_1_5)
# h_1 = hfun_1(index_1_1) + hfun_2(index_1_3) + hfun_3(index_1_5) + 0.5*(hfun_1(index_1_1))*(hfun_3(index_1_5))
# 
# weights_2_1 = c(0, 0, 6, 5, 0, 0 , 0, 0)
# weights_2_1 <- weights_2_1/sqrt(sum(weights_2_1^2))
# weights_2_2 = c(8, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 3 , 2, 5)
# weights_2_2 <- weights_2_2/sqrt(sum(weights_2_2^2))
# weights_2_4 = c(6, 5, 0, 0, 0, 0)
# weights_2_4 <- weights_2_4/sqrt(sum(weights_2_4^2))
# index_2_1 = Z[,1:8] %*% weights_2_1
# index_2_1 = scale(index_2_1)
# index_2_2 = Z[,9:24] %*% weights_2_2
# index_2_2 = scale(index_2_2)
# index_2_4 = Z[,27:32] %*% weights_2_4
# index_2_4 = scale(index_2_4)
# h_2 = hfun_1(index_2_1) + hfun_2(index_2_2) + hfun_3(index_2_4) + 0.5*(hfun_1(index_2_1))*(hfun_3(index_2_4))

# #scenario 3
weights_1_1 = c(0.5, 0, 0, 0, 0, 0 , 0, 0)
# weights_1_1 <- weights_1_1/sqrt(sum(weights_1_1^2))
weights_1_5 = c(0, 0, 1, 0)
# weights_1_5 <- weights_1_5/sqrt(sum(weights_1_5^2))
index_1_1 = Z[,1:8] %*% weights_1_1
index_1_1 = scale(index_1_1)
index_1_5 = Z[,33:36] %*% weights_1_5
index_1_5 = scale(index_1_5)
h_1 = index_1_1 + index_1_5 + 0.5*index_1_1*index_1_5 

weights_2_2 = c(0.6, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0)
# weights_2_1 <- weights_2_1/sqrt(sum(weights_2_1^2))
weights_2_4 = c(0, 0, 0, 0, 0, 1)
# weights_2_4 <- weights_2_4/sqrt(sum(weights_2_4^2))
index_2_2 = Z[,9:24] %*% weights_2_2
index_2_2 = scale(index_2_2)
index_2_4 = Z[,27:32] %*% weights_2_4
index_2_4 = scale(index_2_4)
h_2 = index_2_2 + index_2_4 + 0.5*index_2_2*index_2_4

true_h1 = h_1
true_h2 = h_2

qs = c(0.25, 0.75, 0.50)
# compute_truth = function(Z){
#   true_h1 = (hfun_1(Z[1,1:8] %*% weights_1_1) + hfun_2(Z[1,25:26] %*% weights_1_3) + hfun_3(Z[1,33:36] %*% weights_1_5) + 0.5*(hfun_1(Z[1,1:8] %*% weights_1_1))*hfun_3(Z[1,33:36] %*% weights_1_5))- (hfun_1(Z[2,1:8] %*% weights_1_1) + hfun_2(Z[2,25:26] %*% weights_1_3) + hfun_3(Z[2,33:36] %*% weights_1_5) + 0.5*(hfun_1(Z[2,1:8] %*% weights_1_1))*hfun_3(Z[2,33:36] %*% weights_1_5))
#   true_h2 = (hfun_1(Z[1,1:8] %*% weights_2_1) + hfun_2(Z[1,9:24] %*% weights_2_2) + hfun_3(Z[1,27:32] %*% weights_2_4) + 0.5*(hfun_1(Z[1,1:8] %*% weights_2_1))*hfun_3(Z[1,27:32] %*% weights_2_4))- (hfun_1(Z[2,1:8] %*% weights_2_1) + hfun_2(Z[2,9:24] %*% weights_2_2) + hfun_3(Z[2,27:32] %*% weights_2_4) + 0.5*(hfun_1(Z[2,1:8] %*% weights_2_1))*hfun_3(Z[2,27:32] %*% weights_2_4))
#   list(true_h1, true_h2)
# }
compute_truth = function(Z){
  true_h1 = (hfun_1(Z[1,1:8] %*% weights_1_1) + hfun_3(Z[1,33:36] %*% weights_1_5) + 0.5*(hfun_1(Z[1,1:8] %*% weights_1_1))*hfun_3(Z[1,33:36] %*% weights_1_5))- (hfun_1(Z[2,1:8] %*% weights_1_1) + hfun_3(Z[2,33:36] %*% weights_1_5) + 0.5*(hfun_1(Z[2,1:8] %*% weights_1_1))*hfun_3(Z[2,33:36] %*% weights_1_5))
  true_h2 = (hfun_1(Z[1,9:24] %*% weights_2_2) + hfun_3(Z[1,27:32] %*% weights_2_4) + 0.5*(hfun_1(Z[1,9:24] %*% weights_2_2))*hfun_3(Z[1,27:32] %*% weights_2_4))- (hfun_1(Z[2,9:24] %*% weights_2_2) + hfun_3(Z[2,27:32] %*% weights_2_4) + 0.5*(hfun_1(Z[2,9:24] %*% weights_2_2))*hfun_3(Z[2,27:32] %*% weights_2_4))
  list(true_h1, true_h2)
}


for (j in 1:M) {
  cross_sec 		= rbind(apply(Z, 2, median), apply(Z, 2, median))
  cross_sec[,j] 	= c(quantile(Z[,j], qs[2]), quantile(Z[,j], qs[1]))
  hgrid_cross_sec = newh.postmean(Z = Z, Znew = cross_sec, fit = fit) #returns 75% h2, 25% h2, 75% h1, 25% h1
  mat[1,j] 	= hgrid_cross_sec$postmean[3] - hgrid_cross_sec$postmean[4]
  mat[2,j] 	= hgrid_cross_sec$postmean[1] - hgrid_cross_sec$postmean[2]
  mat_sd[1,j] 	= sqrt(hgrid_cross_sec$postvar[3,3] + hgrid_cross_sec$postvar[4,4] - 2*hgrid_cross_sec$postvar[3,4])
  mat_sd[2,j] 	= sqrt(hgrid_cross_sec$postvar[1,1] + hgrid_cross_sec$postvar[2,2] - 2*hgrid_cross_sec$postvar[1,2])
  cross_sec_truth = compute_truth(cross_sec)
  mat_truth[1,j] = cross_sec_truth[[1]]
  mat_truth[2,j] = cross_sec_truth[[2]]
}

chem_mix = read.csv('preprocessed/data_motpdi.csv')[,c('MBP',	'MIBP',	'MEP2',	'MBZP2'
                                                       ,'MCPP',	'MEHHP',	'MECPP',	'BPA',	'PCB28',	'PCB74'
                                                       ,'PCB99',	'PCB105',	'PCB118',	'PCB138_158',	'PCB146'
                                                       ,'PCB153',	'PCB156',	'PCB170',	'PCB183',	'PCB187'
                                                       ,'PCB194','PCB196_203',	'PCB199',	'PCB206',	'PP_DDE'
                                                       ,'T_NONA',	'BB153',	'PBDE28',	'PBDE47',	'PBDE99',	'PBDE100'
                                                       ,'PBDE153',	'PFOA',	'PFOS',	'PFNA',	'PFHXS')]
chem_names = names(chem_mix)
chem_names <- chem_names[chem_names != "participant_id"]
chem_group = c()
for (name in chem_names){
  if (name %in% c('MBP', 'MIBP', 'MEP2', 'MBZP2', 'MCPP'
                  , 'MEHP', 'MEHHP', 'MECPP', 'BPA')){
    chem_group = c(chem_group, 1)
  }else if (name %in% c('PCB28', 'PCB66', 'PCB74', 'PCB99', 'PCB101', 'PCB105'
                        ,'PCB118','PCB138_158', 'PCB146', 'PCB153', 'PCB156'
                        , 'PCB157', 'PCB167', 'PCB170', 'PCB172', 'PCB177'
                        , 'PCB178', 'PCB183', 'PCB187', 'PCB194', 'PCB195'
                        , 'PCB196_203', 'PCB199', 'PCB206', 'PCB209')){
    chem_group = c(chem_group, 2)
  }else if (name %in% c('B_HCCH', 'PP_DDT', 'PP_DDE', 'OXYCHLOR', 'T_NONA')){
    chem_group = c(chem_group, 3)
  }else if (name %in% c('BB153', 'PBDE28', 'PBDE47', 'PBDE85'
                        , 'PBDE99', 'PBDE100', 'PBDE153', 'PBDE154')){
    chem_group = c(chem_group, 4)
  }else if (name %in% c('PFOA', 'PFOS', 'PFNA', 'PFHXS')){
    chem_group = c(chem_group, 5)
  }
}
group_1_ind = c(1:8)
group_2_ind = c(9:24)
group_3_ind = c(25,26)
group_4_ind = c(27:32)
group_5_ind = c(33:36)

par(mfrow=c(1,2))
#Plot
ylim 		= c(-1, 2)
plot(1:8, mat[1,group_1_ind], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Baseline (h1)"
)
points(1:8, mat_truth[1,group_1_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:8, mat[1,group_1_ind] - 1.96*mat_sd[1,group_1_ind], 1:8, mat[1,group_1_ind] + 1.96*mat_sd[1,group_1_ind], length=0.05, angle=90, code=3)
axis(1, at=1:8, labels=c('Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', 'Z8'))
abline(h=0)

#Plot
ylim 		= c(-1, 2)
plot(1:length(group_1_ind), mat[2,group_1_ind], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Trajectory (h2)"
)
points(1:length(group_1_ind), mat_truth[2,group_1_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_1_ind), mat[2,group_1_ind] - 1.96*mat_sd[2,group_1_ind], 1:length(group_1_ind), mat[2,group_1_ind] + 1.96*mat_sd[2,group_1_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_1_ind), labels=c('Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', 'Z8'))
abline(h=0)

par(mfrow=c(1,2))
#Plot
ylim 		= c(-1, 1)
plot(1:length(group_2_ind), mat[1,group_2_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Baseline (h1)"
)
points(1:length(group_2_ind), mat_truth[1,group_2_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_2_ind), mat[1,group_2_ind] - 1.96*mat_sd[1,group_2_ind], 1:length(group_2_ind), mat[1,group_2_ind] + 1.96*mat_sd[1,group_2_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_2_ind), labels=c('Z9', 'Z10', 'Z11', 'Z12', 'Z13', 'Z14', 'Z15', 'Z16'
                                           , 'Z17', 'Z18', 'Z19', 'Z20', 'Z21', 'Z22', 'Z23', 'Z24'))
abline(h=0)

#Plot
ylim 		= c(-1, 1)
plot(1:length(group_2_ind), mat[2,group_2_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Trajectory (h2)"
)
points(1:length(group_2_ind), mat_truth[2,group_2_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_2_ind), mat[2,group_2_ind] - 1.96*mat_sd[2,group_2_ind], 1:length(group_2_ind), mat[2,group_2_ind] + 1.96*mat_sd[2,group_2_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_2_ind), labels=c('Z9', 'Z10', 'Z11', 'Z12', 'Z13', 'Z14', 'Z15', 'Z16'
                                           , 'Z17', 'Z18', 'Z19', 'Z20', 'Z21', 'Z22', 'Z23', 'Z24'))
abline(h=0)

par(mfrow=c(1,2))
#Plot
ylim 		= c(-2, 4)
plot(1:length(group_3_ind), mat[1,group_3_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Baseline (h1)"
)
points(1:length(group_3_ind), mat_truth[1,group_3_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_3_ind), mat[1,group_3_ind] - 1.96*mat_sd[1,group_3_ind], 1:length(group_3_ind), mat[1,group_3_ind] + 1.96*mat_sd[1,group_3_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_3_ind), labels=c('Z25', 'Z26'))
abline(h=0)

#Plot
ylim 		= c(-3, 4)
plot(1:length(group_3_ind), mat[2,group_3_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Trajectory (h2)"
)
points(1:length(group_3_ind), mat_truth[2,group_3_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_3_ind), mat[2,group_3_ind] - 1.96*mat_sd[2,group_3_ind], 1:length(group_3_ind), mat[2,group_3_ind] + 1.96*mat_sd[2,group_3_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_3_ind), labels=c('Z25', 'Z26'))
abline(h=0)


par(mfrow=c(1,2))
#Plot
ylim 		= c(-2, 4)
plot(1:length(group_4_ind), mat[1,group_4_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Baseline (h1)"
)
points(1:length(group_4_ind), mat_truth[1,group_4_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_4_ind), mat[1,group_4_ind] - 1.96*mat_sd[1,group_4_ind], 1:length(group_4_ind), mat[1,group_4_ind] + 1.96*mat_sd[1,group_4_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_4_ind), labels=c('Z27', 'Z28', 'Z29', 'Z30', 'Z31', 'Z32'))
abline(h=0)

#Plot
ylim 		= c(-3, 4)
plot(1:length(group_4_ind), mat[2,group_4_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Trajectory (h2)"
)
points(1:length(group_4_ind), mat_truth[2,group_4_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_4_ind), mat[2,group_4_ind] - 1.96*mat_sd[2,group_4_ind], 1:length(group_4_ind), mat[2,group_4_ind] + 1.96*mat_sd[2,group_4_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_4_ind), labels=c('Z27', 'Z28', 'Z29', 'Z30', 'Z31', 'Z32'))
abline(h=0)


par(mfrow=c(1,2))
#Plot
ylim 		= c(-2, 4)
plot(1:length(group_5_ind), mat[1,group_5_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Baseline (h1)"
)
points(1:length(group_5_ind), mat_truth[1,group_5_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_5_ind), mat[1,group_5_ind] - 1.96*mat_sd[1,group_5_ind], 1:length(group_5_ind), mat[1,group_5_ind] + 1.96*mat_sd[1,group_5_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_5_ind), labels=c('Z33', 'Z34', 'Z35', 'Z36'))
abline(h=0)

#Plot
ylim 		= c(-3, 4)
plot(1:length(group_5_ind), mat[2,group_5_ind], xaxt="n",
     ylim=ylim,
     pch=15, xlab="Metal", ylab="Main effect of each metal",
     main="Trajectory (h2)"
)
points(1:length(group_5_ind), mat_truth[2,group_5_ind])
# hack: we draw arrows but with very special "arrowheads"
arrows(1:length(group_5_ind), mat[2,group_5_ind] - 1.96*mat_sd[2,group_5_ind], 1:length(group_5_ind), mat[2,group_5_ind] + 1.96*mat_sd[2,group_5_ind], length=0.05, angle=90, code=3)
axis(1, at=1:length(group_5_ind), labels=c('Z33', 'Z34', 'Z35', 'Z36'))
abline(h=0)