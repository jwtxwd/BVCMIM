Sys.setenv(HOME="C:/Users/Wei Jia")
# set.seed(02032024)
set.seed(02025024)
library(cmdstanr)
#library(posterior)
#library(bayesplot)
#color_scheme_set("brightblue")

setwd("C:/Mixture")

mod <- cmdstan_model('stan/BVCMIM_hs.stan')

source("BVCMIM/BVCMIM_hs_Sim_Data_36_scenario_2.R")

simulated_data = list(n= 100,
  q= 2,
  p= 2,
  M= 36,
  T= 5,
  Z= Z,
  Y= drop(Y),
  X= X,
  U= U,
  W= W,
  G = 5)
  # regularized=FALSE,
  # nu=1,
  # global_scale=2,
  # global_df=4,
  # slab_scale=2,
  # slab_df=4)
# mod = mod$compile(cpp_options = list(O=0))
saveRDS(list(n= 100,
  q= 2,
  p= 2,
  M= 36,
  T= 5,
  Z = Z,
  Z_quant= Z,
  Y= drop(Y),
  X= X,
  U= U,
  W= W,
  regularized=FALSE,
  nu=1,
  global_scale=2,
  global_df=4,
  slab_scale=2,
  slab_df=4), 'BVCMIM/Simulated_data/scenario_2.RDS')
# # # 
fit <- mod$sample(
  data = simulated_data,
  seed = 111,
  chains = 1,
  parallel_chains = 1,
  refresh = 50, # print update every 500 iters
  iter_warmup = 1000,
  iter_sampling = 1000
)
fit$save_object('BVCMIM/results/scenario_2.RDS')
# fit = mod$variational(data = simulated_data, seed=123, draws = 4000, grad_samples = 1)
# fit = mod$pathfinder(data = simulated_data, seed=123, draws = 4000)
# fit = readRDS('bvckmr_mcmc.RDS')
Summary_h 	= matrix(NA, nrow=q, ncol=4)

# #scenario 1
# true_h1 = 5*(0.1*Z[, 1] + 0.4*Z[,12] + 0.2*Z[,35] + 0.5*Z[,12]*Z[,35])
# true_h2 = 3*(0.1*Z[, 9] + 0.5*Z[,18]^2 + 0.2*Z[,33])

#scenario 2
true_h1 = 5*(0.1*Z[, 1] + 0.4*Z[,12] + 0.2*Z[,35] + 0.5*Z[,12]*Z[,35])
true_h2 = 3*(0.1*Z[, 2] + 0.5*Z[,18] + 0.2*Z[,33]^2)

# #scenario 3
# true_h1 = 5*(0.1*Z[, 1] + 0.3*Z[, 12] - 0*Z[, 23] + 0.4*Z[,24] + 0.2*Z[,35])
# true_h2 = 3*(0*Z[, 6] + 0.1*Z[, 9] + 0.5*Z[, 23] + 0.2*Z[,27] + 0.2*Z[,33])

# #scenario 4
# true_h1 = 5*(0.1*Z[, 1] + 0.3*Z[, 12] - 0.3*Z[, 23] + 0.4*Z[,24] + 0.2*Z[,35])
# true_h2 = 3*(-0.2*Z[, 2] + 0.1*Z[, 9] + 0.5*Z[, 23] + 0.2*Z[,27] + 0.2*Z[,33])
# hhat_1	= apply(fit$draws('h')[,1,], 3, mean)[1:n]
# hhat_2	= apply(fit$draws('h')[,1,], 3, mean)[as.integer(n+1):as.integer(2*n)]

hhat_1	= apply(fit$draws('h1'), 3, mean)
hhat_2	= apply(fit$draws('h2'), 3, mean)

model_1 = lm(hhat_1~true_h1)
model_2 = lm(hhat_2~true_h2)

Summary_h[1,] = c(summary(model_1)$coef[1,1], summary(model_1)$coef[2,1], summary(model_1)$r.sq, sqrt( sum( (hhat_1 - true_h1)^2) / n))
Summary_h[2,] = c(summary(model_2)$coef[1,1], summary(model_2)$coef[2,1], summary(model_2)$r.sq, sqrt( sum( (hhat_2 - true_h2)^2) / n))
Summary_h