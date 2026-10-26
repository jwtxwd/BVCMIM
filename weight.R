Sys.setenv(HOME="C:/Users/Wei Jia")
# set.seed(02032024)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(plotly)
library(akima)
library(gridExtra)
library(grid)
color_scheme_set("brightblue")

setwd("C:/Mixture")
source("BVCMIM/BVCMIM_Sim_Pred.R")
simulated_data = readRDS('BVCMIM/simulations/scenario_3/Simulated_data/5_rho/scenario_3_5_rho_1.RDS')
fit = readRDS('BVCMIM/simulations/scenario_4/results/5_rho/scenario_4_5_rho_1.RDS')
q = simulated_data$q
M = simulated_data$M 
Z = simulated_data$Z
Y = simulated_data$Y
X = simulated_data$X
U = simulated_data$U
W = simulated_data$W
n = simulated_data$n

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

weights_1_1_25 = as.matrix(apply(fit$draws('weights_1_1'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_1_2_25 = as.matrix(apply(fit$draws('weights_1_2'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_1_3_25 = as.matrix(apply(fit$draws('weights_1_3'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_1_4_25 = as.matrix(apply(fit$draws('weights_1_4'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_1_5_25 = as.matrix(apply(fit$draws('weights_1_5'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_2_1_25 = as.matrix(apply(fit$draws('weights_2_1'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_2_2_25 = as.matrix(apply(fit$draws('weights_2_2'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_2_3_25 = as.matrix(apply(fit$draws('weights_2_3'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_2_4_25 = as.matrix(apply(fit$draws('weights_2_4'), 3, function(x){quantile(x, probs=c(0.025))}))
weights_2_5_25 = as.matrix(apply(fit$draws('weights_2_5'), 3, function(x){quantile(x, probs=c(0.025))}))

weights_1_1_975 = as.matrix(apply(fit$draws('weights_1_1'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_1_2_975 = as.matrix(apply(fit$draws('weights_1_2'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_1_3_975 = as.matrix(apply(fit$draws('weights_1_3'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_1_4_975 = as.matrix(apply(fit$draws('weights_1_4'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_1_5_975 = as.matrix(apply(fit$draws('weights_1_5'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_2_1_975 = as.matrix(apply(fit$draws('weights_2_1'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_2_2_975 = as.matrix(apply(fit$draws('weights_2_2'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_2_3_975 = as.matrix(apply(fit$draws('weights_2_3'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_2_4_975 = as.matrix(apply(fit$draws('weights_2_4'), 3, function(x){quantile(x, probs=c(0.975))}))
weights_2_5_975 = as.matrix(apply(fit$draws('weights_2_5'), 3, function(x){quantile(x, probs=c(0.975))}))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
w_1_1 = range01(weights_1_1_mean)
df_1_1 = data.frame(z=c("z1", "z2", "z3", "z4", "z5", 'z6', 'z7', 'z8'), weight=w_1_1)
w_1_2 = range01(weights_1_2_mean)
df_1_2 = data.frame(z=c("z9", "z10", "z11", "z12", "z13", 'z14', 'z15', 'z16', 'z17', 'z18', 'z19', 'z20', 'z21', 'z22', 'z23', 'z24'), weight=w_1_2)
w_1_3 = range01(weights_1_3_mean)
df_1_3 = data.frame(z=c("z25", "z26"), weight=w_1_3)
w_1_4 = range01(weights_1_4_mean)
df_1_4 = data.frame(z=c("z27", "z28", 'z29', 'z30', 'z31', 'z32'), weight=w_1_4)
w_1_5 = range01(weights_1_5_mean)
df_1_5 = data.frame(z=c("z33", 'z34', 'z35', 'z36'), weight=w_1_5)
w_2_1 = range01(weights_2_1_mean)
df_2_1 = data.frame(z=c("z1", "z2", "z3", "z4", "z5", 'z6', 'z7', 'z8'), weight=w_2_1)
w_2_2 = range01(weights_2_2_mean)
df_2_2 = data.frame(z=c("z9", "z10", "z11", "z12", "z13", 'z14', 'z15', 'z16', 'z17', 'z18', 'z19', 'z20', 'z21', 'z22', 'z23', 'z24'), weight=w_2_2)
w_2_3 = range01(weights_2_3_mean)
df_2_3 = data.frame(z=c("z25", "z26"), weight=w_2_3)
w_2_4 = range01(weights_2_4_mean)
df_2_4 = data.frame(z=c("z27", "z28", 'z29', 'z30', 'z31', 'z32'), weight=w_2_4)
w_2_5 = range01(weights_2_5_mean)
df_2_5 = data.frame(z=c("z33", 'z34', 'z35', 'z36'), weight=w_2_5)
library(ggplot2)
# Basic barplot har
p_1_1<-ggplot(data=df_1_1[order(df_1_1$weight, decreasing = TRUE),], aes(x=z, y=weight)) +
  geom_bar(stat="identity") + coord_flip()
p_1_1
