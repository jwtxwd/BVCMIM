Sys.setenv(HOME="C:/Users/Wei Jia")
# set.seed(02032024)
set.seed(02025024)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(plotly)
library(akima)
color_scheme_set("brightblue")

setwd("C:/Mixture")
source("BVCMIM/BVCMIM_Sim_Pred.R")
simulated_data = readRDS('BVCMIM/Simulated_data/scenario_2.RDS')
fit = readRDS('BVCMIM/results/scenario_2.RDS')
q = simulated_data$q
M = simulated_data$M 
Z = simulated_data$Z
Y = simulated_data$Y
X = simulated_data$X
U = simulated_data$U
W = simulated_data$W
n = simulated_data$n

index_1 = rep(0, 1000)
index_2 = rep(0, 1000)
index_3 = rep(0, 1000)
index_4 = runif(1000, -1.6, 1.6)
index_5 = runif(1000, -1.6, 1.6)

index_1 = runif(1000, -1.6, 1.6)
index_2 = runif(1000, -1.6, 1.6)
index_3 = rep(0, 1000)
index_4 = rep(0, 1000)
index_5 = rep(0, 1000)

index_new = cbind(index_1, index_2, index_3, index_4, index_5)

hgrid_cross_sec = newh_postmean_index(Z, index_new, fit)
# par(mfrow=c(1,2))
# seq(-1.5, 1.5, 0.01)
df = data.frame(h = hgrid_cross_sec$postmean[1:length(index_4)],
                x = round(index_1, 1),
                y = round(index_2, 1))



resolution <- 0.01 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
a <- interp(x=df$x, y=df$y, z=df$h, 
            xo=seq(min(df$x),max(df$x),by=resolution), 
            yo=seq(min(df$y),max(df$y),by=resolution), duplicate="mean")

image(a, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))


