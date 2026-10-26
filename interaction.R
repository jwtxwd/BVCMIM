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
fit = readRDS('BVCMIM/simulations/scenario_3/results/5_rho/scenario_3_5_rho_1.RDS')
q = simulated_data$q
M = simulated_data$M 
Z = simulated_data$Z
Y = simulated_data$Y
X = simulated_data$X
U = simulated_data$U
W = simulated_data$W
n = simulated_data$n

#   xlab("expos1")
qs = c(0.25, 0.75, 0.50)


#############################
hfun_1 <- function(z) z #dnorm(z-0.2,0,0.5) #(0.85*z)^2-0.5 # h1 <- 3/(1+exp(-2*wx1))
hfun_2 <- function(z) z
hfun_3 <- function(z) 1.5/(1+exp(-15*z)) #1.5/(1+exp(-15*z)) -0.25*z
group_1_ind = c(1:8)
group_2_ind = c(9:24)
group_3_ind = c(25,26)
group_4_ind = c(27:32)
group_5_ind = c(33:36)
ind_ls = list(group_1_ind, group_2_ind, group_3_ind, group_4_ind, group_5_ind)

###################### h1
library(plyr)      # for .(...)
vars <- c('1', '2', '3', '4', '5')    # explanatory variables= 
vars <- data.frame(t(expand.grid(vars,vars)))
gg <- do.call(rbind,lapply(vars,function(v){
  v   <- as.character(v)
  index_1_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
  index_1_z = index_1_z[rep(seq_len(nrow(index_1_z)), 50),]
  samp = runif(100, min = -3, max = 3)
  for (idx in eval(parse(text = paste(paste('group_', v[1], '_ind', sep = ''))))){
    index_1_z[,idx] 	= runif(100, min = -3, max = 3)
  }
  index_1_z_75 = index_1_z
  for (idx in eval(parse(text = paste(paste('group_', v[2], '_ind', sep = ''))))){
    index_1_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
  }
  index_1_z_25 = index_1_z
  for (idx in eval(parse(text = paste(paste('group_', v[2], '_ind', sep = ''))))){
    index_1_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
  }
  h1_75 = newh.postmean(Z = Z, Znew = index_1_z_75, fit = fit)$postmean[101:200]
  
  h1_25 = newh.postmean(Z = Z, Znew = index_1_z_25, fit = fit)$postmean[101:200]
  assign('weights', as.matrix(apply(fit$draws(paste('weights_1_', v[1], sep = '')), 3, mean)))
  assign('group_ind', eval(parse(text = paste(paste('group_', v[1], '_ind', sep = '')))))
  df = data.frame(y = c(h1_75, h1_25), x = c(index_1_z[,group_ind] %*% weights, index_1_z[,group_ind] %*% weights), q = c(rep(75, 100), rep(25, 100)))
  return(cbind(H=paste('Index',v[1]),V=paste('Index',v[2]),df))
}))
gg     <- data.frame(gg)                     # ggplot needs a data frame
labels <- aggregate(cbind(x,y)~H+V,gg,mean)  # labels for the diagonals

png(filename="C:/Mixture/BVCMIM/plots/s3_h1_ia.png", height=700,width=800)
ggplot(gg, aes(x, y))+
  geom_smooth(data = gg[gg$H!=gg$V,], aes(col = factor(q)), se=F)+
  # geom_tile(subset=.(as.numeric(H) < as.numeric(V)),aes(x,y,fill=pred),height=1,width=1)+
  geom_text(data=labels[labels$H==labels$V,],aes(x,y,label=H),size=8)+
  facet_grid(V~H,scales="free")+
  # scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
  # scale_fill_gradientn(colours=colorRampPalette(c("red","yellow","darkgreen"))(100))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(
    colour = "Quantile")
dev.off()

######################## h2
vars <- c('1', '2', '3', '4', '5')    # explanatory variables= 
vars <- data.frame(t(expand.grid(vars,vars)))
gg <- do.call(rbind,lapply(vars,function(v){
  v   <- as.character(v)
  index_1_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
  index_1_z = index_1_z[rep(seq_len(nrow(index_1_z)), 50),]
  samp = runif(100, min = -3, max = 3)
  for (idx in eval(parse(text = paste(paste('group_', v[1], '_ind', sep = ''))))){
    index_1_z[,idx] 	= runif(100, min = -3, max = 3)
  }
  index_1_z_75 = index_1_z
  for (idx in eval(parse(text = paste(paste('group_', v[2], '_ind', sep = ''))))){
    index_1_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
  }
  index_1_z_25 = index_1_z
  for (idx in eval(parse(text = paste(paste('group_', v[2], '_ind', sep = ''))))){
    index_1_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
  }
  h1_75 = newh.postmean(Z = Z, Znew = index_1_z_75, fit = fit)$postmean[1:100]
  
  h1_25 = newh.postmean(Z = Z, Znew = index_1_z_25, fit = fit)$postmean[1:100]
  assign('weights', as.matrix(apply(fit$draws(paste('weights_2_', v[1], sep = '')), 3, mean)))
  assign('group_ind', eval(parse(text = paste(paste('group_', v[1], '_ind', sep = '')))))
  df = data.frame(y = c(h1_75, h1_25), x = c(index_1_z[,group_ind] %*% weights, index_1_z[,group_ind] %*% weights), q = c(rep(75, 100), rep(25, 100)))
  return(cbind(H=paste('Index',v[1]),V=paste('Index',v[2]),df))
}))
gg     <- data.frame(gg)                     # ggplot needs a data frame
labels <- aggregate(cbind(x,y)~H+V,gg,mean)  # labels for the diagonals

png(filename="C:/Mixture/BVCMIM/plots/s3_h2_ia.png", height=700,width=800)
ggplot(gg, aes(x, y))+
  geom_smooth(data = gg[gg$H!=gg$V,], aes(col = factor(q)), se=F)+
  # geom_tile(subset=.(as.numeric(H) < as.numeric(V)),aes(x,y,fill=pred),height=1,width=1)+
  geom_text(data=labels[labels$H==labels$V,],aes(x,y,label=H),size=8)+
  facet_grid(V~H,scales="free")+
  # scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
  # scale_fill_gradientn(colours=colorRampPalette(c("red","yellow","darkgreen"))(100))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  labs(
    colour = "Quantile")
dev.off()
############ h1
######## index 1 ########
#### with index 2
index_1_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_1_z = index_1_z[rep(seq_len(nrow(index_1_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_1_ind){
  index_1_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_1_z_75 = index_1_z
for (idx in group_2_ind){
  index_1_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_1_z_25 = index_1_z
for (idx in group_2_ind){
  index_1_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_1_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_1_z_25, fit = fit)$postmean[101:200]
weights_1_1_mean = as.matrix(apply(fit$draws('weights_1_1'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_1_z[,group_1_ind] %*% weights_1_1_mean, index_1_z[,group_1_ind] %*% weights_1_1_mean), q = c(rep(75, 100), rep(25, 100)))
fig_1_2 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)
###### with index 3
index_1_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_1_z = index_1_z[rep(seq_len(nrow(index_1_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_1_ind){
  index_1_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_1_z_75 = index_1_z
for (idx in group_3_ind){
  index_1_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_1_z_25 = index_1_z
for (idx in group_3_ind){
  index_1_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_1_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_1_z_25, fit = fit)$postmean[101:200]
weights_1_1_mean = as.matrix(apply(fit$draws('weights_1_1'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_1_z[,group_1_ind] %*% weights_1_1_mean, index_1_z[,group_1_ind] %*% weights_1_1_mean), q = c(rep(75, 100), rep(25, 100)))
fig_1_3 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)
######## with index 4
index_1_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_1_z = index_1_z[rep(seq_len(nrow(index_1_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_1_ind){
  index_1_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_1_z_75 = index_1_z
for (idx in group_4_ind){
  index_1_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_1_z_25 = index_1_z
for (idx in group_4_ind){
  index_1_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_1_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_1_z_25, fit = fit)$postmean[101:200]
weights_1_1_mean = as.matrix(apply(fit$draws('weights_1_1'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_1_z[,group_1_ind] %*% weights_1_1_mean, index_1_z[,group_1_ind] %*% weights_1_1_mean), q = c(rep(75, 100), rep(25, 100)))
fig_1_4 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)
############ with index 5 
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
  geom_smooth(aes(col = factor(q)), se=F)
  # facet_grid(variable2 ~ variable1) +
  # ggtitle("h(expos1 | quantiles of expos2)") +
  # xlab("expos1")

######## index 2 ########
#### with index 1
index_2_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_2_z = index_2_z[rep(seq_len(nrow(index_2_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_2_ind){
  index_2_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_2_z_75 = index_2_z
for (idx in group_1_ind){
  index_2_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_2_z_25 = index_2_z
for (idx in group_1_ind){
  index_2_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_2_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_2_z_25, fit = fit)$postmean[101:200]
weights_1_2_mean = as.matrix(apply(fit$draws('weights_1_2'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_2_z[,group_2_ind] %*% weights_1_2_mean, index_2_z[,group_2_ind] %*% weights_1_2_mean), q = c(rep(75, 100), rep(25, 100)))
fig_2_1 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)
###### with index 3
index_2_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_2_z = index_2_z[rep(seq_len(nrow(index_2_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_2_ind){
  index_2_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_2_z_75 = index_2_z
for (idx in group_3_ind){
  index_2_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_2_z_25 = index_2_z
for (idx in group_3_ind){
  index_2_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_2_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_2_z_25, fit = fit)$postmean[101:200]
weights_1_2_mean = as.matrix(apply(fit$draws('weights_1_2'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_2_z[,group_2_ind] %*% weights_1_2_mean, index_2_z[,group_2_ind] %*% weights_1_2_mean), q = c(rep(75, 100), rep(25, 100)))
fig_2_3 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)
######## with index 4
index_2_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_2_z = index_2_z[rep(seq_len(nrow(index_2_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_2_ind){
  index_2_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_2_z_75 = index_2_z
for (idx in group_4_ind){
  index_2_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_2_z_25 = index_2_z
for (idx in group_4_ind){
  index_2_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_2_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_2_z_25, fit = fit)$postmean[101:200]
weights_1_2_mean = as.matrix(apply(fit$draws('weights_1_2'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_2_z[,group_2_ind] %*% weights_1_2_mean, index_2_z[,group_2_ind] %*% weights_1_2_mean), q = c(rep(75, 100), rep(25, 100)))
fig_2_4 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)
############ with index 5 
index_2_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_2_z = index_2_z[rep(seq_len(nrow(index_2_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_2_ind){
  index_2_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_2_z_75 = index_2_z
for (idx in group_5_ind){
  index_2_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_2_z_25 = index_2_z
for (idx in group_5_ind){
  index_2_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_2_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_2_z_25, fit = fit)$postmean[101:200]
weights_1_2_mean = as.matrix(apply(fit$draws('weights_1_2'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_2_z[,group_2_ind] %*% weights_1_2_mean, index_2_z[,group_2_ind] %*% weights_1_2_mean), q = c(rep(75, 100), rep(25, 100)))
fig_2_5 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

######## index 3 ########
#### with index 1
index_3_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_3_z = index_3_z[rep(seq_len(nrow(index_3_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_3_ind){
  index_3_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_3_z_75 = index_3_z
for (idx in group_1_ind){
  index_3_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_3_z_25 = index_3_z
for (idx in group_1_ind){
  index_3_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_3_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_3_z_25, fit = fit)$postmean[101:200]
weights_1_3_mean = as.matrix(apply(fit$draws('weights_1_3'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_3_z[,group_3_ind] %*% weights_1_3_mean, index_3_z[,group_3_ind] %*% weights_1_3_mean), q = c(rep(75, 100), rep(25, 100)))
fig_3_1 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 2
index_3_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_3_z = index_3_z[rep(seq_len(nrow(index_3_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_3_ind){
  index_3_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_3_z_75 = index_3_z
for (idx in group_2_ind){
  index_3_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_3_z_25 = index_3_z
for (idx in group_2_ind){
  index_3_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_3_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_3_z_25, fit = fit)$postmean[101:200]
weights_1_3_mean = as.matrix(apply(fit$draws('weights_1_3'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_3_z[,group_3_ind] %*% weights_1_3_mean, index_3_z[,group_3_ind] %*% weights_1_3_mean), q = c(rep(75, 100), rep(25, 100)))
fig_3_2 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 4
index_3_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_3_z = index_3_z[rep(seq_len(nrow(index_3_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_3_ind){
  index_3_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_3_z_75 = index_3_z
for (idx in group_4_ind){
  index_3_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_3_z_25 = index_3_z
for (idx in group_4_ind){
  index_3_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_3_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_3_z_25, fit = fit)$postmean[101:200]
weights_1_3_mean = as.matrix(apply(fit$draws('weights_1_3'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_3_z[,group_3_ind] %*% weights_1_3_mean, index_3_z[,group_3_ind] %*% weights_1_3_mean), q = c(rep(75, 100), rep(25, 100)))
fig_3_4 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 5
index_3_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_3_z = index_3_z[rep(seq_len(nrow(index_3_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_3_ind){
  index_3_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_3_z_75 = index_3_z
for (idx in group_5_ind){
  index_3_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_3_z_25 = index_3_z
for (idx in group_5_ind){
  index_3_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_3_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_3_z_25, fit = fit)$postmean[101:200]
weights_1_3_mean = as.matrix(apply(fit$draws('weights_1_3'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_3_z[,group_3_ind] %*% weights_1_3_mean, index_3_z[,group_3_ind] %*% weights_1_3_mean), q = c(rep(75, 100), rep(25, 100)))
fig_3_5 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

######## index 4 ########
#### with index 1
index_4_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_4_z = index_4_z[rep(seq_len(nrow(index_4_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_4_ind){
  index_4_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_4_z_75 = index_4_z
for (idx in group_1_ind){
  index_4_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_4_z_25 = index_4_z
for (idx in group_1_ind){
  index_4_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_4_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_4_z_25, fit = fit)$postmean[101:200]
weights_1_4_mean = as.matrix(apply(fit$draws('weights_1_4'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_4_z[,group_4_ind] %*% weights_1_4_mean, index_4_z[,group_4_ind] %*% weights_1_4_mean), q = c(rep(75, 100), rep(25, 100)))
fig_4_1 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)


#### with index 2
index_4_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_4_z = index_4_z[rep(seq_len(nrow(index_4_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_4_ind){
  index_4_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_4_z_75 = index_4_z
for (idx in group_2_ind){
  index_4_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_4_z_25 = index_4_z
for (idx in group_2_ind){
  index_4_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_4_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_4_z_25, fit = fit)$postmean[101:200]
weights_1_4_mean = as.matrix(apply(fit$draws('weights_1_4'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_4_z[,group_4_ind] %*% weights_1_4_mean, index_4_z[,group_4_ind] %*% weights_1_4_mean), q = c(rep(75, 100), rep(25, 100)))
fig_4_2 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 2
index_4_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_4_z = index_4_z[rep(seq_len(nrow(index_4_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_4_ind){
  index_4_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_4_z_75 = index_4_z
for (idx in group_2_ind){
  index_4_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_4_z_25 = index_4_z
for (idx in group_2_ind){
  index_4_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_4_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_4_z_25, fit = fit)$postmean[101:200]
weights_1_4_mean = as.matrix(apply(fit$draws('weights_1_4'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_4_z[,group_4_ind] %*% weights_1_4_mean, index_4_z[,group_4_ind] %*% weights_1_4_mean), q = c(rep(75, 100), rep(25, 100)))
fig_4_2 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 3
index_4_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_4_z = index_4_z[rep(seq_len(nrow(index_4_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_4_ind){
  index_4_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_4_z_75 = index_4_z
for (idx in group_3_ind){
  index_4_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_4_z_25 = index_4_z
for (idx in group_3_ind){
  index_4_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_4_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_4_z_25, fit = fit)$postmean[101:200]
weights_1_4_mean = as.matrix(apply(fit$draws('weights_1_4'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_4_z[,group_4_ind] %*% weights_1_4_mean, index_4_z[,group_4_ind] %*% weights_1_4_mean), q = c(rep(75, 100), rep(25, 100)))
fig_4_3 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 5
index_4_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_4_z = index_4_z[rep(seq_len(nrow(index_4_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_4_ind){
  index_4_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_4_z_75 = index_4_z
for (idx in group_5_ind){
  index_4_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_4_z_25 = index_4_z
for (idx in group_5_ind){
  index_4_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_4_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_4_z_25, fit = fit)$postmean[101:200]
weights_1_4_mean = as.matrix(apply(fit$draws('weights_1_4'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_4_z[,group_4_ind] %*% weights_1_4_mean, index_4_z[,group_4_ind] %*% weights_1_4_mean), q = c(rep(75, 100), rep(25, 100)))
fig_4_5 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

######## index 5 ########
#### with index 1
index_5_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_5_z = index_5_z[rep(seq_len(nrow(index_5_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_5_ind){
  index_5_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_5_z_75 = index_5_z
for (idx in group_1_ind){
  index_5_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_5_z_25 = index_5_z
for (idx in group_1_ind){
  index_5_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_5_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_5_z_25, fit = fit)$postmean[101:200]
weights_1_5_mean = as.matrix(apply(fit$draws('weights_1_5'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_5_z[,group_5_ind] %*% weights_1_5_mean, index_5_z[,group_5_ind] %*% weights_1_5_mean), q = c(rep(75, 100), rep(25, 100)))
fig_5_1 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 2
index_5_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_5_z = index_5_z[rep(seq_len(nrow(index_5_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_5_ind){
  index_5_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_5_z_75 = index_5_z
for (idx in group_2_ind){
  index_5_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_5_z_25 = index_5_z
for (idx in group_2_ind){
  index_5_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_5_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_5_z_25, fit = fit)$postmean[101:200]
weights_1_5_mean = as.matrix(apply(fit$draws('weights_1_5'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_5_z[,group_5_ind] %*% weights_1_5_mean, index_5_z[,group_5_ind] %*% weights_1_5_mean), q = c(rep(75, 100), rep(25, 100)))
fig_5_2 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 3
index_5_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_5_z = index_5_z[rep(seq_len(nrow(index_5_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_5_ind){
  index_5_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_5_z_75 = index_5_z
for (idx in group_3_ind){
  index_5_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_5_z_25 = index_5_z
for (idx in group_3_ind){
  index_5_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_5_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_5_z_25, fit = fit)$postmean[101:200]
weights_1_5_mean = as.matrix(apply(fit$draws('weights_1_5'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_5_z[,group_5_ind] %*% weights_1_5_mean, index_5_z[,group_5_ind] %*% weights_1_5_mean), q = c(rep(75, 100), rep(25, 100)))
fig_5_3 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)

#### with index 4
index_5_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_5_z = index_5_z[rep(seq_len(nrow(index_5_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_5_ind){
  index_5_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_5_z_75 = index_5_z
for (idx in group_4_ind){
  index_5_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_5_z_25 = index_5_z
for (idx in group_4_ind){
  index_5_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h1_75 = newh.postmean(Z = Z, Znew = index_5_z_75, fit = fit)$postmean[101:200]

h1_25 = newh.postmean(Z = Z, Znew = index_5_z_25, fit = fit)$postmean[101:200]


weights_1_5_mean = as.matrix(apply(fit$draws('weights_1_5'), 3, mean))
df = data.frame(y = c(h1_75, h1_25), x = c(index_5_z[,group_5_ind] %*% weights_1_5_mean, index_5_z[,group_5_ind] %*% weights_1_5_mean), q = c(rep(75, 100), rep(25, 100)))
fig_5_4 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)), se=F)  + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    # axis.title.y = element_blank
  )

empty_plot = ggplot(mtcars, aes(wt, mpg))

# Nothing to see here!

empty_plot()

fig <- subplot(empty, fig_1_2, fig_1_3, fig_1_4, fig_1_5, fig_2_1, empty_plot(), fig_2_3, fig_2_4, fig_2_5, 
               fig_3_1, fig_3_2, empty_plot(), fig_3_4, fig_3_5, 
               fig_4_1, fig_4_2, fig_4_3, empty_plot(), fig_4_5,
               fig_5_1, fig_5_2, fig_5_3, fig_5_4, empty_plot(), nrows = 5, shareY = TRUE) %>% 
  layout(title = list(text = "Multiple Subplots with Shared Y-Axes"))
         # plot_bgcolor='#e5ecf6', 
         # xaxis = list( 
         #   zerolinecolor = '#ffff', 
         #   zerolinewidth = 2, 
         #   gridcolor = 'ffff'), 
         # yaxis = list( 
         #   zerolinecolor = '#ffff', 
         #   zerolinewidth = 2, 
         #   gridcolor = 'ffff')) 
fig
############### h2
index_2_z = rbind(apply(Z, 2, median), apply(Z, 2, median))
index_2_z = index_2_z[rep(seq_len(nrow(index_2_z)), 50),]
samp = runif(100, min = -3, max = 3)
for (idx in group_3_ind){
  index_2_z[,idx] 	= runif(100, min = -3, max = 3)
}
index_2_z_75 = index_2_z
for (idx in group_4_ind){
  index_2_z_75[,idx] = rep(quantile(Z[,idx], qs[2]), 100)
}
index_2_z_25 = index_2_z
for (idx in group_4_ind){
  index_2_z_25[,idx] = rep(quantile(Z[,idx], qs[1]), 100)
}
h2_75 = newh.postmean(Z = Z, Znew = index_2_z_75, fit = fit)$postmean[1:100]

h2_25 = newh.postmean(Z = Z, Znew = index_2_z_25, fit = fit)$postmean[1:100]
weights_2_2_mean = as.matrix(apply(fit$draws('weights_2_2'), 3, mean))
weights_2_3_mean = as.matrix(apply(fit$draws('weights_2_3'), 3, mean))
weights_2_4_mean = as.matrix(apply(fit$draws('weights_2_4'), 3, mean))
df = data.frame(y = c(h2_75, h2_25), x = c(index_2_z[,group_3_ind] %*% weights_2_3_mean, index_2_z[,group_3_ind] %*% weights_2_3_mean), q = c(rep(75, 100), rep(25, 100)))
fig_2_4 = ggplot(df, aes(x, y)) +
  geom_smooth(aes(col = factor(q)))+
  # facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")


fig1 <- ggplotly(fig_2_4)
fig2 <- plot_ly(x = c(20, 30, 40), y = c(5, 5, 5), type = 'scatter', mode = 'lines+markers') 
fig3 <- plot_ly(x = c(2, 3, 4), y = c(600, 700, 800), type = 'scatter', mode = 'lines+markers') 
fig4 <- plot_ly(x = c(4000, 5000, 6000), y = c(7000, 8000, 9000), type = 'scatter', mode = 'lines+markers') 
fig <- subplot(fig1, fig2, fig3, fig4, nrows = 2, shareY = TRUE) %>% 
  layout(title = list(text = "Multiple Subplots with Shared Y-Axes"))
         # plot_bgcolor='#e5ecf6', 
         # xaxis = list( 
         #   zerolinecolor = '#ffff', 
         #   zerolinewidth = 2, 
         #   gridcolor = 'ffff'), 
         # yaxis = list( 
         #   zerolinecolor = '#ffff', 
         #   zerolinewidth = 2, 
         #   gridcolor = 'ffff')) 
fig