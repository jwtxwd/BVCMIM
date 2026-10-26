########################
# Data simulation code #
########################
setwd("C:/Users/Wei Jia/OneDrive - University of Cincinnati/Mixture")
library(MASS); library(stats)
# set.seed(02025024)
chem_mix = read.csv('preprocessed/data_menmdi.csv')[,c('MBP',	'MIBP',	'MEP2',	'MBZP2'
                                                       ,'MCPP',	'MEHHP',	'MECPP',	'BPA',	'PCB28',	'PCB74'
                                                       ,'PCB99',	'PCB105',	'PCB118',	'PCB138_158',	'PCB146'
                                                       ,'PCB153',	'PCB156',	'PCB170',	'PCB183',	'PCB187'
                                                       ,'PCB194','PCB196_203',	'PCB199',	'PCB206',	'PP_DDE'
                                                       ,'T_NONA',	'BB153',	'PBDE28',	'PBDE47',	'PBDE99',	'PBDE100'
                                                       ,'PBDE153',	'PFOA',	'PFOS',	'PFNA',	'PFHXS')]
Sigma_Z = matrix(cor(chem_mix), nrow = 36)
data_chem_long = read.csv('preprocessed/data_menmdi_long.csv')
data_chem_long = data_chem_long[order(data_chem_long$participant_id), ]
# Sigma_Z = matrix(cor(chem_mix[,c(1,9,13,21,33,3,16,25,28,35)]), nrow = 10)
# Sigma_Z = matrix(cor(chem_mix[,c(1,9,13,21,33)]), nrow = 5)
n 				= 149 # Number of subjects
q 				= 2 # Random slope and random intercept
p 				= 3 # Number of confounders
M				  = 36 # Number of metals
age 			= c(0,1,2)
# age = c(0,1,2)
T				  = length(age) # Number of time points
N				  = n*T

#############################


# Z 				= matrix(mvrnorm(n = n, mu=c(rep(0,M)), Sigma = matrix(data = c(1,0.2,0.3,0.45,0.15, 0.2,1,0.4,0.45,0.3, 0.3,0.4,1,0.3,0.15, 0.45,0.45,0.3,1,0.2, 0.15,0.3,0.15,0.2,1), nrow=M)), byrow=FALSE, nrow=n)

Z 				= as.matrix(chem_mix)

# X = as.matrix(data_chem_long[,c('ga', 'b_wght')])
# X 				= apply(X, 2, scale)
gender=rep(0, length(data_chem_long[,c('GENDER')]))
for (i in 1:length(data_chem_long[,c('GENDER')])){
  if (data_chem_long[, c('GENDER')][i]=='Male'){
    gender[i]=1
  }else{gender[i]=0}
}
# X = as.matrix(cbind(X, as.matrix(gender)))
X = as.matrix(gender)
year=rep(0, length(data_chem_long[,c('age')]))
for (i in 1:length(data_chem_long[,c('age')])){
  if (data_chem_long[, c('age')][i]=='menmdi_M12'){
    year[i]=0
  }else if (data_chem_long[, c('age')][i]=='menmdi_M24') {
     year[i]=1
  }else{year[i]=2}
}
X = as.matrix(cbind(X, as.matrix(year)))
intercept = rep(1, length(data_chem_long[,c('age')]))
X = as.matrix(cbind(X, as.matrix(intercept)))
beta			= rep(0.12,p)

b				= c( t(mvrnorm(n = n, mu=rep(0,q), Sigma = matrix(data = c(0.25,0.005,0.005,0.05), nrow=q))))

res.sd = 1

# Make U matrix

U				= matrix(data=rep(0,T*2*n), nrow=T)
U[,1]			= rep(1, T)
U[,2]			= age

counter 		= 1
while(counter < n) {
    U_i					= matrix(data=rep(0,T*2*n), nrow=T)
    U_i[,counter*2+1]	= rep(1, T)
    U_i[,counter*2+2]	= age
    U 					= rbind(U, U_i)
    counter 				= counter+1
}

Z_quant = Z
#Create quantile index
for (m in 1:M){
  quants = quantile(Z[,m])
  for (i in 1:n){
    if (Z[i,m]<quants[2]){
      Z_quant[i,m] = 0
    }else if (Z[i, m]<quants[3]){
      Z_quant[i, m] = 1
    }else if (Z[i, m]<quants[4]){
      Z_quant[i, m] = 2
    }else {Z_quant[i, m] = 3}
  }
}

# #scenario 3
hfun_1 <- function(z) 1.5/(1+exp(-15*z))
weights_1_1 = c(0.5, 0.6, 0.5, 0, 0, 0 , 0, 0)
# weights_1_1 <- weights_1_1/sqrt(sum(weights_1_1^2))
weights_1_5 = c(0, 0, 0.5, 0.5)
# weights_1_5 <- weights_1_5/sqrt(sum(weights_1_5^2))
index_1_1 = Z[,1:8] %*% weights_1_1
index_1_1 = scale(index_1_1)
index_1_5 = Z[,33:36] %*% weights_1_5
index_1_5 = scale(index_1_5)
h_1 = index_1_1 + hfun_1(index_1_5) + 2*index_1_1*hfun_1(index_1_5) 

weights_2_2 = c(0.6, 0.4, 0.5, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0)
# weights_2_1 <- weights_2_1/sqrt(sum(weights_2_1^2))
weights_2_4 = c(0, 0, 0, 0, 0.4, 0.2)
# weights_2_4 <- weights_2_4/sqrt(sum(weights_2_4^2))
index_2_2 = Z[,9:24] %*% weights_2_2
index_2_2 = scale(index_2_2)
index_2_4 = Z[,27:32] %*% weights_2_4
index_2_4 = scale(index_2_4)
h_2 = index_2_2 + hfun_1(index_2_4) + 2*index_2_2*hfun_1(index_2_4)
# Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + rep((Z[,1:8] %*% weights_1_1 + Z[,33:36] %*% weights_1_5 + 0.5*(Z[,1:8] %*% weights_1_1)*(Z[,33:36] %*% weights_1_5)), each=T) + matrix(rep((Z[,1:8] %*% weights_2_1 + Z[,27:32] %*% weights_2_4 + 0.5*(Z[,1:8] %*% weights_2_1)*(Z[,27:32] %*% weights_2_4)), each=T)) * rep(age, n) + U %*% b
Y = as.matrix(as.numeric(data_chem_long[,'outcome']))
# Make W matrix

W				= matrix(data=rep(0,T*2*n), nrow=T)
W[,1]			= rep(1, T)
W[,n+1]			= age

counter 		= 2
while(counter <= n) {
	W_i				= matrix(data=rep(0,T*2*n), nrow=T)
	W_i[,counter]	= rep(1, T)
	W_i[,n+counter]	= age
	W 				= rbind(W, W_i)
	counter 		= counter+1
}

WTW = t(W) %*% W
XTX	= t(X) %*% X
UTU = t(U) %*% U
