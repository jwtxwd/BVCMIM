########################
# Data simulation code #
########################

library(MASS); library(stats)
set.seed(02025024)
chem_mix = read.csv('preprocessed/data_motpdi.csv')[,c('MBP',	'MIBP',	'MEP2',	'MBZP2'
                                                       ,'MCPP',	'MEHHP',	'MECPP',	'BPA',	'PCB28',	'PCB74'
                                                       ,'PCB99',	'PCB105',	'PCB118',	'PCB138_158',	'PCB146'
                                                       ,'PCB153',	'PCB156',	'PCB170',	'PCB183',	'PCB187'
                                                       ,'PCB194','PCB196_203',	'PCB199',	'PCB206',	'PP_DDE'
                                                       ,'T_NONA',	'BB153',	'PBDE28',	'PBDE47',	'PBDE99',	'PBDE100'
                                                       ,'PBDE153',	'PFOA',	'PFOS',	'PFNA',	'PFHXS')]
Sigma_Z = matrix(cor(chem_mix), nrow = 36)
# Sigma_Z = matrix(cor(chem_mix[,c(1,9,13,21,33,3,16,25,28,35)]), nrow = 10)
# Sigma_Z = matrix(cor(chem_mix[,c(1,9,13,21,33)]), nrow = 5)
n 				= 100 # Number of subjects
q 				= 2 # Random slope and random intercept
p 				= 2 # Number of confounders
M				  = 36 # Number of metals
age 			= c(-2:2)
# age = c(0,1,2)
T				  = length(age) # Number of time points
N				  = n*T

#############################


# Z 				= matrix(mvrnorm(n = n, mu=c(rep(0,M)), Sigma = matrix(data = c(1,0.2,0.3,0.45,0.15, 0.2,1,0.4,0.45,0.3, 0.3,0.4,1,0.3,0.15, 0.45,0.45,0.3,1,0.2, 0.15,0.3,0.15,0.2,1), nrow=M)), byrow=FALSE, nrow=n)
Z 				= matrix(mvrnorm(n = n, mu=c(rep(0,M)), Sigma = Sigma_Z), byrow=FALSE, nrow=n)

Z 				= apply(Z, 2, scale)

X				= cbind( matrix( rep(rnorm(n, mean=0, sd=1), each=T), byrow=FALSE ), matrix( rep(scale(matrix(sample(1:2, n, replace = TRUE))), each=T) ))
beta			= rep(1,p)

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

# #scenario 1
# Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + rep(5*(0.1*Z_quant[, 1] + 0.3*Z_quant[, 12] + 0.4*Z_quant[,24] + 0.2*Z_quant[,35]), each=T) + matrix(rep(3*(0.1*Z_quant[, 9] + 0.5*Z_quant[, 23] + 0.2*Z_quant[,27] + 0.2*Z_quant[,33]), each=T)) * rep(age, n) + U %*% b
#scenario 2
# Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + rep(5*(0.1*Z_quant[, 1] + 0.3*Z_quant[, 12] - 0.3*Z_quant[, 23] + 0.4*Z_quant[,24] + 0.2*Z_quant[,35]), each=T) + matrix(rep(3*(-0.2*Z_quant[, 6] + 0.1*Z_quant[, 9] + 0.5*Z_quant[, 23] + 0.2*Z_quant[,27] + 0.2*Z_quant[,33]), each=T)) * rep(age, n) + U %*% b
#scenario 3
# Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + rep(5*(0.1*Z[, 1] + 0.3*Z[, 12] + 0.4*Z[,24] + 0.2*Z[,35]), each=T) + matrix(rep(3*(0.1*Z[, 9] + 0.5*Z[, 23] + 0.2*Z[,27] + 0.2*Z[,33]), each=T)) * rep(age, n) + U %*% b
# #scenario 4
Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + rep(5*(0.1*Z[, 1] + 0.3*Z[, 12] - 0.3*Z[, 23] + 0.4*Z[,24] + 0.2*Z[,35]), each=T) + matrix(rep(3*(-0.2*Z[, 2] + 0.1*Z[, 9] + 0.5*Z[, 23] + 0.2*Z[,27] + 0.2*Z[,33]), each=T)) * rep(age, n) + U %*% b

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
