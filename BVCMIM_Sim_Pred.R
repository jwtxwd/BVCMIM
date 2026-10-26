#####################################
# Prediction code for new exposures #
# Postmean estimation functions     #
#####################################
library(magic)
library(MASS)
# rbf = function(x, y, eta, rho){
#   K = matrix(0, dim(x)[1], dim(y)[1])
#   for (i in 1:dim(x)[1]){
#     for (j in 1:dim(y)[1]){
#       K[i,j] = eta*exp(-rho*sum((x[i,]-y[j,])^2))
#       # K[i,j] = (x[i,] %*% y[j,]+1)^2;
#     }
#   }
#   return(K)
# }
rbf = function(x, y, eta, rho){
  K = matrix(0, dim(x)[1], dim(y)[1])
  for (i in 1:dim(x)[1]){
    for (j in 1:dim(y)[1]){
      # K[i,j] = eta*exp(-rho*sum((x[i,]-y[j,])^2))
      K[i,j] = eta*exp(-((x[i,]-y[j,])^2)%*%rho)
      # K[i,j] = (x[i,] %*% y[j,]+1)^2;
    }
  }
  return(K)
}
# newh.postmean <- function(Z, Znew, fit){
#   Z_ = as.matrix(rbind(Z, Znew))
#   n = nrow(Z) + nrow(Znew)
#   eta_1	    = as.data.frame(fit$draws('eta_1'))
#   rho_1 = as.data.frame(fit$draws('rho_1'))
#   weights_1 = as.data.frame(fit$draw())
#   eta_2	    = as.data.frame(fit$draws('eta_2'))
#   rho_2 = as.data.frame(fit$draws('rho_2'))
#   
#   h_df = c()
#   for (i in 1:1000){
#     list_K				  = list()
#     list_K[[1]] 	= (rbf(Z_,Z_,eta_1[i,],rho_1[i,]) + diag(n)*1e-7)
#     list_K[[2]] = (rbf(Z_,Z_,eta_2[i,],rho_2[i,]) + diag(n)*1e-7)
#     
#     Sigma_h 				= do.call(adiag, list_K)
#   
#     h = mvrnorm(n=2*n, mu = rep(0, 2*n), Sigma = Sigma_h)
#     h_df = cbind(h_df, h)
#     
#   }
#   return(h_df)
# }
newh.postmean <- function(Z, Znew, fit) {


  if(is.null(dim(Znew))) Znew = matrix(Znew, nrow=1)
  # if(class(Znew)  == "data.frame") Znew <- data.matrix(Znew)
  # Znew            = as.matrix(Znew)
  pnew            = nrow(Znew)

  Sigsq_mean	    = mean(fit$draws('sigma_sq'))
  # Tau_mean	      = as.matrix(apply(fit$draws('tau_sq'), 2, mean))
  b_mean		      = as.matrix(apply(fit$draws('b'), 3, mean))
  Beta_mean	      = as.matrix(apply(fit$draws('Beta'), 3, mean))
  h_mean		      = as.matrix(apply(fit$draws('h'), 3, mean))
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
  eta_1_mean = as.matrix(apply(fit$draws('eta_1'), 3, mean))
  eta_2_mean = as.matrix(apply(fit$draws('eta_2'), 3, mean))
  rho_1_mean = as.matrix(apply(fit$draws('rho_1'), 3, mean))
  rho_2_mean = as.matrix(apply(fit$draws('rho_2'), 3, mean))
  Y			          = Y
  X 			        = X
  U			          = U
  W               = W
  
  index_1_1 = Z[,1:8] %*% weights_1_1_mean
  index_1_1_norm = (index_1_1-mean(index_1_1))/sd(index_1_1)
  index_1_2 = Z[,9:24] %*% weights_1_2_mean
  index_1_2_norm = (index_1_2-mean(index_1_2))/sd(index_1_2)
  index_1_3 = Z[,25:26] %*% weights_1_3_mean
  index_1_3_norm = (index_1_3-mean(index_1_3))/sd(index_1_3)
  index_1_4 = Z[,27:32] %*% weights_1_4_mean
  index_1_4_norm = (index_1_4-mean(index_1_4))/sd(index_1_4)
  index_1_5 = Z[,33:36] %*% weights_1_5_mean
  index_1_5_norm = (index_1_5-mean(index_1_5))/sd(index_1_5)
  index_1 = cbind(index_1_1_norm, index_1_2_norm, index_1_3_norm, index_1_4_norm, index_1_5_norm)
  
  index_2_1 = Z[,1:8] %*% weights_2_1_mean
  index_2_1_norm = (index_2_1-mean(index_2_1))/sd(index_2_1)
  index_2_2 = Z[,9:24] %*% weights_2_2_mean
  index_2_2_norm = (index_2_2-mean(index_2_2))/sd(index_2_2)
  index_2_3 = Z[,25:26] %*% weights_2_3_mean
  index_2_3_norm = (index_2_3-mean(index_2_3))/sd(index_2_3)
  index_2_4 = Z[,27:32] %*% weights_2_4_mean
  index_2_4_norm = (index_2_4-mean(index_2_4))/sd(index_2_4)
  index_2_5 = Z[,33:36] %*% weights_2_5_mean
  index_2_5_norm = (index_2_5-mean(index_2_5))/sd(index_2_5)
  index_2 = cbind(index_2_1_norm, index_2_2_norm, index_2_3_norm, index_2_4_norm, index_2_5_norm)
  
  index_1_1_new = Znew[,1:8] %*% weights_1_1_mean
  index_1_1_new_norm = (index_1_1_new-mean(index_1_1))/sd(index_1_1)
  index_1_2_new = Znew[,9:24] %*% weights_1_2_mean
  index_1_2_new_norm = (index_1_2_new-mean(index_1_2))/sd(index_1_2)
  index_1_3_new = Znew[,25:26] %*% weights_1_3_mean
  index_1_3_new_norm = (index_1_3_new-mean(index_1_3))/sd(index_1_3)
  index_1_4_new = Znew[,27:32] %*% weights_1_4_mean
  index_1_4_new_norm = (index_1_4_new-mean(index_1_4))/sd(index_1_4)
  index_1_5_new = Znew[,33:36] %*% weights_1_5_mean
  index_1_5_new_norm = (index_1_5_new-mean(index_1_5))/sd(index_1_5)
  index_1_new = cbind(index_1_1_new_norm, index_1_2_new_norm, index_1_3_new_norm, index_1_4_new_norm, index_1_5_new_norm)
  
  index_2_1_new = Znew[,1:8] %*% weights_2_1_mean
  index_2_1_new_norm = (index_2_1_new-mean(index_2_1))/sd(index_2_1)
  index_2_2_new = Znew[,9:24] %*% weights_2_2_mean
  index_2_2_new_norm = (index_2_2_new-mean(index_2_2))/sd(index_2_2)
  index_2_3_new = Znew[,25:26] %*% weights_2_3_mean
  index_2_3_new_norm = (index_2_3_new-mean(index_2_3))/sd(index_2_3)
  index_2_4_new = Znew[,27:32] %*% weights_2_4_mean
  index_2_4_new_norm = (index_2_4_new-mean(index_2_4))/sd(index_2_4)
  index_2_5_new = Znew[,33:36] %*% weights_2_5_mean
  index_2_5_new_norm = (index_2_5_new-mean(index_2_5))/sd(index_2_5)
  index_2_new = cbind(index_2_1_new_norm, index_2_2_new_norm, index_2_3_new_norm, index_2_4_new_norm, index_2_5_new_norm)

  list.K				  = list()

  # poly 				    = polydot(degree=2, offset=1)
  # poly 				    = rbfdot(sigma=1)

  list.K[[1]] 	= (rbf(index_1,index_1,eta_1_mean,rho_1_mean) + diag(n)*1e-7)
  list.K[[2]] = (rbf(index_2,index_2,eta_2_mean,rho_2_mean) + diag(n)*1e-7)

  Sigma.h 				= do.call(adiag, list.K)

  Sigma.h.new = "[<-"(matrix(0, (q*n+q*pnew), (q*n+q*pnew)), 1:nrow(Sigma.h), 1:ncol(Sigma.h), value = Sigma.h) # adding rows/cols of zeroes to square matrix

  list.G        = list()

  list.G[[1]] 	= (rbf(rbind(index_1, index_1_new),rbind(index_1, index_1_new),eta_1_mean,rho_1_mean) + diag(n+dim(Znew)[1])*1e-7)
  list.G[[2]] = (rbf(rbind(index_2, index_2_new),rbind(index_2, index_2_new),eta_2_mean,rho_2_mean) + diag(n+dim(Znew)[1])*1e-7)
  # h vector is ordered c(h1, h2, h2new, h1new)
  G.1.11 = list.G[[1]][ 1:n, 1:n ]
  G.1.22 = list.G[[1]][ (n+1):(n+pnew), (n+1):(n+pnew) ]
  G.1.21 = list.G[[1]][ (n+1):(n+pnew), 1:n ]
  G.1.12 = list.G[[1]][ 1:n, (n+1):(n+pnew) ]

  list.G.final = list()
  list.G.final[[1]] = G.1.11
  list.G.final[[2]] = list.G[[2]]
  list.G.final[[3]] = G.1.22
  cov.bf = do.call(adiag, list.G.final)

  cov.bf[(q*n + pnew + 1):(q*n + q*pnew), 1:n] = G.1.21
  cov.bf[1:n, (q*n + pnew + 1):(q*n + q*pnew)] = G.1.12

  Sigma11 = cov.bf[1:(n*q), 1:(n*q)]
  Sigma12 = cov.bf[1:(n*q), (n*q + 1):(n*q + pnew*q)]
  Sigma21 = cov.bf[(n*q + 1):(n*q + pnew*q), 1:(n*q)]
  Sigma22 = cov.bf[(n*q + 1):(n*q + pnew*q), (n*q + 1):(n*q + pnew*q)]

  WTW         = t(W) %*% W
  Sigma.h.inv = solve(Sigma.h + 1e-7*diag(nrow(Sigma.h)))
  
  # Sigma.h.inv = solve(Sigma.h + diag(nrow(Sigma.h)))
  m.h         = solve( WTW / Sigsq_mean + Sigma.h.inv) %*% t(W) %*% (Y - X %*% Beta_mean - U %*% b_mean) * (1/Sigsq_mean)
  s.h         = solve( WTW / Sigsq_mean + Sigma.h.inv)
  Sigma11.inv = Sigma.h.inv

  Mu.hnew.updated = Sigma21 %*% Sigma11.inv %*% m.h
  Sigma.hnew.updated = Sigma22 - Sigma21 %*% Sigma11.inv %*% Sigma12 + Sigma21 %*% Sigma11.inv %*% s.h %*% Sigma11.inv %*% Sigma12
  list(postmean=drop(Mu.hnew.updated), postvar=drop(Sigma.hnew.updated))
}

newh_postmean_index <- function(index_1, index_2, index_1_new, index_2_new, fit) {
  
  
  # if(is.null(dim(index_new))) index_new = matrix(index_new, nrow=1)
  # if(class(Znew)  == "data.frame") Znew <- data.matrix(Znew)
  # Znew            = as.matrix(Znew)
  pnew            = 2
  
  Sigsq_mean	    = mean(fit$draws('sigma_sq'))
  # Tau_mean	      = as.matrix(apply(fit$draws('tau_sq'), 2, mean))
  b_mean		      = as.matrix(apply(fit$draws('b'), 3, mean))
  Beta_mean	      = as.matrix(apply(fit$draws('Beta'), 3, mean))
  h_mean		      = as.matrix(apply(fit$draws('h'), 3, mean))
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
  eta_1_mean = as.matrix(apply(fit$draws('eta_1'), 3, mean))
  eta_2_mean = as.matrix(apply(fit$draws('eta_2'), 3, mean))
  rho_1_mean = as.matrix(apply(fit$draws('rho_1'), 3, mean))
  rho_2_mean = as.matrix(apply(fit$draws('rho_2'), 3, mean))
  Y			          = Y
  X 			        = X
  U			          = U
  W               = W

  
  list.K				  = list()
  
  # poly 				    = polydot(degree=2, offset=1)
  # poly 				    = rbfdot(sigma=1)
  
  list.K[[1]] 	= (rbf(index_1,index_1,eta_1_mean,rho_1_mean) + diag(n)*1e-7)
  list.K[[2]] = (rbf(index_2,index_2,eta_2_mean,rho_2_mean) + diag(n)*1e-7)
  
  Sigma.h 				= do.call(adiag, list.K)
  
  Sigma.h.new = "[<-"(matrix(0, (q*n+q*pnew), (q*n+q*pnew)), 1:nrow(Sigma.h), 1:ncol(Sigma.h), value = Sigma.h) # adding rows/cols of zeroes to square matrix
  
  list.G        = list()
  
  list.G[[1]] 	= (rbf(rbind(index_1, index_1_new),rbind(index_1, index_1_new),eta_1_mean,rho_1_mean) + diag(102)*1e-7)
  list.G[[2]] = (rbf(rbind(index_2, index_2_new),rbind(index_2, index_2_new),eta_2_mean,rho_2_mean) + diag(102)*1e-7)
  
  # h vector is ordered c(h1, h2, h2new, h1new)
  G.1.11 = list.G[[1]][ 1:n, 1:n ]
  G.1.22 = list.G[[1]][ (n+1):(n+pnew), (n+1):(n+pnew) ]
  G.1.21 = list.G[[1]][ (n+1):(n+pnew), 1:n ]
  G.1.12 = list.G[[1]][ 1:n, (n+1):(n+pnew) ]
  
  list.G.final = list()
  list.G.final[[1]] = G.1.11
  list.G.final[[2]] = list.G[[2]]
  list.G.final[[3]] = G.1.22
  cov.bf = do.call(adiag, list.G.final)
  
  cov.bf[(q*n + pnew + 1):(q*n + q*pnew), 1:n] = G.1.21
  cov.bf[1:n, (q*n + pnew + 1):(q*n + q*pnew)] = G.1.12
  
  Sigma11 = cov.bf[1:(n*q), 1:(n*q)]
  Sigma12 = cov.bf[1:(n*q), (n*q + 1):(n*q + pnew*q)]
  Sigma21 = cov.bf[(n*q + 1):(n*q + pnew*q), 1:(n*q)]
  Sigma22 = cov.bf[(n*q + 1):(n*q + pnew*q), (n*q + 1):(n*q + pnew*q)]
  
  WTW         = t(W) %*% W
  Sigma.h.inv = solve(Sigma.h + 1e-7*diag(nrow(Sigma.h)))
  # Sigma.h.inv = solve(Sigma.h + diag(nrow(Sigma.h)))
  m.h         = solve( WTW / Sigsq_mean + Sigma.h.inv) %*% t(W) %*% (Y - X %*% Beta_mean - U %*% b_mean) * (1/Sigsq_mean)
  s.h         = solve( WTW / Sigsq_mean + Sigma.h.inv)
  
  Sigma11.inv = Sigma.h.inv
  
  Mu.hnew.updated = Sigma21 %*% Sigma11.inv %*% m.h
  Sigma.hnew.updated = Sigma22 - Sigma21 %*% Sigma11.inv %*% Sigma12 + Sigma21 %*% Sigma11.inv %*% s.h %*% Sigma11.inv %*% Sigma12
  
  list(postmean=drop(Mu.hnew.updated), postvar=drop(Sigma.hnew.updated))
}