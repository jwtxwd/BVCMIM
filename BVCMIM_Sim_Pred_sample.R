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
newh_postmean_sample <- function(Z, Znew, fit) {


  if(is.null(dim(Znew))) Znew = matrix(Znew, nrow=1)
  # if(class(Znew)  == "data.frame") Znew <- data.matrix(Znew)
  # Znew            = as.matrix(Znew)
  pnew            = nrow(Znew)

  Sigsq_mean	    = mean(fit$draws('sigma_sq'))
  # Tau_mean	      = as.matrix(apply(fit$draws('tau_sq'), 2, mean))
  b_mean		      = as.matrix(apply(fit$draws('b'), 3, mean))
  Beta_mean	      = as.matrix(apply(fit$draws('Beta'), 3, mean))
  h_mean		      = as.matrix(apply(fit$draws('h'), 3, mean))
  weights_1_1 = as.matrix(as.data.frame(fit$draws('weights_1_1')))
  weights_1_2 = as.matrix(as.data.frame(fit$draws('weights_1_2')))
  weights_1_3 = as.matrix(as.data.frame(fit$draws('weights_1_3')))
  weights_1_4 = as.matrix(as.data.frame(fit$draws('weights_1_4')))
  weights_1_5 = as.matrix(as.data.frame(fit$draws('weights_1_5')))
  weights_2_1 = as.matrix(as.data.frame(fit$draws('weights_2_1')))
  weights_2_2 = as.matrix(as.data.frame(fit$draws('weights_2_2')))
  weights_2_3 = as.matrix(as.data.frame(fit$draws('weights_2_3')))
  weights_2_4 = as.matrix(as.data.frame(fit$draws('weights_2_4')))
  weights_2_5 = as.matrix(as.data.frame(fit$draws('weights_2_5')))
  eta_1 = as.matrix(as.data.frame(fit$draws('eta_1')))
  eta_2 = as.matrix(as.data.frame(fit$draws('eta_2')))
  rho_1 = as.matrix(as.data.frame(fit$draws('rho_1')))
  rho_2 = as.matrix(as.data.frame(fit$draws('rho_2')))
  Y			          = Y
  X 			        = X
  U			          = U
  W               = W
  
  h1_new = c()
  h2_new = c()
  for (i in 1:1000){
    index_1_1 = Z[,1:8] %*% t(t(weights_1_1[i,]))
    index_1_1_norm = (index_1_1-mean(index_1_1))/sd(index_1_1)
    index_1_2 = Z[,9:24] %*% t(t(weights_1_2[i,]))
    index_1_2_norm = (index_1_2-mean(index_1_2))/sd(index_1_2)
    index_1_3 = Z[,25:26] %*% t(t(weights_1_3[i,]))
    index_1_3_norm = (index_1_3-mean(index_1_3))/sd(index_1_3)
    index_1_4 = Z[,27:32] %*% t(t(weights_1_4[i,]))
    index_1_4_norm = (index_1_4-mean(index_1_4))/sd(index_1_4)
    index_1_5 = Z[,33:36] %*% t(t(weights_1_5[i,]))
    index_1_5_norm = (index_1_5-mean(index_1_5))/sd(index_1_5)
    index_1 = cbind(index_1_1_norm, index_1_2_norm, index_1_3_norm, index_1_4_norm, index_1_5_norm)
    
    index_2_1 = Z[,1:8] %*% t(t(weights_2_1[i,]))
    index_2_1_norm = (index_2_1-mean(index_2_1))/sd(index_2_1)
    index_2_2 = Z[,9:24] %*% t(t(weights_2_2[i,]))
    index_2_2_norm = (index_2_2-mean(index_2_2))/sd(index_2_2)
    index_2_3 = Z[,25:26] %*% t(t(weights_2_3[i,]))
    index_2_3_norm = (index_2_3-mean(index_2_3))/sd(index_2_3)
    index_2_4 = Z[,27:32] %*% t(t(weights_2_4[i,]))
    index_2_4_norm = (index_2_4-mean(index_2_4))/sd(index_2_4)
    index_2_5 = Z[,33:36] %*% t(t(weights_2_5[i,]))
    index_2_5_norm = (index_2_5-mean(index_2_5))/sd(index_2_5)
    index_2 = cbind(index_2_1_norm, index_2_2_norm, index_2_3_norm, index_2_4_norm, index_2_5_norm)
    
    index_1_1_new = Znew[,1:8] %*% t(t(weights_1_1[i,]))
    index_1_1_new_norm = (index_1_1_new-mean(index_1_1))/sd(index_1_1)
    index_1_2_new = Znew[,9:24] %*% t(t(weights_1_2[i,]))
    index_1_2_new_norm = (index_1_2_new-mean(index_1_2))/sd(index_1_2)
    index_1_3_new = Znew[,25:26] %*% t(t(weights_1_3[i,]))
    index_1_3_new_norm = (index_1_3_new-mean(index_1_3))/sd(index_1_3)
    index_1_4_new = Znew[,27:32] %*% t(t(weights_1_4[i,]))
    index_1_4_new_norm = (index_1_4_new-mean(index_1_4))/sd(index_1_4)
    index_1_5_new = Znew[,33:36] %*% t(t(weights_1_5[i,]))
    index_1_5_new_norm = (index_1_5_new-mean(index_1_5))/sd(index_1_5)
    index_1_new = cbind(index_1_1_new_norm, index_1_2_new_norm, index_1_3_new_norm, index_1_4_new_norm, index_1_5_new_norm)
    
    index_2_1_new = Znew[,1:8] %*% t(t(weights_2_1[i,]))
    index_2_1_new_norm = (index_2_1_new-mean(index_2_1))/sd(index_2_1)
    index_2_2_new = Znew[,9:24] %*% t(t(weights_2_2[i,]))
    index_2_2_new_norm = (index_2_2_new-mean(index_2_2))/sd(index_2_2)
    index_2_3_new = Znew[,25:26] %*% t(t(weights_2_3[i,]))
    index_2_3_new_norm = (index_2_3_new-mean(index_2_3))/sd(index_2_3)
    index_2_4_new = Znew[,27:32] %*% t(t(weights_2_4[i,]))
    index_2_4_new_norm = (index_2_4_new-mean(index_2_4))/sd(index_2_4)
    index_2_5_new = Znew[,33:36] %*% t(t(weights_2_5[i,]))
    index_2_5_new_norm = (index_2_5_new-mean(index_2_5))/sd(index_2_5)
    index_2_new = cbind(index_2_1_new_norm, index_2_2_new_norm, index_2_3_new_norm, index_2_4_new_norm, index_2_5_new_norm)
    
    list.K				  = list()
    list.K[[1]] 	= (rbf(rbind(index_1, index_1_new),rbind(index_1, index_1_new),eta_1[i],rho_1[i,]) + diag(n+dim(Znew)[1])*1e-7)
    list.K[[2]] = (rbf(rbind(index_2, index_2_new),rbind(index_2, index_2_new),eta_2[i],rho_2[i,]) + diag(n+dim(Znew)[1])*1e-7)
    
    h1 = mvrnorm(mu=rep(0, n+dim(Znew)[1]), Sigma = list.K[[1]])
    h2 = mvrnorm(mu=rep(0, n+dim(Znew)[1]), Sigma = list.K[[2]])
    
    h1_new = cbind(h1_new, h1)
    h2_new = cbind(h2_new, h2)
    
    # print(i)
    
  }
  list(postmean=c(mean(h1_new[101,]), mean(h1_new[102,])), postvar=c(sd(h1_new[101,]), sd(h2_new[102,])), h1_new=h1_new, h2_new=h2_new)
  # list(h1_new=h1_new, h2_new=h2_new)
}

