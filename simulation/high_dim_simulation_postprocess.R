rm(list=ls())
load("../results/high_dim_simulation_feasible.RData")
res_feasible <- res
load("../results/high_dim_simulation_infeasible.RData")
res_infeasible <- res

#extract beta_error
beta_list <- vector("list", 20)
for(i in 1:length(res)){
  idx1 <- which(sapply(res_feasible[[i]], length) > 1)
  idx2 <- which(sapply(res_infeasible[[i]], length) > 1)

  idx <- intersect(idx1, idx2)

  tmp1 <- res_feasible[[i]][idx]
  tmp2 <- res_infeasible[[i]][idx]

  beta_list[[i]] <- sapply(1:length(tmp1), function(x){
    c(tmp1[[x]]$beta_error1, tmp1[[x]]$beta_error2,
      tmp2[[x]]$beta_error3)
  })
}

#format in matrices
beta_mat_list <- vector("list", 2)
for(i in 1:2){
  mat <- matrix(0, nrow = 10, ncol = 3)
  for(j in 1:10){
    mat[j,] <- apply(beta_list[[(i-1)*10+j]], 1, median)
  }
  beta_mat_list[[i]] <- mat
}

par(mfrow = c(1,2))
col_vec <- c(1,2,3)
plot(NA, xlim = range(paramMat[,"n"]), ylim = range(unlist(beta_mat_list)),
     main = "Sum of Beta L2 squared difference\n(Identity covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:3){
  points(paramMat[1:10,"n"], beta_mat_list[[1]][,i], col = col_vec[i], pch = 16, cex = 1)
  lines(paramMat[1:10,"n"], beta_mat_list[[1]][,i], col = col_vec[i], lwd = 2)
}

plot(NA, xlim = range(paramMat[,"n"]), ylim = range(unlist(beta_mat_list)),
     main = "Sum of Beta L2 squared difference\n(Toeplitz covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:3){
  points(paramMat[1:10,"n"], beta_mat_list[[2]][,i], col = col_vec[i], pch = 16, cex = 1)
  lines(paramMat[1:10,"n"], beta_mat_list[[2]][,i], col = col_vec[i], lwd = 2)
}


############

#extract hausdorff_error
haus_list <- vector("list", 20)
for(i in 1:length(res)){
  idx1 <- which(sapply(res_feasible[[i]], length) > 1)
  idx2 <- which(sapply(res_infeasible[[i]], length) > 1)

  idx <- intersect(idx1, idx2)

  tmp1 <- res_feasible[[i]][idx]
  tmp2 <- res_infeasible[[i]][idx]

  haus_list[[i]] <- sapply(1:length(tmp1), function(x){
    c(tmp1[[x]]$haus1, tmp1[[x]]$haus2,
      tmp2[[x]]$haus3)
  })
}

#format in matrices
haus_mat_list <- vector("list", 2)
for(i in 1:2){
  mat <- matrix(0, nrow = 10, ncol = 3)
  for(j in 1:10){
    mat[j,] <- apply(haus_list[[(i-1)*10+j]], 1, median)/paramMat[j,"n"]
  }
  haus_mat_list[[i]] <- mat
}

par(mfrow = c(1,2))
col_vec <- c(1,2,3)
plot(NA, xlim = range(paramMat[,"n"]), ylim = range(unlist(haus_mat_list)),
     main = "Hausdorff distance\n(Identity covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:3){
  points(paramMat[1:10,"n"], haus_mat_list[[1]][,i], col = col_vec[i], pch = 16, cex = 1)
  lines(paramMat[1:10,"n"], haus_mat_list[[1]][,i], col = col_vec[i], lwd = 2)
}

plot(NA, xlim = range(paramMat[,"n"]), ylim = range(unlist(haus_mat_list)),
     main = "Hausdorff distance\n(Toeplitz covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:3){
  points(paramMat[1:10,"n"], haus_mat_list[[2]][,i], col = col_vec[i], pch = 16, cex = 1)
  lines(paramMat[1:10,"n"], haus_mat_list[[2]][,i], col = col_vec[i], lwd = 2)
}

############

# NEED TO UPDATE

# extract lambda
# remember that glmnet solves the version with 1/n, so we need to scale up
lambda_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)
  tmp <- res[[i]][idx]
  lambda_list[[i]] <- sapply(tmp, function(x){x$lambda})
}
lambda_mat <- matrix(sapply(lambda_list, median), ncol = 2)
lambda_mat <- apply(lambda_mat, 2, function(x){x*paramMat[1:10,"n"]})
plot(NA, xlim = range(paramMat[,"n"]), ylim = range(lambda_mat),
     main = "Oracle Lambda",
     xlab = "n", ylab = "Lambda value")
points(paramMat[1:10,"n"], lambda_mat[,1], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], lambda_mat[,1], col = "black", lwd = 2)
points(paramMat[1:10,"n"], lambda_mat[,2], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], lambda_mat[,2], col = "black", lwd = 2, lty = 2)

####################

# extract tau
tau_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)
  tmp <- res[[i]][idx]
  tau_list[[i]] <- sapply(tmp, function(x){x$tau})
}
tau_mat <- matrix(sapply(tau_list, median), ncol = 2)
plot(NA, xlim = range(paramMat[,"n"]), ylim = range(tau_mat),
     main = "Oracle tau",
     xlab = "n", ylab = "tau value")
points(paramMat[1:10,"n"], tau_mat[,1], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], tau_mat[,1], col = "black", lwd = 2)
points(paramMat[1:10,"n"], tau_mat[,2], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], tau_mat[,2], col = "black", lwd = 2, lty = 2)


