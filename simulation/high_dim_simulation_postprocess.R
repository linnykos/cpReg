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

plot(NA, xlim = range(paramMat[,"n"]), ylim = range(beta_mat),
     main = "Sum of Beta L2 squared difference",
     xlab = "n", ylab = "Error")
points(paramMat[1:10,"n"], beta_mat[,1], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], beta_mat[,1], col = "black", lwd = 2)
points(paramMat[1:10,"n"], beta_mat[,2], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], beta_mat[,2], col = "black", lwd = 2, lty = 2)

############

#extract hausdorff_error
haus_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)
  tmp <- res[[i]][idx]
  haus_list[[i]] <- sapply(tmp, function(x){x$haus})
}
haus_mat <- matrix(sapply(haus_list, median), ncol = 2)
haus_mat <- apply(haus_mat, 2, function(x){x/paramMat[1:10,"n"]})
plot(NA, xlim = range(paramMat[,"n"]), ylim = range(haus_mat),
     main = "Hausdorff distance divided by n",
     xlab = "n", ylab = "Error")
points(paramMat[1:10,"n"], haus_mat[,1], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], haus_mat[,1], col = "black", lwd = 2)
points(paramMat[1:10,"n"], haus_mat[,2], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], haus_mat[,2], col = "black", lwd = 2, lty = 2)

############

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


