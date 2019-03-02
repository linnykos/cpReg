rm(list=ls())
load("../results/high_dim_simulation_1000.RData")

#extract beta_error
beta_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)
  tmp <- res[[i]][idx]
  beta_list[[i]] <- sapply(tmp, function(x){x$beta_error})
}
beta_mat <- matrix(sapply(beta_list, median), ncol = 2)
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


