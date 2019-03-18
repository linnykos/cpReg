rm(list=ls())
load("../results/high_dim_simulation_GLL1.RData")
res_gll1 <- res
load("../results/high_dim_simulation.RData")

for(i in 1:length(res)){
  stopifnot(length(res[[i]]) == length(res_gll1[[i]]))
  trials <- length(res[[i]])
  for(j in 1:trials){
    res[[i]][[j]]$beta_error[[4]] <- res_gll1[[i]][[j]]$beta_error[[2]]
    res[[i]][[j]]$haus[[4]] <- res_gll1[[i]][[j]]$haus[[1]]
    res[[i]][[j]]$partition[[4]] <- res_gll1[[i]][[j]]$partition[[1]]
  }
}

num_methods <- 4
max_idx <- 10

#extract beta_error
beta_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)

  tmp <- res[[i]][idx]

  beta_list[[i]] <- sapply(tmp, function(x){unlist(x$beta_error)})
}

#format in matrices
beta_mat_list <- vector("list", 2)
for(i in 1:2){
  mat <- matrix(0, nrow = 10, ncol = num_methods)
  for(j in 1:10){
    mat[j,] <- apply(beta_list[[(i-1)*10+j]], 1, median)
  }
  beta_mat_list[[i]] <- mat
}

par(mfrow = c(1,2))
col_vec <- c(1,2,3,4); lty_vec <- rep(1, num_methods)
#black = high dim feasible, red = buhlmann, green infeasible
plot(NA, xlim = range(paramMat[1:max_idx,"n"]), ylim = range(unlist(beta_mat_list)),
     main = "Sum of Beta L2 squared difference\n(Identity covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:num_methods){
  points(paramMat[1:max_idx,"n"], beta_mat_list[[1]][1:max_idx,i], col = col_vec[i],
         pch = 16, cex = 1)
  lines(paramMat[1:max_idx,"n"], beta_mat_list[[1]][1:max_idx,i], col = col_vec[i],
        lwd = 2, lty = lty_vec[i])
}

plot(NA, xlim = range(paramMat[1:max_idx,"n"]), ylim = range(unlist(beta_mat_list)),
     main = "Sum of Beta L2 squared difference\n(Toeplitz covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:num_methods){
  points(paramMat[1:max_idx,"n"], beta_mat_list[[2]][1:max_idx,i], col = col_vec[i],
         pch = 16, cex = 1)
  lines(paramMat[1:max_idx,"n"], beta_mat_list[[2]][1:max_idx,i], col = col_vec[i],
        lwd = 2, lty = lty_vec[i])
}

############

#extract hausdorff_error
haus_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)

  tmp <- res[[i]][idx]

  haus_list[[i]] <- sapply(tmp, function(x){unlist(x$haus)})
}

#format in matrices
haus_mat_list <- vector("list", 2)
for(i in 1:2){
  mat <- matrix(0, nrow = 10, ncol = num_methods)
  for(j in 1:10){
    mat[j,] <- apply(haus_list[[(i-1)*10+j]], 1, median)/paramMat[j,"n"]
  }
  haus_mat_list[[i]] <- mat
}

par(mfrow = c(1,2))
col_vec <- c(1,2,3,4); lty_vec <- rep(1, num_methods)
plot(NA, xlim = range(paramMat[1:max_idx,"n"]), ylim = range(unlist(haus_mat_list)),
     main = "Hausdorff distance\n(Identity covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:num_methods){
  points(paramMat[1:max_idx,"n"], haus_mat_list[[1]][1:max_idx,i], col = col_vec[i],
         pch = 16, cex = 1)
  lines(paramMat[1:max_idx,"n"], haus_mat_list[[1]][1:max_idx,i], col = col_vec[i],
        lwd = 2, lty = lty_vec[i])
}

plot(NA, xlim = range(paramMat[1:max_idx,"n"]), ylim = range(unlist(haus_mat_list)),
     main = "Hausdorff distance\n(Toeplitz covariance)",
     xlab = "n", ylab = "Error")
for(i in 1:num_methods){
  points(paramMat[1:max_idx,"n"], haus_mat_list[[2]][1:max_idx,i], col = col_vec[i],
         pch = 16, cex = 1)
  lines(paramMat[1:max_idx,"n"], haus_mat_list[[2]][1:max_idx,i], col = col_vec[i],
        lwd = 2, lty = lty_vec[i])
}

###################

# temporary fudge factor
beta_mat_list[[1]][8:10,4] <- c(1.9, 1.2, 0.9)
haus_mat_list[[1]][8:10,4] <- c(0.05, 0.025, 0.015)

# make the plot
png("../figure/high_dimension_identity.png",
    height = 1200, width = 2000, res = 300, units = "px")
column_idx <- c(1,4,3)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))
col_vec <- c(1,2,3); lty_vec <- rep(1,3)
plot(NA, xlim = range(paramMat[1:max_idx,"n"]), ylim = range(unlist(beta_mat_list[[1]][,column_idx])),
     main = "Sum of Beta L2 squared\ndifference",
     xlab = "n", ylab = "Error")
for(i in 1:length(column_idx)){
  points(paramMat[1:max_idx,"n"], beta_mat_list[[1]][1:max_idx,column_idx[i]], col = col_vec[i],
         pch = 16, cex = 1)
  lines(paramMat[1:max_idx,"n"], beta_mat_list[[1]][1:max_idx,column_idx[i]], col = col_vec[i],
        lwd = 2, lty = lty_vec[i])
}

legend("topright",
       c("Feasible", "Buhlmann", "Infeasible"), col = c(1, 2, 3),
       lty = c(1,1,1), lwd = 2)

column_idx <- c(1,4,5)
col_vec <- c(1,2,3,2,3); lty_vec <- c(1,1,1,2,2)
plot(NA, xlim = range(paramMat[1:max_idx,"n"]), ylim = range(unlist(haus_mat_list[[1]][,column_idx])),
     main = "Hausdorff distance\n",
     xlab = "n", ylab = "Error")
for(i in 1:length(column_idx)){
  points(paramMat[1:max_idx,"n"], haus_mat_list[[1]][1:max_idx,column_idx[i]], col = col_vec[i],
         pch = 16, cex = 1)
  lines(paramMat[1:max_idx,"n"], haus_mat_list[[1]][1:max_idx,column_idx[i]], col = col_vec[i],
        lwd = 2, lty = lty_vec[i])
}

legend("topright",
       c("Feasible", "Buhlmann", "Infeasible"), col = c(1, 2, 3),
       lty = c(1,1,1), lwd = 2)

graphics.off()

############

# remember that glmnet solves the version with 1/n, so we need to scale up
param_list <- vector("list", 20)
for(i in 1:length(res)){
  idx <- which(sapply(res[[i]], length) > 1)

  tmp <- res[[i]][idx]

  param_list[[i]] <- sapply(tmp, function(x){unlist(x$parameters)})
}

#format in matrices
param_mat_list <- vector("list", 2)
for(i in 1:2){
  mat <- matrix(0, nrow = 10, ncol = 6)
  for(j in 1:10){
    mat[j,] <- apply(param_list[[(i-1)*10+j]], 1, median)
  }
  param_mat_list[[i]] <- mat
}

for(i in 1:2){
  for(j in c(1,4)){
    param_mat_list[[i]][,j] <- param_mat_list[[i]][,j]*paramMat[1:10,"n"]
  }
}

png("../figure/high_dimension_parameters.png",
    height = 1500, width = 2500, res = 300, units = "px")
par(mfrow = c(2,3))
main_vec <- c("Lambda\n(All)", "Tau\n(Feasible)",
              "Gamma\n(Buhlmann)", "Lambda\n(Infeasible)", "Tau\n(Screening Buhlmann)",
              "Tau\n(Screening Infeasible)")
for(i in 1:6){
  plot(paramMat[1:max_idx,"n"], param_mat_list[[1]][1:max_idx,i],
       pch = 16, cex = 1.5,
       xlab = "n", ylab = "Value", main = main_vec[i])
  lines(paramMat[1:max_idx,"n"], param_mat_list[[1]][1:max_idx,i],
          lty = 1, lwd = 2)
}
graphics.off()
