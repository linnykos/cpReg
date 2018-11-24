rm(list=ls())
load("../results/SEPP_reproduce_simulation.RData")

# compute the MSE
mse_mat <- sapply(res, function(x){
  sapply(x, function(y){
    sum((y$A - A)^2)
  })
})
mse_vec <- apply(mse_mat, 2, mean)

l0_mat <- sapply(res, function(x){
  sapply(x, function(y){
    idx1 <- which(y$A != 0)
    idx2 <- which(A != 0)
    length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
  })
})
l0_vec <- apply(l0_mat, 2, mean)

lambda_mat <- sapply(res, function(x){
  sapply(x, function(y){
    y$lambda
  })
})
lambda_vec <- apply(lambda_mat, 2, mean)

##########

# plotting

png("../figure/SEPP_reproduce.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3))
plot(paramMat[,1], mse_vec, pch = 16, xlab = "T", ylab = "MSE")
plot(paramMat[,1], l0_vec, pch = 16, xlab = "T", ylab = "Jaccard Index")
plot(paramMat[,1], lambda_vec, pch = 16, xlab = "T", ylab = "Lambda")
graphics.off()

