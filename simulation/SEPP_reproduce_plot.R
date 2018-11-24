rm(list=ls())
load("../results/SEPP_reproduce_simulation.RData")

# compute the MSE
mse_mat1 <- sapply(res, function(x){
  sapply(x, function(y){
    sum((y$fit1$A - A)^2)
  })
})
mse_vec1 <- apply(mse_mat1, 2, mean)

mse_mat2 <- sapply(res, function(x){
  sapply(x, function(y){
    sum((y$fit2$A - A)^2)
  })
})
mse_vec2 <- apply(mse_mat2, 2, mean)

##

l0_mat1 <- sapply(res, function(x){
  sapply(x, function(y){
    idx1 <- which(y$fit1$A != 0)
    idx2 <- which(A != 0)
    length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
  })
})
l0_vec1 <- apply(l0_mat1, 2, mean)

l0_mat2 <- sapply(res, function(x){
  sapply(x, function(y){
    idx1 <- which(y$fit2$A != 0)
    idx2 <- which(A != 0)
    length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
  })
})
l0_vec2 <- apply(l0_mat2, 2, mean)

##

lambda_mat1 <- sapply(res, function(x){
  sapply(x, function(y){
    y$fit1$lambda
  })
})
lambda_vec1 <- apply(lambda_mat1, 2, mean)

lambda_mat2 <- sapply(res, function(x){
  sapply(x, function(y){
    y$fit2$lambda
  })
})
lambda_vec2 <- apply(lambda_mat2, 2, mean)

##########

# plotting

png("../figure/SEPP_reproduce.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3))
plot(paramMat[,1], mse_vec1, pch = 16, xlab = "T", ylab = "MSE",
     ylim = range(c(mse_vec1, mse_vec2)))
points(paramMat[,1], mse_vec2, pch = 16, col = "red")

plot(paramMat[,1], l0_vec1, pch = 16, xlab = "T", ylab = "Jaccard Index",
     ylim = range(c(l0_vec1, l0_vec2)))
points(paramMat[,1], l0_vec2, pch = 16, col = "red")

plot(paramMat[,1], lambda_vec1, pch = 16, xlab = "T", ylab = "Lambda",
     ylim = range(c(lambda_vec1, lambda_vec2)))
points(paramMat[,1], lambda_vec2, pch = 16, col = "red")

graphics.off()

