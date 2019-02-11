rm(list=ls())
load("low_dim_simulation.RData")

# gather the results
res_mat <- vector("list", 20)
for(i in 1:20){
  idx <- which(sapply(res[[i]], length) > 1)
  res_mat[[i]] <- matrix(unlist(res[[i]][idx]), nrow = 4)
}

# construct the vectors
beta_error <- matrix(0, nrow = 4, ncol = 10)
for(i in 1:10){
  beta_error[1,i] <- median(res_mat[[i]][1,], na.rm = T)
  beta_error[2,i] <- median(res_mat[[i]][2,], na.rm = T)

  beta_error[3,i] <- median(res_mat[[i+10]][1,], na.rm = T)
  beta_error[4,i] <- median(res_mat[[i+10]][2,], na.rm = T)
}

hausdorff_error <- matrix(0, nrow = 4, ncol = 10)
for(i in 1:10){
  hausdorff_error[1,i] <- median(res_mat[[i]][3,], na.rm = T)
  hausdorff_error[2,i] <- median(res_mat[[i]][4,], na.rm = T)

  hausdorff_error[3,i] <- median(res_mat[[i+10]][3,], na.rm = T)
  hausdorff_error[4,i] <- median(res_mat[[i+10]][4,], na.rm = T)
}

png("../figure/low_dimensional_beta.png",
    height = 1500, width = 2000, res = 300, units = "px")
plot(NA, xlim = range(paramMat[1:10,"n"]), ylim = range(beta_error),
     main = "Sum of Beta L2 squared difference",
     xlab = "n", ylab = "Error")
points(paramMat[1:10,"n"], beta_error[1,], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], beta_error[1,], col = "black", lwd = 2)
points(paramMat[1:10,"n"], beta_error[2,], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], beta_error[2,], col = "black", lwd = 2, lty = 2)
points(paramMat[1:10,"n"], beta_error[3,], col = "red", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], beta_error[3,], col = "red", lwd = 2)
points(paramMat[1:10,"n"], beta_error[4,], col = "red", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], beta_error[4,], col = "red", lwd = 2, lty = 2)
legend("topright", c("Potts - Identity", "Potts - Toeplitz", "SGL - Identity",
                    "SGL - Toeplitz"), cex = 0.8, col = c("black", "red", "black", "red"),
       lty = c(1,1,2,2))
graphics.off()

png("../figure/low_dimensional_hausdorff.png",
    height = 1500, width = 2000, res = 300, units = "px")
plot(NA, xlim = range(paramMat[1:10,"n"]), ylim = range(hausdorff_error),
     main = "Hausdorff distance", xlab = "n", ylab = "Error")
points(paramMat[1:10,"n"], hausdorff_error[1,], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], hausdorff_error[1,], col = "black", lwd = 2)
points(paramMat[1:10,"n"], hausdorff_error[2,], col = "black", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], hausdorff_error[2,], col = "black", lwd = 2, lty = 2)
points(paramMat[1:10,"n"], hausdorff_error[3,], col = "red", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], hausdorff_error[3,], col = "red", lwd = 2)
points(paramMat[1:10,"n"], hausdorff_error[4,], col = "red", pch = 16, cex = 1)
lines(paramMat[1:10,"n"], hausdorff_error[4,], col = "red", lwd = 2, lty = 2)
legend("topleft", c("Potts - Identity", "Potts - Toeplitz", "SGL - Identity",
                        "SGL - Toeplitz"), cex = 0.8, col = c("black", "red", "black", "red"),
       lty = c(1,1,2,2))
graphics.off()
