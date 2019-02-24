screening <- function(fit, tau, M = 100, delta = 1, verbose = F){
  wbs(fit, data_length_func = nrow,
      compute_cusum_func = .compute_cusum,
      tau = tau, M = M, delta = delta, verbose = verbose)
}

hausdorff <- function(set1, set2, one.sided = FALSE){
  #handle corner cases
  if(length(set1) == 0 | length(set2) == 0) return(NA)
  if(length(set2) == 1) set2 = c(set2, set2)
  if(length(set1) == 1) set1 = c(set1, set1)

  dist.mat = sapply(set1, function(i){abs(i-set2)})
  if(class(dist.mat) != "matrix") dist.mat = as.matrix(dist.mat, nrow = 1)
  dist.vecx = apply(dist.mat, 2, min)

  if(!one.sided) dist.vecy = apply(dist.mat, 1, min) else dist.vecy = 0

  max(dist.vecx, dist.vecy)
}

################

.compute_cusum <- function(fit, s, e, b){
  stopifnot(s < b, b < e, s >= 0, e <= nrow(fit))

  factor1 <- sqrt((e-b)/((e-s)*(b-s)))
  factor2 <- sqrt((b-s)/((e-s)*(e-b)))
  .l2norm(factor1*colSums(fit[(s+1):b,,drop = F]) - factor2*colSums(fit[(b+1):e,,drop = F]))
}

.l2norm <- function(x){
  sqrt(sum(x^2))
}
