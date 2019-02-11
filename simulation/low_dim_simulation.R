rm(list=ls())
library(simulation)
library(cpReg)
source("../simulation/SGL_solver.R")

paramMat <- as.matrix(expand.grid(round(exp(seq(log(20), log(1000), length.out = 10))), c(1,2)))
colnames(paramMat) <- c("n", "X_type")

##############

edge_mutation <- function(lev, n=200){
  mn = rep(0,n)
  mn[seq(from=n-40+1, to=n)] = lev
  return(mn)
}

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}
true_jumps <- c(100, 140)

#########

test_func_closure <- function(contrast){
  function(y, fit = NA, jump = NA){
    as.numeric(contrast %*% y)
  }
}
declutter_func <- function(x){selectiveModel::declutter(x, sign_vec = rep(1, length(x)),
                                                        how_close = 2)$jump_vec}

#############

rule <- function(vec){
  middle_mutation(lev = vec["SnR"], n = n) + stats::rnorm(n)
}

criterion <- function(dat, vec, y){
  if(vec["method"] == 1){
    fit_method <- function(x){binseginf::bsfs(x, numSteps = vec["ksteps"])}
  } else {
    fit_method <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = vec["ksteps"])}
  }
  fit <- fit_method(dat)

  sign_mat <- binseginf::jump_sign(fit)
  if(vec["decluttered"] == 0){
    tmp <- unlist(lapply(true_jumps, function(x){x + c(-2:2)}))
    cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                              how_close = 0,
                                              desired_jumps = tmp)
  } else {
    cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                              how_close = 2,
                                              desired_jumps = true_jumps)
  }

  res <- rep(NA, 3*vec["ksteps"])
  len <- length(cluster_list$jump_vec)
  res[1:len] <- cluster_list$jump_vec
  names(res) <- c(paste0("Jump ", 1:vec["ksteps"]), paste0("Direction ", 1:vec["ksteps"]),
                  paste0("Pvalue ", 1:vec["ksteps"]))

  for(i in 1:len){
    if(cluster_list$target_bool[i]){
      set.seed(10*y)
      contrast <- selectiveModel:::contrast_from_cluster(cluster_list, n, i)
      test_func <- test_func_closure(contrast)
      if(cluster_list$sign_mat["sign:-1",i] == 0){
        direction <- 1
      } else if(cluster_list$sign_mat["sign:+1",i] == 0){
        direction <- -1
      } else {
        direction <- NA
      }

      tmp <- selectiveModel::selected_model_inference(dat, fit_method = fit_method,
                                                      test_func = test_func,
                                                      declutter_func = declutter_func,
                                                      num_samp = num_samp,
                                                      direction = direction,
                                                      ignore_jump = i, sigma = vec["sigma"],
                                                      verbose = F, param = list(burn_in = burn_in,
                                                                                lapse = 1))
      res[i+vec["ksteps"]] <- direction
      res[i+2*vec["ksteps"]] <- tmp$pval
    }
  }
  res
}

# set.seed(1); criterion(rule(paramMat[6,]), paramMat[6,], 1)
# set.seed(1); criterion(rule(paramMat[12,]), paramMat[12,], 1)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = paramMat[,"trials"],
                                        cores = 15, as_list = F,
                                        filepath = "main_powercurve_onejump_fl_tmp.RData",
                                        verbose = T)
save.image("main_powercurve_onejump_fl.RData")
