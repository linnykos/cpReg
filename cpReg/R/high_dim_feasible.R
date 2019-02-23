high_dim_feasible_estimate <- function(X, y, lambda, tau, M = 100){
  n <- length(y); tree <- .create_node(1, n)
  cp <- c()

  for(steps in 1:numSteps){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])

      res <- .find_breakpoint(y, leaf$start, leaf$end)

      leaf$breakpoint <- res$breakpoint; leaf$cusum <- res$cusum
    }

    node.name <- .find_leadingBreakpoint(tree)
    node.selected <- data.tree::FindNode(tree, node.name)
    node.selected$active <- steps
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
  }

  obj <- structure(list(tree = tree, y.fit = y.fit, numSteps = numSteps), class = "bsFs")
  cp <- jumps(obj)
  leaves <- .enumerate_splits(tree)
  cp.sign <- sign(as.numeric(sapply(leaves, function(x){
    data.tree::FindNode(tree, x)$cusum})))
  obj <- structure(list(tree = tree, y.fit = y.fit, numSteps = numSteps, cp = cp,
                        cp.sign=cp.sign, y=y), class = "bsFs")

}

#########

.create_node <- function(start, end, breakpoint = NA, cusum = NA, active = NA){
  node <- data.tree::Node$new(paste0(start, "-", end))

  node$start <- start
  node$end <- end
  node$breakpoint <- breakpoint
  node$cusum <- cusum
  node$active <- active

  is_valid(node)

  node
}

.get_leaves_names <- function(tree){
  leaves <- tree$leaves
  vec <- sapply(leaves, function(x){x$name})
  names(vec) <- NULL
  sort(vec)
}


.find_breakpoint <- function(X, y, lambda, start, end){
  ## stopifnot(!any(duplicated(y)))
  if(start > end) stop("start must be smaller than or equal to end")
  if(start == end) return(list(breakpoint = start, cusum = 0))

  breakpoint <- seq(from = start, to = end - 1, by = 1)
  cusum.vec <- sapply(breakpoint, .regression_cusum, X = X, y = y,
                      start = start, end = end,
                      lambda = lambda)

  idx <- which.max(apply(cusum.vec, 2, .l2norm))
  list(breakpoint = breakpoint[idx], cusum = cusum.vec[idx])
}

.find_leadingBreakpoint <- function(tree){
  leaves.names <- .get_leaves_names(tree)
  cusum.vec <- sapply(leaves.names, function(x){
    data.tree::FindNode(tree, x)$cusum
  })

  leaves.names[which.max(abs(cusum.vec))]
}

.split_node <- function(node){
  if(is.na(node$breakpoint)) stop("node does not have a set breakpoint yet")
  if(node$breakpoint >= node$end) stop("node breakpoint must be less than end")

  left <- .create_node(node$start, node$breakpoint)
  right <- .create_node(node$breakpoint + 1, node$end)

  list(left = left, right = right)
}

.enumerate_splits <- function(tree){
  names(sort(tree$Get("active")))
}

.regression_cusum <- function(X, y, start, end, lambda, breakpoint){
  X1 <- X[start:breakpoint,,drop = F]
  y1 <- y[start:breakpoint]
  X2 <- X[(breakpoint+1):end,,drop = F]
  y2 <- y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, lambda)
  beta2 <- .lasso_regression(X2, y2, lambda)

  sqrt((breakpoint - start)*(end-breakpoint)/(end-start))*(beta1 - beta2)
}

.lasso_regression <- function(X, y, lambda){
  glmnet_res <- glmnet::glmnet(X, y, intercept = F)
  as.numeric(glmnet::cv.glmnet(glmnet_res, s = lambda))[-1]
}

