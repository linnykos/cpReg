high_dim_feasible_estimate <- function(X, y, lambda, tau, M = 100,
                                       delta = 10){
  n <- length(y); tree <- .create_node(1, n)
  cp <- c()
  steps <- 1

  while(TRUE){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])

      if(is.na(leaf$breakpoint) & leaf$end - leaf$start > 2*delta){
        res <- .find_breakpoint2(X, y, lambda, leaf$start, leaf$end,
                                 delta)

        leaf$breakpoint <- res$breakpoint; leaf$cusum <- res$cusum
      }
    }

    node.name <- .find_leadingBreakpoint(tree)
    if(length(node.name) == 0) break()
    node.selected <- data.tree::FindNode(tree, node.name)

    if(node.selected$cusum < tau) break()

    node.selected$active <- steps
    steps <- steps + 1
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
  }

  structure(list(tree = tree, numSteps = steps), class = "HDF")
}

#########

#' Checks S3 object validity
#'
#' Generic function that fails noisily with a stop message if the object is
#' invalid. Otherwise, nothing happens.
#'
#' @param  obj  The object to check
#' @return void
#' @export
is_valid <- function(obj) UseMethod("is_valid")

#' Check whether tree is valid or not
#'
#' @param obj The tree (of class Node)
#'
#' @return TRUE if valid
#' @export
is_valid.Node <- function(obj){
  if(obj$start > obj$end) stop("the start must be less or equal to end")
  if(!is.na(obj$breakpoint) & (obj$start > obj$breakpoint & obj$end < obj$breakpoint))
    stop("breakpoint must be between start and end (inclusive)")

  TRUE
}

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


.find_breakpoint2 <- function(X, y, lambda, start, end, delta){
  if(start > end) stop("start must be smaller than or equal to end")
  if(start == end) return(list(breakpoint = start, cusum = 0))
  stopifnot(end - start > 2*delta)

  breakpoint <- seq(from = start + delta, to = end - 1 - delta, by = 1)
  cusum_vec <- sapply(breakpoint, function(x){
    .regression_cusum(X = X, y = y, start = start, end = end,
                      lambda = lambda, breakpoint = x)})

  cusum_vec <- apply(cusum_vec, 2, .l2norm)

  idx <- which.max(cusum_vec)
  list(breakpoint = breakpoint[idx], cusum = cusum_vec[idx])
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
  as.numeric(glmnet::coef.glmnet(glmnet_res, s = lambda))[-1]
}

