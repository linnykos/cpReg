#' High dimensional estimator - Buhlmann
#'
#' @param X \code{n} by \code{d} matrix
#' @param y length \code{n} vector
#' @param lambda numeric
#' @param k numeric
#' @param delta numeric
#' @param max_candidates numeric
#' @param max_changepoints numeric
#' @param verbose boolean
#'
#' @return list containing \code{partition} and \code{coef_list}
#' @export
high_dim_buhlmann_estimate <- function(X, y, lambda, k,
                                       delta = 10, max_candidates = NA,
                                       max_changepoints = NA,
                                       verbose = F){
  stopifnot(nrow(X) == length(y))
  data <- list(X = X, y = y)

  #initialization
  n <- length(y); tree <- .create_node(0, n)
  cp <- c()

  for(steps in 1:k){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])

      res <- .find_breakpoint(data, c(leaf$start, leaf$end), delta = delta,
                              max_candidates = max_candidates,
                              data_length_func = function(x){nrow(x$X)},
                              compute_cusum_func = .compute_regression_buhlmann,
                              verbose = verbose, lambda = lambda)

      leaf$breakpoint <- res$b; leaf$val <- res$val
    }

    node.name <- .find_leadingBreakpoint(tree)
    node.selected <- data.tree::FindNode(tree, node.name)
    node.selected$active <- steps
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
  }

  partition <- c(0, .jumps(tree), n)
  coef_list <- .refit_high_dim(X, y, lambda, partition/n)

  list(partition = partition, coef_list = coef_list)
}

##################

.compute_regression_buhlmann <- function(data, start, end, breakpoint, lambda){
  n <- nrow(data$X)
  X1 <- data$X[(start+1):breakpoint,,drop = F]
  y1 <- data$y[(start+1):breakpoint]
  X2 <- data$X[(breakpoint+1):end,,drop = F]
  y2 <- data$y[(breakpoint+1):end]

  beta1 <- .lasso_regression(X1, y1, lambda/sqrt(breakpoint-start))
  beta2 <- .lasso_regression(X2, y2, lambda/sqrt(end-breakpoint))

  -1*(as.numeric(.l2norm(X1%*%beta1 - y1)^2) +
        as.numeric(.l2norm(X2%*%beta2 - y2)^2))
}

.create_node <- function(start, end, breakpoint = NA, val = NA, active = NA){
  node <- data.tree::Node$new(paste0(start, "-", end))

  node$start <- start
  node$end <- end
  node$breakpoint <- breakpoint
  node$val <- val
  node$active <- active

  .is_valid(node)

  node
}

.is_valid <- function(obj){
  if(obj$start > obj$end) stop("the start must be less or equal to end")
  if(!is.na(obj$breakpoint) & (obj$start > obj$breakpoint & obj$end < obj$breakpoint))
    stop("breakpoint must be between start and end (inclusive)")

  TRUE
}

.get_leaves_names <- function(tree){
  leaves <- tree$leaves
  vec <- sapply(leaves, function(x){x$name})
  names(vec) <- NULL
  sort(vec)
}

.find_leadingBreakpoint <- function(tree){
  leaves_names <- .get_leaves_names(tree)
  val_vec <- sapply(leaves_names, function(x){
    data.tree::FindNode(tree, x)$val
  })

  leaves_names[which.max(abs(val_vec))]
}

.split_node <- function(node){
  if(is.na(node$breakpoint)) stop("node does not have a set breakpoint yet")
  if(node$breakpoint >= node$end) stop("node breakpoint must be less than end")

  left <- .create_node(node$start, node$breakpoint)
  right <- .create_node(node$breakpoint, node$end)

  list(left = left, right = right)
}

.enumerate_splits <- function(tree){
  names(sort(tree$Get("active")))
}

.jumps <- function(obj){
  leaves <- .enumerate_splits(obj)
  sort(as.numeric(sapply(leaves, function(x){data.tree::FindNode(obj, x)$breakpoint})))
}

