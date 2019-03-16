cp_fixedstep <- function(data, data_length_func,
                         compute_cusum_func,
                         K,
                         delta = 1,
                         max_candidates = NA,
                         verbose = F, ...){
  stopifnot(K < data_length_func(data))

  #initialization
  n <- data_length_func(data); tree <- .create_node(0, n)

  for(steps in 1:K){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])

      res <- .find_breakpoint(data, c(leaf$start, leaf$end), delta = delta,
                              max_candidates = max_candidates,
                              data_length_func = data_length_func,
                              compute_cusum_func = compute_cusum_func,
                              verbose = verbose, ...)

      leaf$breakpoint <- res$b; leaf$val <- res$val
    }

    node.name <- .find_leadingBreakpoint(tree)
    node.selected <- data.tree::FindNode(tree, node.name)
    node.selected$active <- steps
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
  }

  c(0, .jumps(tree), n)
}

###########


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

