SGL_solver <- function(X, y, lambda = NA){
  p <- ncol(X); n <- nrow(X)
  if(is.na(lambda)) lambda <- 2/n^0.5

  beta <- CVXR::Variable(n,p)
  D <- matrix(0, nrow=n, ncol=n)
  diag(D) <- 1; diag(D[,-1]) <- -1
  D <- D[-n,]

  prob <- CVXR::Problem(CVXR::Minimize((CVXR::cvxr_norm(y-CVXR::diag(X%*%beta),p=2)^2)/n   +
                                         lambda*CVXR::cvxr_norm(D%*%beta,p=1)))

  result <- CVXR::solve(prob, solver="SCS")
  result$getValue(beta)
}