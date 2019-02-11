set.seed(10)
dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))

dat$x <- dat$X
library(ChangePointCalc)
zz <- ChangePointCalc::SGLmain(dat, gamma = 1, standardize = F) #?? not sure what's going on


library(CVXR)
X <- dat$X; y <- dat$y
p <- ncol(X); n <- nrow(X); X <- t(X)
lambda=2/(n)^0.5
beta <- Variable(n,p)
D=matrix(0,nrow=n, ncol=n)
diag(D)=1
diag(D[,-1])=-1
D=D[-n,]
# prob <- Problem( Minimize( (cvxr_norm(y-diag( beta%*%X),p=2)^2)/n   +
#                              0.5*lambda*cvxr_norm(D%*%beta,p=1)+
#                              0.5*lambda*cvxr_norm(D%*%beta,p=2)))
prob <- Problem( Minimize( (cvxr_norm(y-diag( beta%*%X),p=2)^2)/n   +
                             lambda*cvxr_norm(D%*%beta,p=1)))

result <- solve(prob,solver="SCS")
beta.hat=result$getValue(beta) #looks like it roughly works, good enough
beta.hat

screening(beta.hat, tau = 5)
