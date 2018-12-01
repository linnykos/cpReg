rm(list=ls())
#setwd("C:/Users/daren/Desktop/SEPP code/")
library(ncvreg)
set.seed(1)

sd=1

#dimension
p = 4

# number of repeatition of each n
B=10

#number of different n
data.index=5

#res.repeat is the matrix to record ress
res.repeat=matrix(0, ncol=data.index,nrow=B)
for ( nnr in 1: data.index){

  # data size
  n = 30*nnr+60

  for ( bb in 1:B){


    print(nnr*100+bb )

    #design matrix
    # the package force user to have a slope, and there is no way to remove it
    #so the design matrix need to have a slope, but the corresponding beta is set to be 0.
    x = matrix(rnorm(n*(p-1)),nrow=n)
    #xx = matrix(rnorm(n*(p-1)),nrow=n)
    x=cbind(rep(1,n),x)
    #xx=cbind(rep(1,n),xx)

    #Beta: the change points are n/3 and 2n/3 as Beta suggests.
    Beta1=c(0,rep(-1,3),rep(0,p-4))

    Beta2= c(0,rep(1,2),rep(0,p-3))
    Beta3= c(0,rep(3,3),rep(0,p-4))
    Beta=c(rep(Beta1,n/3),rep(Beta2,n/3),rep(Beta3,n/3))

    # y is the data
    y=rep(0,n)

    #yy is testing set. Without using the  testing set (so that the effective data size is larger than n),
    # it is extremely unclear how to choose the tuning parameter
    yy=rep(0,n)

    # X is the matrix to generate y only
    X=matrix(0,nrow=n,ncol=n*p)
    #XX=X
    for( i in 1:n){
      X[i,((i-1)*p+1):(i*p)]=x[i,]
      #XX[i,((i-1)*p+1):(i*p)]=x[i,]
    }

    # generate y
    for( i in 1:n){
      y[i]=X[i,]%*%Beta+rnorm(n=1,mean=0,sd=sd)
      yy[i]=X[i,]%*%Beta+rnorm(n=1,mean=0,sd=sd)
    }

    #X.tilde is the matrix to fit in the function ncvreg
    X.tilde=matrix(nrow=0, ncol=n*p)
    #XX.tilde=matrix(nrow=0, ncol=n*p)

    for( i in 1:n){
      X.tilde=rbind(X.tilde, c(rep(x[i,],i), rep(0, n*p-p*i)))
      #XX.tilde=rbind(XX.tilde, c(rep(xx[i,],i), rep(0, n*p-p*i)))
    }

    #the coefficient of fit estimate theta_i =beta_i -beta_i
    fit <- ncvreg(X=X.tilde[,-1], y=y, penalty=c("SCAD"), alpha=1,
                  lambda.min = 0.0001, nlambda=100, eps=1e-5,
                  max.iter = 25000)

    #validation using testing data
    MM=length(fit$beta[1,])
    res.cross=rep(0,MM)
    for ( m in 1:MM){
      res.cross[m]= sum((yy-X.tilde%*%fit$beta[,m])^2)
    }
    mm=which.min(res.cross)
    #MM
    #mm
    #beta.hat = DD%*%fit$beta[seq(1,n*p,4)+1,which.min(res.cross)]


    #plot(beta.hat)
    #points(Beta[seq(1,480,4)+1])



    # find beta hat from theta hat
    DD=matrix(0,nrow=n,ncol=n)
    DD[lower.tri(DD, diag = T)]=1


    res.temp=0
    for( j in 1:p){
      beta.hat = DD%*%fit$beta[seq(1,n*p,p)+j-1,which.min(res.cross)]
      res.temp= res.temp+norm((beta.hat -Beta[seq(1,n*p,p)+j-1]),type="2")^2/n

    }
    res.repeat[bb,nnr]=res.temp
    #print(which.min(res.cross))
  }
}
res.scad =res.repeat
save(res.scad,file="scad.RData")

colMeans(res.scad)
apply(res.scad,2,sd)

#plot(DD%*%fit$beta[seq(1,n*p,p)+p-1,which.min(res.cross)],col="red",ylim=c(-10,10))
#points(Beta[seq(1,n*p,p)+p-1])
