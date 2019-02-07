#Examples
rm(list=ls())
setwd("C:/Users/daren/Desktop/simulations/")

library(CVXR)
library(genlasso)
set.seed(100)

p=40
#lambda.factor=10

B.rep=30
data.index=5
lambda.par=1
potts.res=matrix(0,nrow=B.rep, ncol=data.index)
potts.haus=potts.res

fused.res=potts.res
fused.haus=fused.res



lasso.res=function(y,x,p){
  
  D=matrix(0, nrow=p,ncol=p)
  diag(D)=1
  n.temp=length(y)
  lambda=lambda.par*sqrt(n.temp*log(p) )
  out = genlasso(y,x,D)
  beta.temp=coef(out, lambda=lambda )$beta
  res.temp=y- predict.genlasso(out, Xnew=x,lambda=lambda)$fit
  return(0.5*sum(res.temp^2)+  sum( abs( beta.temp ) ) * lambda )
}



haus.dist=function(a,b){
  N=length(a)
  M=length(b)
  dist.matrix=matrix(0,nrow=N, ncol=M)
  for ( i in 1:N){
    for ( j in 1:M){
      dist.matrix[i,j]=abs(a[i]-b[j])
      
    }
    
  }
  
  c1=max(apply(dist.matrix,2,min))
  c2=max(apply(t(dist.matrix),2,min))
  
  return(max(c1,c2))}

for ( nnr in 1: data.index){
  
  n = 60+nnr*60
  change.control=c(0,n/3,2*n/3,n)
  Beta1=c(rep(-1,3),rep(0,p-3)) 
  
  Beta2= c(rep(1,3),rep(0,p-3))
  Beta3= c(rep(3,3),rep(0,p-3))
  
  Beta1=t(matrix(rep(Beta1,n/3),ncol=n/3,nrow=p))
  Beta2=t(matrix(rep(Beta2,n/3),ncol=n/3,nrow=p))
  Beta3=t(matrix(rep(Beta3,n/3),ncol=n/3,nrow=p))
  
  Beta=rbind(Beta1,Beta2,Beta3)
  for ( bb in 1:B.rep){
    print(nnr*100+bb )
    
    X=matrix(rnorm(n*p),nrow=p)
    y=rep(0,n)
    
    
    for( i in 1:n){
      
      y[i]=sum(X[,i]*Beta[i,])+rnorm(n=1,mean=0,sd=1)
      #yy[i]=XX[i,]%*%Beta+rnorm(n=1,mean=0,sd=1)
    }
    
    
    D=matrix(0,nrow=n, ncol=n)
    diag(D)=1
    diag(D[,-1])=-1
    D=D[-n,]
    #end of generate data  
    
    
    

    change.true=c(0,n/3,2*n/3,n)
   
   
    
    #Fused
    beta <- Variable(n,p)


    
    lambda=2/(n)^0.5
    prob <- Problem( Minimize( (cvxr_norm(y-diag( beta%*%X),p=2)^2)/n   +
                                0.5*lambda*cvxr_norm(D%*%beta,p=1)+  
                               0.5*lambda*cvxr_norm(D%*%beta,p=2)))
    
    result <- solve(prob,solver="SCS")
    beta.hat=result$getValue(beta)
    
    print(norm( beta.hat-Beta,type="2")^2 /n)
    
    fused.res[bb,nnr]=norm( beta.hat-Beta,type="2")^2 /n
    
    change.estimate = c(0,which(apply(abs(D%*% beta.hat ) ,1,function(x) norm(x, type="2"))>0.5),n)
    
    fused.haus[bb,nnr] =haus.dist(change.true,change.estimate)/n
    
    
    #potts functional
    B=rep(0,n)
    partition=B
    gamma=6*log(n*p)
    for (r in 1:n){
      B[r]=1/0
      for (l in 1: r)  {
        if ((r-l) <(n/2)){
          ress= ifelse((r-l)>10 ,lasso.res( y=y[l:r],x= t(X[,l:r]),p=p),0)
          b = ifelse(l-1==0,-1*gamma,B[l-1] ) + gamma+ ress
          
          if (b < B[r] ){
            B[r]=b
            partition[r]=l-1
          }
        }
      }
      
    }
    
    change=c(n)
    while(min(change)!= 0){
      change=c(change,partition[min(change)]) }
    change=sort(change)
    
    change
    beta.hat.dp=matrix(nrow=0,ncol=p)
    D=matrix(0, nrow=p,ncol=p)
    diag(D)=1
    for( j in 1:(length(change)-1)){
      n.temp=change[j+1]-change[j]
      coef.temp=coef (genlasso(y[(change[j]+1):change[j+1]],t(X[,(change[j]+1):change[j+1]]) ,D),
                      lambda=lambda.par*sqrt(n.temp*log(p) ) )$beta
      #print(coef.temp)
      beta.hat.dp= rbind(beta.hat.dp,
                     t(matrix(rep(coef.temp, n.temp), nrow=p)))
      
    }
    
    print(norm( beta.hat.dp-Beta,type="2")^2 /n)
    
    potts.res[bb,nnr]=norm( beta.hat.dp-Beta,type="2")^2 /n
    
    potts.haus[bb,nnr] =haus.dist(change.true,change)/n
    
    print(change)
    
    
  }
  
 
    
    
  }

colMeans(fused.res)


colMeans(fused.haus)
colMeans(potts.haus)
colMeans(potts.res)
save.image()
N= 60+c(1:5)*60

MSE.plot (Data=fused.res,the.color="blue",dimension=40, N=N,ylim=c(0,2), lty=1, 
          title="coefficient estimation ", ylabel = "MSE")
MSE.points (Data=potts.res,the.color="black", N=N,lty=1)

legend(y=2,x=300,legend=c("SGL","Potts functional"),lty=c(1,1),
       cex=1, col=c(  "blue", "black"))

MSE.plot (Data=fused.haus,the.color="blue",dimension=40, N=N,ylim=c(0,0.2), lty=1, 
          title="change point localization ", ylabel = "Haus")
MSE.points (Data=potts.haus,the.color="black", N=N,lty=1)

legend(y=0.2,x=300,legend=c("SGL","Potts functional"),lty=c(1,1),
       cex=1, col=c(  "blue", "black"))

