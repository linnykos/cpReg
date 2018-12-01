# Dynamic programing

rm(list=ls())
library("genlasso")

setwd("C:/Users/daren/Desktop/SEPP code/")


set.seed(100)

sd=1
S=3

p = 30
B.rep=30
data.index=5
res.repeat=matrix(0,nrow=B.rep, ncol=data.index)
lambda.factor=0.3
res.control=res.repeat

D=matrix(0, nrow=p,ncol=p)
diag(D)=1

lasso.res=function(y,x,p,sd,D){
  n.temp=length(y)
  lambda=lambda.factor*sd*sqrt(n.temp*log(p) )
  out = genlasso(y,x,D)
  beta.temp=coef(out, lambda=lambda )$beta
  res.temp=y- predict.genlasso(out, Xnew=x,lambda=lambda)$fit
  return(0.5*sum(res.temp^2)+  sum( abs( beta.temp ) ) * lambda )
}



for ( nnr in 1: data.index){
  
  n = 30*nnr+60
  change.control=c(0,n/3,2*n/3,n)
  gamma=2*S*log(n*p)
  
  for ( bb in 1:B.rep){
    print(nnr*100+bb )
    x = matrix(rnorm(n*(p-1)),nrow=n)
    x=cbind(rep(1,n),x)
    
        #Beta =[Beta_1,Beta_2,Beta_3]
    Beta1=c(0,rep(-1,3),rep(0,p-4)) 
    
    Beta2= c(0,rep(1,3),rep(0,p-4))
    Beta3= c(0,rep(3,3),rep(0,p-4))
    Beta=c(rep(Beta1,n/3),rep(Beta2,n/3),rep(Beta3,n/3))
    y=rep(0,n)
    X=matrix(0,nrow=n,ncol=n*p)
    
    for( i in 1:n){
      X[i,((i-1)*p+1):(i*p)]=x[i,]
      #XX[i,((i-1)*p+1):(i*p)]=xx[i,]
    }
    for( i in 1:n){
      
      y[i]=X[i,]%*%Beta+rnorm(n=1,mean=0,sd=sd)
      #yy[i]=XX[i,]%*%Beta+rnorm(n=1,mean=0,sd=1)
    }
    
B=rep(0,n)
partition=B

for (r in 1:n){
  B[r]=1/0
 for (l in 1: r)  {
   if ((r-l) <(n/2)){
   ress= ifelse((r-l)>(p/3) ,lasso.res( y=y[l:r],x= x[l:r,],p,sd=sd,D=D),0)
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
beta.hat.dp=c()

for( j in 1:(length(change)-1)){
  n.temp=change[j+1]-change[j]
  coef.temp=coef (genlasso(y[(change[j]+1):change[j+1]],x[(change[j]+1):change[j+1],] ,D),
                  lambda=lambda.factor*sd*sqrt(n.temp*log(p) ) )$beta
  beta.hat.dp= c(beta.hat.dp,
                 rep(coef.temp, n.temp))

}

beta.control=c()
for( j in 1:(length(change.control)-1)){
  n.temp=change.control[j+1]-change.control[j]
  coef.temp=coef ( genlasso(y[(change.control[j]+1):change.control[j+1]],x[(change.control[j]+1):change.control[j+1],] ,D),
                   lambda=lambda.factor*sd*sqrt(n.temp*log(p) ) )$beta
  beta.control= c(beta.control,
                 rep(coef.temp,n.temp ))
}

#plot(beta.hat.dp[seq(1,p*n,3)],ylim=c(0,5))
#points(Beta[seq(1,p*n,3)],col="red")
res.repeat[bb,nnr]=(norm(Beta-beta.hat.dp,type="2")^2)/n
res.control[bb,nnr]=(norm(Beta-beta.control,type="2")^2)/n

#print(length(Beta)==length(beta.control))
#if(length(Beta)!=length(beta.control)){
 # print(beta.control)
#}
print(res.repeat)
print(res.control)
  }
}
#res.repeat
#res.control

res.potts= res.repeat
res.oracle= res.control
save(res.potts, res.oracle, file = "potts.RData")

colMeans(res.potts)
colMeans(res.oracle)
apply(res.repeat,2,sd)
sd(res.repeat[1:B,])
res.repeat/B.rep
res.control/B.rep

beta.temp=coef(out, lambda=0.1*sqrt(n.temp*log(p)))$beta
res.temp=y- predict.genlasso(out, Xnew=x,lambda=0.1*sqrt(n.temp*log(p)))$fit
0.5*sum(res.temp^2)+  sum( abs( beta.temp ) ) * sqrt(n.temp*log(p)) 
