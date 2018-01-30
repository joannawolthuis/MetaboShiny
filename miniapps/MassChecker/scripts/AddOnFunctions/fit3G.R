fit3G_2 <- function(x,y,sig1,sig2,sig3,mu1,scale1,mu2,scale2,mu3,scale3,useBounds){
  
  f = function(p){
    d = p[2]*dnorm(x,mean=p[1],sd=sig1) + p[4]*dnorm(x,mean=p[3],sd=sig2) + p[6]*dnorm(x,mean=p[5],sd=sig3)
    sum((d-y)^2)
  }

  if (useBounds){
    lower = c(x[1],0,x[1],0,x[1],0)
    upper = c(x[length(x)],Inf,x[length(x)],Inf,x[length(x)],Inf)
      
    optim(c(mu1,scale1,mu2,scale2,mu3,scale3),f,control=list(maxit=10000),method="L-BFGS-B",lower=lower,upper=upper)
    
  } else {
    optim(c(mu1,scale1,mu2,scale2,mu3,scale3),f,control=list(maxit=10000))
  }
  
}
