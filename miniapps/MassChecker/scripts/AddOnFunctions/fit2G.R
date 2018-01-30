fit2G_2 <- function(x,y,sig1,sig2,mu1,scale1,mu2,scale2,useBounds){
  
  f = function(p){
    d = p[2]*dnorm(x,mean=p[1],sd=sig1) + p[4]*dnorm(x,mean=p[3],sd=sig2)
    sum((d-y)^2)
  }
  
  if (useBounds){
    lower = c(x[1],0,x[1],0)
    upper = c(x[length(x)],Inf,x[length(x)],Inf)
    
    if (is.null(mu2) && is.null(scale2) && is.null(sig2)){
      sig2=sig1
      optim(c(as.numeric(mu1),
              as.numeric(scale1),
              as.numeric(mu1),
              as.numeric(scale1)),
              f,control=list(maxit=10000),method="L-BFGS-B",lower=lower,upper=upper)
    } else {
      optim(c(as.numeric(mu1),
              as.numeric(scale1),
              as.numeric(mu2),
              as.numeric(scale2)),
              f,control=list(maxit=10000),method="L-BFGS-B",lower=lower,upper=upper)
    }
      
  } else {
    if (is.null(mu2) && is.null(scale2) && is.null(sig2)){
      sig2=sig1
      optim(c(as.numeric(mu1),
              as.numeric(scale1),
              as.numeric(mu1),
              as.numeric(scale1)),
              f,control=list(maxit=10000))
    } else{
      optim(c(as.numeric(mu1),
              as.numeric(scale1),
              as.numeric(mu2),
              as.numeric(scale2)),
              f,control=list(maxit=10000))
    }  
  }
}
