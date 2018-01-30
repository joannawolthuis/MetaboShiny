fitG_2 <- function(x,y,sig,mu,scale,useBounds) {
  
    f = function(p) {
      d = p[2]*dnorm(x,mean=p[1],sd=sig)
      sum((d-y)^2)
    }
  
    if (useBounds){
      lower = c(x[1],0,x[1],0)
      upper = c(x[length(x)],Inf,x[length(x)],Inf)
      
      optim(c(as.numeric(mu), as.numeric(scale)),
              f,control=list(maxit=10000),method="L-BFGS-B",lower=lower,upper=upper)
    } else {  
      #optim(c(mu,scale),f)
      optim(c(as.numeric(mu),as.numeric(scale)),f,control=list(maxit=10000))
    }  
}
