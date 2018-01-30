optimizeGauss <- function(x,y,sigma,mu) {

  f = function(p,x,y,sigma,mu) {
    curve = p*dnorm(x,mu,sigma)
    return((max(curve)-max(y))^2)
  }
  
  rval = optimize(f, c(0, 100000), tol = 0.0001,x,y,sigma,mu)
  
  return(rval$minimum)
}
