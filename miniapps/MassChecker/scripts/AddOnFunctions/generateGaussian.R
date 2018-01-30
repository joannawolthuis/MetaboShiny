generateGaussian <- function(x,y,resol,plot,scanmode,int.factor,width,height) {

  factor=1.5
  index = which(y==max(y))
  x=x[index]
  y=y[index]
  mu = x
  fwhm = getFwhm(mu,resol)
  x.p = c(mu-factor*fwhm, x, mu+factor*fwhm)
  y.p = c(0, y, 0)
  
#   if (plot) dir.create("./results/plots",showWarnings = FALSE)
#   if (plot) dir.create("./results/plots/Gaussian_fit",showWarnings = FALSE)

  if (plot) {
    if (scanmode=="positive"){
      plot_label="pos_fit.png"
    } else {
      plot_label="neg_fit.png"
    }
  }
  
  mz.range = x.p[length(x.p)] - x.p[1]
  x2 = seq(x.p[1],x.p[length(x.p)],length=mz.range*int.factor)
  sigma = getSD(x.p,y.p)
  scale = optimizeGauss(x.p,y.p,sigma,mu)
  
  if (plot) {
    #CairoPNG(filename=paste("./results/Gaussian_fit",paste(sampname, mu, plot_label, sep="_"), sep="/"), width, height)
    
    plot(x.p,y.p,xlab="m/z",ylab="I", ylim=c(0,1.5*max(y)))
    lines(x2,scale*dnorm(x2,mu,sigma), col="green")
    half_max = max(scale*dnorm(x2,mu,sigma))*0.5
    lines(c(mu - 0.5*fwhm, mu + 0.5*fwhm),c(half_max,half_max),col="orange")
    abline(v = mu, col="green")
    h=c(paste("mean =", mu, sep=" "))
    legend("topright", legend=h)
    
    #dev.off()
  }
  
#  area = sum(scale*dnorm(x2,mu,sigma))
#  area = max(scale*dnorm(x2,mu,sigma))
  area = getArea(mu,resol,scale,sigma,int.factor)
  
  return(list("mean"=mu, "area"=area, "min"=x2[1] , "max"=x2[length(x2)]))

}