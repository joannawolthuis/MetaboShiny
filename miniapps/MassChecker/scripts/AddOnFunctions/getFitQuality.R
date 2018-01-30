getFitQuality <- function(x,y,muFirst,muLast,resol,scale=NULL,sigma=NULL,sumFit=NULL){
  # muFirst=sort(c(p5[1],p5[3],p5[5],p5[7]))[1]
  # muLast=sort(c(p5[1],p5[3],p5[5],p5[7]))[4]

  # muFirst=p2[1]
  # muLast=p2[1]
  # scale=p2[2]
  
  # x,y,sort(c(p5[1],p5[3],p5[5],p5[7]))[1],sort(c(p5[1],p5[3],p5[5],p5[7]))[4],resol,sumFit=sumFit
  
  if (is.null(sumFit)){
    # fwhm = getFwhm(muFirst,resol)
    # intRangeMin = muFirst - 1.2*fwhm
    # intRangeMax = muFirst + 1.2*fwhm
    # tmp = which((x>intRangeMin) & (x<intRangeMax))
    
    # # mu outside range!!!
    # if (length(tmp)==0) return(list("fq_new"=1, "x_int"=x, "y_int"=y))
    # 
    # x_int = x[tmp]
    # y_int = y[tmp] 
    
    x_int = x
    y_int = y
    
    # fq_new = mean(abs((scale*dnorm(x_int,muFirst,sigma)) - y_int)/(scale*dnorm(x_int,muFirst,sigma)))
    fq_new = mean(abs((scale*dnorm(x_int,muFirst,sigma)) - y_int)/rep((max(scale*dnorm(x_int,muFirst,sigma))/2),length(x_int)))
    
    # lines(x_int,scale*dnorm(x_int,muFirst,sigma), col="blue")
    # fq_new = abs(sum(y) - sum(p2[2]*dnorm(x,p2[1],sigma)))/sum(y)
    
  } else {
    # fwhm = getFwhm(muFirst,resol)
    # intRangeMin = muFirst - 1.2*fwhm
    # fwhm = getFwhm(muLast,resol)
    # intRangeMax = muLast + 1.2*fwhm
    # 
    # tmp = which((x>intRangeMin) & (x<intRangeMax))
    # 
    # if (length(tmp)==0){
    #   # outside range
    #   fq_new = 1
    #   x_int = x
    #   y_int = y
    #   
    # } else {
    
      # x_int = x[tmp]
      # y_int = y[tmp]
      # sumFit_int = sumFit[tmp] 
        
      sumFit_int = sumFit
      y_int = y
      x_int = x
      
      # fq_new = mean(abs(sumFit_int - y_int)/sumFit_int)
      fq_new = mean(abs(sumFit_int - y_int)/rep(max(sumFit_int)/2,length(sumFit_int)))
    # }    
  }
  
  # Caused by dividing by 0
  if (is.nan(fq_new)) fq_new=1 
  
  return(list("fq_new"=fq_new, "x_int"=x_int, "y_int"=y_int))  
}
