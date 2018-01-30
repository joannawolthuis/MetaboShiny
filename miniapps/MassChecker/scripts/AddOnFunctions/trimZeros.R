trimZeros <- function(x, y) {
  tmp = which(y==0)
  if (length(tmp)!=0){
    y = y[-tmp]
    x = x[-tmp]
  }
  return(list(x,y))
}
