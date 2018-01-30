checkOverlap <- function(range1,range2){
  if (length(intersect(range1,range2))==2) {
    # Overlap
    # message("Overlap, smaller range is used")
    if (length(range1) >= length(range2)){
      range1=range1[-length(range1)]  
    } else {
      range2=range2[-1]
    }
  } else if (length(intersect(range1,range2))==3){
    # message("Overlap, smaller range is used")
    range1=range1[-length(range1)]  
    range2=range2[-1]
  }
  return(list("range1"=range1,"range2"=range2))
}