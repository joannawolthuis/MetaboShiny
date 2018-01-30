elementInfo <- function(name, elements = NULL) { # from Rdisop function .getElement  
  if (!is.list(elements) || length(elements)==0 ) {
    elements <- initializePSE() }
  if (name=="CH3OH+H"){rex<-"^CH3OH\\+H$"}else{rex <- paste ("^",name,"$", sep="")}
  elements [[grep (rex, sapply (elements, function(x) {x$name}))]]
}