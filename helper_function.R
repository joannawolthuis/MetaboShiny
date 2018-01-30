listFunctions_inner <- function(
  name, 
  do.recursive=FALSE,
  .do.verbose=FALSE,
  .buffer=new.env()
){
  ..name  <- "listFunctions_inner"
  if (!is.character(name) | missing(name)) {
    stop(paste(..name, " // expecting 'name' of class 'character'", sep=""))
  }
  name.0 <- name
  if (tryCatch(is.function(get(name)), error=function(e) FALSE)) {
    # PROCESS FUNCTIONS       
    if (.do.verbose) {
      message(paste(..name, " // processing function: '", name, "'", sep=""))
    } 
    # Get the function's code:
    code <- deparse(get(name))  
    # break code up into sections preceding left brackets:
    left.brackets <- c(unlist(strsplit(code, split="[[:space:]]*\\(")))  
    out <- sort(unique(unlist(lapply(left.brackets, function (x) {
      # Split up according to anything that can't be in a function name.
      # split = not alphanumeric, not '_', and not '.'
      words <- c(unlist(strsplit(x, split="[^[:alnum:]_.]")))
      
      last.word <- tail(words, 1)
      last.word.is.function <- tryCatch(is.function(get(last.word)),
                                        error=function(e) return(FALSE))
      out <- last.word[last.word.is.function]
      return(out)
    }))))
    if (do.recursive){           
      # funs.checked: We need to keep track of which functions 
      # we've checked to avoid infinite loops.
      .buffer$funs.checked <- c(.buffer$funs.checked, name)
      funs.next <- out[!(out %in% .buffer$funs.checked)]        
      if (length(funs.next)) {
        out <- sort(unique(unlist(c(out, do.call(c,
                                                 lapply(funs.next, function(x) {
                                                   if (x == ".Primitive") {
                                                     return(NULL)
                                                   }
                                                   listFunctions_inner(
                                                     name=x, 
                                                     do.recursive=TRUE,
                                                     .buffer=.buffer
                                                   )
                                                 })
        )))))            
      }
    } 
    out <- sort(unique(unlist(out)))
  } else {
    # PROCESS NAMESPACES
    if (.do.verbose) {
      message(paste(..name, " // processing namespace: '", name, "'", sep=""))
    }
    name    <- paste("package", ":", name, sep="")
    if (!name %in% search()) {
      stop(paste(..name, " // invalid namespace: '", name.0, "'"))
    }
    # KEEP AS REFERENCE       
    #        out <- ls(name)
    funlist <- lsf.str(name)
    out     <- head(funlist, n=length(funlist))
  }
  out
}