library(shiny)
library(shinyFiles)
library(xcms)
library(gsubfn)

CheckFile<-function(file, type=".mzXML"){
  print(paste("Loading File:", file, sep=""))
  xr <- NULL
  attempt <- 1
  while( is.null(xr) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      xr <- xcmsRaw(file, profstep=0)
    )
  }
  for(i in 1:length(xr@scanindex)){
    scan<-getScan(xr, scan=i)
    if(is.unsorted(scan[,"mz"]) == TRUE){
      newfile<-sub(type, "-Fixed.mzdata", file, ignore.case=TRUE, fixed=TRUE)
      write.mzdata(xr, newfile)
      file.copy(file, sub(type, ".OLD", file, ignore.case=TRUE))
      unlink(file)
      return("corrupt")
      
    }
    if(i == length(xr@scanindex)){
      return("ok")
    }
  }
}

fadeImageButton <- function(inputId, img.path=NA,value=FALSE) {
  # ---------------
  if(is.na(img)) stop("Please enter an image name (place in www folder please)")
  # ---------------
  tagList(
    tags$div(class = "form-group shiny-input-container",
             if(value){
               tags$input(type="checkbox", class="fadebox", id=inputId)
             } else{
               tags$input(type="checkbox", class="fadebox", id=inputId, checked="")
             },
             tags$label(class="btn", "for"=inputId,
                        tags$img(src=img.path, id="btnLeft"))
    )
  )
}

sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)
