library(shiny)
library(shinyFiles)
library(shinyBS)
library(xcms)
library(gsubfn)
library(data.table)
library(DT)
library(xlsx)
library(doSNOW)
library(tools)
library(MassSpecWavelet)

# ---------------

CheckFile<-function(file, type=".mzXML"){
  print(paste("Loading File:", 
              file, 
              sep=""))
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
      newfile<-sub(type, "-Fixed.mzdata", 
                   file, 
                   ignore.case=TRUE, 
                   fixed=TRUE)
      write.mzdata(xr, newfile)
      file.copy(file, sub(type, 
                          ".OLD", 
                          file, 
                          ignore.case=TRUE))
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

scriptdir <<- "./scripts"
scripts = list.files(path = scriptdir, 
                     pattern = "\\.R",
                     full.names = T,
                     recursive = T)
lapply(scripts, function(x){
  print(basename(x))
  source(x)
})

## DEFAULTS ##
cores <- 1
nrepl <- 3
trim <- 0.1
ppm = 2
dimsThresh <- 100
resol <- 140000
thresh_pos <- 2000
thresh_neg <- 2000
modes <- c("pos", "neg")

# -----------------

sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)