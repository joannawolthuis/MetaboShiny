library(devtools)
library(shiny)

# --- install base base packages ---

install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    install.packages(package)
  }
}

# ----------------------------------

base.packs <- c("pacman", "shiny", "DT", "data.table", "shinyFiles")

for(package in base.packs){
  install.if.not(package)
}

# bioconductor 

source("https://bioconductor.org/biocLite.R")
biocInstaller::biocLite(suppressUpdates = T)

wdir <<- "/Users/jwolthuis/Google Drive/MetaboShiny"
setwd(wdir)

# ---------------------------------
suppressWarnings(
  runApp(".")
)

runApp(".")
