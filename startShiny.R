
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

library(devtools)
library(shiny)

# bioconductor 

source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = T)

wdir <<- "C:/Users/joby/Software/MetaboShiny"
setwd(wdir)
wdir
# ---------------------------------

runApp(".",launch.browser = T)
