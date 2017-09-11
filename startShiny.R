library(devtools)

# --- install base base packages ---

install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    install.packages(package)
  }
}

# ----------------------------------

base.packs <- c("pacman", "shiny", "Bioconductor")

for(package in base.packs){
  install.if.not(package)
}

library(pacman)
library(shiny)

# bioconductor 

source("https://bioconductor.org/biocLite.R")
biocLite()

# ---------------------------------

runApp(".")