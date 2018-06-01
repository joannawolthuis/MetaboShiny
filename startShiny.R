library(devtools)
library(limma)

# MetaboAnalyst Packages
metanr_packages <- function(){
  
  metr_pkgs <- c("Rserve", "RColorBrewer", "xtable", "som", "ROCR", "RJSONIO", "gplots", "e1071", "caTools", "igraph", "randomForest", "Cairo", "pls", "pheatmap", "lattice", "rmarkdown", "knitr", "data.table", "pROC", "Rcpp", "caret", "ellipse",
                 "scatterplot3d", "impute", "pcaMethods", "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", "KEGGgraph", "preprocessCore", "genefilter", "SSPA", "sva")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs) != 0){
    
    #source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs) < 1)){
    print("No new packages added...")
  }
}

# --- install base base packages ---

install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    if(package == "MetaboAnalystR"){
      metanr_packages()
      # Step 1: Install devtools
      install.packages("devtools")
      # Step 2: Install MetaboAnalystR with documentation
      devtools::install_github("xia-lab/MetaboAnalystR", build_vignettes=TRUE)
    }else{
      install.packages(package)
    }
  }
}

# - - - - - - - - - - - - - - - - -

base.packs <- c("pacman", "shiny", "DT", "data.table", "shinyFiles")

for(package in base.packs){
  install.if.not(package)
}

# bioconductor 

#source("https://bioconductor.org/biocLite.R")
biocInstaller::biocLite(suppressUpdates = T)

wdir <<- "/Users/jwolthuis/Google Drive/MetaboShiny"
setwd(wdir)

# ---------------------------------

shiny::runApp(".")

