# sudo R CMD javareconf # FIXES EVERYTHING!!! JUST NEED TO USE ADMIN MODE...
#install.packages("rJava", type="source")


# --- install base base packages ---

install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    if(package %in% c("MetaboAnalystR", "BatchCorrMetabolomics")){
      metanr_packages()
      # Step 1: Install devtools
      install.if.not("devtools")
      # Step 2: Install MetaboAnalystR with documentation
      gitfolder <- switch(package, MetaboAnalystR = "xia-lab/MetaboAnalystR",
                          BatchCorrMetabolomics = "rwehrens/BatchCorrMetabolomics")
      devtools::install_github(gitfolder)#, build_vignettes=TRUE)
    }else{
      install.packages(package)
    }
  }
}

metanr_packages <- function(){
  # - - - - - - - - - - - - - -
  metr_pkgs <- c("Rserve", "BatchCorrMetabolomics", "RColorBrewer", "xtable", "som", "ROCR", "RJSONIO", "gplots", "e1071", "caTools", "igraph", "randomForest", "Cairo", "pls", "pheatmap", "lattice", "rmarkdown", "knitr", "data.table", "pROC", "Rcpp", "caret", "ellipse",
                 "scatterplot3d", "impute", "rhandsontable", "pcaMethods", "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", "KEGGgraph", "preprocessCore", "genefilter", "SSPA", "sva")
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs) != 0){
    
    source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs) < 1)){
    print("No new packages added...")
  }
}

# -- IF YOU HAVE PROBLEMS STARTING UP --
# requires openssl, libcurl, libxt, libxml2, cairo, java installed

base.packs <- c("httr", "curl", "git2r", "devtools", 
                "pacman", "gsubfn", "shiny", "DT", 
                "R.utils", "data.table", "shinyFiles", 
                "shinyBS", "rhandsontable", "XML", 
                "MetaboAnalystR", "BatchCorrMetabolomics", 
                "colorRamps", "enviPat", "shinyalert",
                "shinyWidgets", "colourpicker")

base.packs <- c("xlsx", "rJava", "randomForest", "rgl", "rcdk", "xcms")

for(package in base.packs){
  install.if.not(package)
}

# please put your working dir here 

wdir <<- "/Users/jwolthuis/Google Drive/MetaboShiny"
wdir <<- "/home/joanna/Code/MetaboShiny"

setwd(wdir)

# ---------------------------------

shiny::runApp(".", launch.browser = T)
