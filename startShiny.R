# - - helloo - - 
"
- - -
This script is currently used to start MetaboShiny.
It takes care of installing packages necessary for 
MetaboShiny to even start.
- - -
"


#' Function to install packages, either through regular method or through downloading from git directly
#' @param package package name to install, either CRAN or bioconductor
install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    if(package %in% c("MetaboAnalystR", "BatchCorrMetabolomics", "showtext")){
      metanr_packages() # Installs MetaboAnalyst-specific packages
      # Step 1: Install devtools
      install.if.not("devtools")
      # Step 2: Install MetaboAnalystR with documentation
      gitfolder <- switch(package, MetaboAnalystR = "xia-lab/MetaboAnalystR",
                          BatchCorrMetabolomics = "rwehrens/BatchCorrMetabolomics",
                          showtext = "yixuan/showtext")
      devtools::install_github(gitfolder)#, build_vignettes=TRUE)
    }else{
      install.packages(package)
    }
  }
}


# installs packages that metaboanalyst needs to run
metanr_packages <- function(){
  # - - - - - - - - - - - - - -
  metr_pkgs <- c("Rserve", "BatchCorrMetabolomics", "RColorBrewer", 
                 "xtable", "som", "ROCR", "RJSONIO", "gplots", 
                 "e1071", "caTools", "igraph", "randomForest", "Cairo", 
                 "pls", "pheatmap", "lattice", "rmarkdown", "knitr", 
                 "data.table", "pROC", "Rcpp", "caret", "ellipse",
                 "scatterplot3d", "impute", "rhandsontable", "pcaMethods", 
                 "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", "KEGGgraph", 
                 "preprocessCore", "genefilter", "SSPA", "sva", "showtext", "wordcloud2")
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs) != 0){
    source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
    print(paste0(new_pkgs, " added..."))
  }
  if((length(new_pkgs) < 1)){
    print("No new packages added...")
  }
}

# -- IF YOU HAVE PROBLEMS STARTING UP --
# requires openssl, libcurl, libxt, 
# libxml2, libnetcdf, cairo, java, mesa-common, 
# libgit2, (LIKELY -devel versions)
# libgl1, libglu1, libpng installed

# packages needed to start up
base.packs <- c("httr", "curl", "git2r", "devtools", 
                "pacman", "gsubfn", "shiny", "DT", 
                "R.utils", "data.table", "shinyFiles", 
                "shinyBS", "rhandsontable", "XML", 
                "MetaboAnalystR", "BatchCorrMetabolomics", 
                "colorRamps", "enviPat", "shinyalert",
                "shinyWidgets", "colourpicker", "here",
                "ECharts2Shiny", "shinyjqui", "later")

# install the base packages needed to start up
for(package in base.packs){
  install.if.not(package)
}

# set working directory to where the rstudio file currently is
wdir <<- dirname(rstudioapi::getSourceEditorContext()$path) # TODO: make this not break when not running from rstudio
setwd(wdir)

# go run it! :-)
shiny::runApp(".")#, launch.browser = T)

