# - - helloo - -
"
This script is currently used to start MetaboShiny.
It takes care of installing packages necessary for
MetaboShiny to even start.
"

#' Function to install packages, either through regular method or through downloading from git directly
#' @param package package name to install, either CRAN or bioconductor
install.if.not <- function(packages){
  for(package in packages){
    if(package %in% rownames(installed.packages())){
      NULL
      #print(paste("Already installed base package", package))
    }else{
      if(package %in% c("MetaboAnalystR", "BatchCorrMetabolomics", "MetaDBparse")){
        # Step 1: Install devtools
        install.if.not("devtools")
        # Step 2: Install MetaboAnalystR with documentation
        gitfolder <- switch(package, 
                            MetaboAnalystR = "xia-lab/MetaboAnalystR",
                            BatchCorrMetabolomics = "rwehrens/BatchCorrMetabolomics",
                            MetaDBparse = "UMCUGenetics/MetaDBparse",
                            showtext = "yixuan/showtext")
        devtools::install_github(gitfolder)
      }else{
        BiocManager::install(package)
      }
    }
  }
}

# Install R packages that are required
# TODO: add further package if you need!
needed.packages <- c("BiocManager", "shiny", "shinydashboard", "httr", "curl", "git2r", "devtools",
                     "pacman", "gsubfn", "DT", "R.utils", "data.table", "shinyFiles",
                     "shinyBS", "rhandsontable", "XML", "colorRamps", "enviPat", "shinyalert",
                     "shinyWidgets", "colourpicker", "here", "ECharts2Shiny", "shinyjqui",
                     "later", "shinycssloaders", "qdapDictionaries", "sysfonts", "showtext",
                     "wordcloud2", "Rserve", "RColorBrewer", "xtable", "som", "ROCR",
                     "RJSONIO", "gplots", "e1071", "caTools", "igraph", "randomForest",
                     "Cairo", "pls", "lattice", "rmarkdown", "knitr", "pROC", "Rcpp",
                     "caret", "ellipse", "scatterplot3d", "impute", "pcaMethods",
                     "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", "KEGGgraph",
                     "preprocessCore", "genefilter", "SSPA", "sva", "DBI", "RSQLite",
                     "ggplot2", "minval", "plotly", "pbapply", "sqldf", "plyr", "ChemmineR",
                     "stringr", "heatmaply", "reshape2", "xlsx", "pheatmap", "rJava",
                     "KEGGREST", "manhattanly", "rgl", "glmnet", "TSPred", "VennDiagram",
                     "rcdk", "SPARQL", "webchem", "WikidataQueryServiceR", "openxlsx",
                     "doParallel", "missForest", "InterpretMSSpectrum", "tm", "RISmed",
                     "qdap", "extrafont", "gmp", "shadowtext", "fgsea", "Rmisc", 
                     "BioMedR", "MetaDBparse")

missing.packages <- setdiff(needed.packages,rownames(installed.packages()))

if("BiocManager" %in% missing.packages){
  install.packages('BiocManager')
}

if(length(missing.packages)>0){
  install.if.not(missing.packages)
}

# make metaboshiny_storage dir in home first..
# docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
# with autorun
# docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
# docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
# current instructions
#Rshiny app to analyse untargeted metabolomics data! BASH INSTRUCTIONS: STEP 1: mydir=~"/MetaboShiny" #or another of your choice | STEP 2: mkdir $mydir | STEP 3: docker run -p 8080:8080 -v $mydir:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /start.sh

# packages needed to start up
git.packages <<- c("MetaboAnalystR",
                   "BatchCorrMetabolomics")

# install the base packages needed to start up
for(package in git.packages){
  install.if.not(package)
}

library(httr)

# rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'

runmode <- if(file.exists(".dockerenv")) 'docker' else 'local'

options('unzip.unzip' = getOption("unzip"),
        'download.file.extra' = switch(runmode, 
                                       docker="--insecure",
                                       local=""),  # bad but only way to have internet in docker...
        'download.file.method' = 'curl',
        width = 1200, height=800)

switch(runmode,
       local = {
         wdir <<- dirname(rstudioapi::getSourceEditorContext()$path) # TODO: make this not break when not running from rstudio
         setwd(wdir)
         shiny::runApp(".",launch.browser = T)
       },
       docker = {
         shiny::runApp(".",
                       port = 8080,
                       host = "0.0.0.0",
                       launch.browser = F) # not possible from within docker for now...
       })
