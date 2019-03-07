# - - helloo - - 
"
This script is currently used to start MetaboShiny.
It takes care of installing packages necessary for 
MetaboShiny to even start.
"


#' Function to install packages, either through regular method or through downloading from git directly
#' @param package package name to install, either CRAN or bioconductor
install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    if(package %in% c("MetaboAnalystR", "BatchCorrMetabolomics")){
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

install.packages('BiocManager')

# Install R packages that are required
# TODO: add further package if you need!

needed.packages <- c("shiny", "shinydashboard", "httr", "curl", "git2r", "devtools", 
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
                     "qdap", "extrafont", "gmp", "shadowtext")

missing.packages <- setdiff(needed.packages,rownames(installed.packages()))

print(missing.packages)

if(length(missing.packages)>0){
  BiocManager::install(missing.packages)
}

## used to create the docker run statements
# groupies <- split(unique(to.install), ceiling(seq_along(unique(to.install))/10))
# for(group in groupies){
#   packages <- paste0("'",paste0(group, collapse="', '"), "'")
#   cat(gsubfn::fn$paste('RUN R -e "BiocManager::install(c($packages))"\n'))
# }
## docker internet
# docker run --net=host -it metaboshiny/test1
## docker with volume mounted
# docker run -p 8080:8080 --net=host -it --mount src='~/MetaboShiny',target='/databases/,type=bind metaboshiny/test1

# attempt 1:
# make metaboshiny_storage dir in home first..
# docker run -p 8080:8080 --mount src=~/MetaboShiny_storage,target=/databases/,type=bind --rm -it metaboshiny/test1 sh

# packages needed to start up
git.packages <<- c("MetaboAnalystR", "BatchCorrMetabolomics")

# install the base packages needed to start up
for(package in git.packages){
  install.if.not(package)
}

# go run it! :-)
shiny::runApp(".",
              port = 8080, 
              host = "0.0.0.0", 
              launch.browser = FALSE)
