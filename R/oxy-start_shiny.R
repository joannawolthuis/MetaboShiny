#' @title Start MetaboShiny
#' @description Function to start the MetaboShiny app
#' @param port Port to host the app on, Default: 8080
#' @param inBrowser Open in browser automatically?, Default: F
#' @param debug Run in debug mode?, Default: F
#' @param runmode Run locally or on server? (deprecated, ignore), Default: 'local'
#' @seealso 
#'  \code{\link[shiny]{runApp}}
#' @rdname start.metshi
#' @export 
#' @importFrom shiny runApp
start_metshi <- function(port=8080, inBrowser=F, 
                         debug=F, runmode="local"){
  
  options("download.file.method" = "libcurl")
  
  ## make metaboshiny_storage dir in home first..
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
  # NEWEST
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny:dev /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny:dev Rscript -e "MetaboShiny::start_metshi()"
  packages = installed.packages()
  options(repos = getOption("repos")["CRAN"])
  
  
  if(!("MetaboAnalystR" %in% rownames(packages))){
    pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz",
                     "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", 
                     "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"))
    # windows only
    if(.Platform$OS.type == "windows") {
      install.packages("Rserve", repos = "https://cran.r-project.org", type="win.binary")
    }
    urlPackage <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz"
    install.packages(urlPackage, repos=NULL, type="source") 
    
    pacman::p_load(c('Rserve', 'randomForest', 'impute', 'pcaMethods', 'globaltest', 'GlobalAncova', 'Rgraphviz', 'preprocessCore', 'genefilter', 'SSPA', 
                     'sva', 'limma', 'caret', 'KEGGgraph', 'siggenes', 'BiocParallel', 'MSnbase'),character.only = T)
    urlPackage <- "http://cran.nexr.com/src/contrib/Archive/locfit/locfit_1.5-9.tar.gz"
    install.packages(urlPackage, repos=NULL, type="source") 
    
    pacman::p_load("ctc", "locfit")
    pacman::p_load(fgsea)
    
    pacman::p_load(char=c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","httr",
                          "gdata", "glasso", "huge",
                          "robustbase","qqconf",
                          "ppcor","crmn","plotly"))
    
    devtools::install_github("xia-lab/MetaboAnalystR",
                             quiet = F, ref = "0d61192c")
  }
  if(!("MetaDBparse" %in% rownames(packages))){
    pacman::p_load(char = c('anytime', 'rJava', 'ggraph', 'tidygraph', 'WikidataR', 'xmlparsedata'))
    devtools::install_github("lvaudor/glitter", "674418b")
    pacman::p_load(char=c('rJava', 'rvest',
                     'rcdk', 'webchem', 
                     'KEGGREST', 
                     'ChemmineR', 'Rdisop'))
    pacman::p_load(rcdk)
    devtools::install_github("joannawolthuis/MetaDBparse",quiet = F, upgrade = F)
  }
  if(!("MetaboShiny" %in% rownames(packages))){
    urlPackage="https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-8.6.tar.gz"
    install.packages(urlPackage, repos=NULL, type="source") 
    pacman::p_load(char = c("car", "pathview", "DEoptimR", "fpc"))
    remotes::install_github("deepanshu88/shinyDarkmode")
    devtools::install_github("joannawolthuis/MetaboShiny", "dev", quiet = T, upgrade = F)
  }
  if(!("ggVennDiagram" %in% rownames(packages))){
    devtools::install_github("joannawolthuis/ggVennDiagram",quiet = F)
  }
  if(!("BatchCorrMetabolomics" %in% rownames(packages))){
    devtools::install_github("rwehrens/ChemometricsWithR", quiet = T, upgrade = F)
    devtools::install_github("rwehrens/BatchCorrMetabolomics", quiet = T, upgrade = F)
  }
  
  #library(httr)
  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'
  requireNamespace("MetaboAnalystR", quietly = TRUE)
  requireNamespace("ggVennDiagram", quietly = TRUE)
  requireNamespace("MetaDBparse", quietly = TRUE)

  options('unzip.unzip' = getOption("unzip"),
          'download.file.extra' = switch(runmode, 
                                         docker="--insecure",
                                         local=""),  # bad but only way to have internet in docker...
          'download.file.method' = 'curl',
          width = 1200, height=800)
  

  appdir = system.file(package = "MetaboShiny")
  
  runmode <<- runmode
  
  shiny::runApp(
    appDir = appdir,
    launch.browser = inBrowser,
    port = 8080, # DONT REMOVE OR DOCKER STOPS WORKING
    host = "0.0.0.0", # DONT REMOVE OR DOCKER STOPS WORKING
    display.mode = if(debug) "showcase" else "normal")
}

