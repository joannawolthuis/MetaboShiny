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
  ## with autorun
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny R -e "MetaboShiny::start.metshi(inBrowser=F, port=8080, runmode='docker')"
  # NEWEST
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny:latest /bin/bash
  
  packages = installed.packages()
  if(!("MetaboAnalystR" %in% rownames(packages))){
    devtools::install_github("xia-lab/MetaboAnalystR",quiet = T, upgrade=F)
  }
  if(!("MetaDBparse" %in% rownames(packages))){
    devtools::install_github("UMCUGenetics/MetaDBparse",quiet = T, upgrade=F)
  }
  if(!("MetaboShiny" %in% rownames(packages))){
    devtools::install_github("UMCUGenetics/MetaboShiny",quiet = T, upgrade=F)
  }
  if(!("ggVennDiagram" %in% rownames(packages))){
    devtools::install_github("joannawolthuis/ggVennDiagram",quiet = T, upgrade=F)
  }
  if(!("BatchCorrMetabolomics" %in% rownames(packages))){
    devtools::install_github("rwehrens/BatchCorrMetabolomics",quiet = T, upgrade=F)
  }
  
  #library(httr)
  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'
  requireNamespace("MetaboAnalystR", quietly = TRUE)
  requireNamespace("ggVennDiagram", quietly = TRUE)
  requireNamespace("MetaDBparse", quietly = TRUE)
  requireNamespace("BatchCorrMetabolomics", quietly = TRUE)
    
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

