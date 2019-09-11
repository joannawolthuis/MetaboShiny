start.metshi <- function(port=8080, inBrowser=F){
  
  options("download.file.method" = "libcurl")
  
  ## make metaboshiny_storage dir in home first..
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
  ## with autorun
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R

  library(httr)
  
  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'
  
  runmode <- if(file.exists(".dockerenv")) 'docker' else 'local'
  
  options('unzip.unzip' = getOption("unzip"),
          'download.file.extra' = switch(runmode, 
                                         docker="--insecure",
                                         local=""),  # bad but only way to have internet in docker...
          'download.file.method' = 'curl',
          width = 1200, height=800)
  
  appdir = system.file(package = "MetaboShiny")

  switch(runmode,
         local = {
           shiny::runApp(
             appDir = appdir,
             launch.browser = inBrowser)
         },
         docker = {
           shiny::runApp(appdir,
                         port = port,
                         host = "0.0.0.0",
                         launch.browser = T)
         })
}

