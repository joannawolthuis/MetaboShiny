start.metshi <- function(port=8080, inBrowser=F, 
                         debug=F, runmode="local"){
  
  options("download.file.method" = "libcurl")
  
  ## make metaboshiny_storage dir in home first..
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
  ## with autorun
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R

  library(httr)
  
  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'
  
  options('unzip.unzip' = getOption("unzip"),
          'download.file.extra' = switch(runmode, 
                                         docker="--insecure",
                                         local=""),  # bad but only way to have internet in docker...
          'download.file.method' = 'curl',
          width = 1200, height=800)
  
  appdir = system.file(package = "MetaboShiny")

  switch(runmode,
         local = {
           MetaboShiny::start_orca(9091)
           
           shiny::runApp(
             appDir = appdir,
             launch.browser = inBrowser,
             display.mode = if(debug) "showcase" else "normal")
         },
         docker = {
           plotly::orca_serve(port = 9091)
           
           shiny::runApp(appdir,
                         port = port,
                         host = "0.0.0.0",
                         launch.browser = T,
                         display.mode = if(debug) "showcase" else "normal")
         })
}

