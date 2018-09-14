 load.necessarities <- function(shiny.tab){
 NULL
   }
#   switch(shiny.tab,
#          setup = {
#            library(pacman)
#            library(data.table)},
#          database = {
#            library(RSQLite)
#            library(gsubfn)
#            library(DBI)
#            library(parallel)
#            library(XML)
#            library(minval)
#            library(curl)
#            library(enviPat)
#            library(SPARQL)
#            data(isotopes, package = "enviPat")
#            if(any(is.na(session_cl))){
#              session_cl <<- makeCluster(detectCores())
#              clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
#                "isotopes",
#                "subform.joanna", 
#                "mergeform.joanna",
#                "multiform.joanna",
#                "check.ded.joanna",
#                "data.table",
#                "rbindlist",
#                "isopattern",
#                "keggFind",
#                "keggGet",
#                "kegg.charge",
#                "regexpr",
#                "regmatches"
#              ))
#            }
#          },
#          upload = {
#            library(RSQLite)
#            library(DBI)
#            library(reshape2)
#            library(data.table)
#            library(xlsx)
#          },
#          document = {
#            library(plotly)
#            library(data.table)
#            library(reshape2)},
#          filter = {
#            sourceAll(file.path("backend", 
#                                "scripts", 
#                                "metaboanalyst"))
#            library(plotly)
#            library(data.table)
#            library(Cairo)
#            library(preprocessCore)},
#          analysis = {
#            sourceAll(file.path("backend", 
#                                "scripts", 
#                                "metaboanalyst"))
#            library(plotly)
#            library(heatmaply)
#            library(data.table)},
#          options = {
#            library(shinyFiles)
#            library(colorRamps)
#            library(grDevices)
#            library(colourpicker)
#          })
# }

get.package.table <- function(shiny.tab){
  if(shiny.tab == "setup"){
    status <- sapply(global$constants$packages, FUN=function(package){
      if(package %in% rownames(installed.packages())) "Yes" else "No"
    })
    version <- sapply(global$constants$packages, FUN=function(package){
      if(package %in% rownames(installed.packages())){packageDescription(package)$Version} else ""
    })
    result <-data.table(
      Package = global$constants$packages,
      Installed = status,
      Version = version
    )
  } else result <- NULL
  # --- return ---
  result
}
