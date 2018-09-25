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
