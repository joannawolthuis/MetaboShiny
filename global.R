# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)


# ------------------------

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

getOptions <- function(file.loc){
  opt_conn <- file(file.loc)
  # ----------------
  options_raw <<- readLines(opt_conn)
  close(opt_conn)
  # --- list-ify ---
  options <- list()
  for(line in options_raw){
    split  <- (strsplit(line, ' = '))[[1]]
    options[[split[[1]]]] = split[[2]]
  }
  # --- return ---
  options 
}

setOption <- function(file.loc, key, value){
  opt_conn <- file(file.loc)
  # -------------------------
  options <- getOptions(file.loc)
  # --- add new or change ---
  options[[key]] = value
  # --- list-ify ---
  new_options <- lapply(seq_along(options), FUN=function(i){
    line <- paste(names(options)[i], options[i], sep=" = ")
    line
  })
  print(new_options)
  writeLines(opt_conn, text = unlist(new_options))
  close(opt_conn)
}

# --- beta stuff ---

mainmode <- if(exists("dataSet")) dataSet$shinymode else("stat")

# === SOURCE OWN CODE ===

sourceAll <- function(where, 
                      which=c("general", "stats", "time", "enrich_path", "power_roc", "utils")){
  p_load(compiler)
  # ----------------------------
  print("sourcing R code ... ");
  for(i in 1:length(which)){
    script.loc <- file.path(where, which[i])
    print(script.loc)
    files <- list.files(file.path(where, which[i]),full.names=TRUE, pattern=".R$");
    print(files)
    for(f in files){
      print(f)
      source(f)
    }
  }
  return("TRUE");
}

# --------------------------
options <- getOptions(".conf")

sourceDir("backend/scripts/joanna")
load("backend/umcfiles/adducts/AdductTableWKZ.RData")
# === LOAD LIBRARIES ===

library(RSQLite)
library(DBI)
library(reshape2)
library(data.table)
library(xlsx)
library(plotly)
library(preprocessCore)
library(heatmaply)
library(shinyFiles)
library(colorRamps)
library(grDevices)
library(colourpicker)
library(pacman)
library(RSQLite)
library(gsubfn)
library(DBI)
library(parallel)
library(XML)
library(minval)
library(curl)
library(enviPat)
library(SPARQL)
library(KEGGREST)
sourceAll(file.path("backend", 
                    "scripts", 
                    "joanna"))
data(isotopes, package = "enviPat")
if(!exists("session_cl")){
  session_cl <<- makeCluster(detectCores())
  clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
    "isotopes",
    "subform.joanna", 
    "mergeform.joanna",
    "multiform.joanna",
    "check.ded.joanna",
    "data.table",
    "rbindlist",
    "isopattern",
    "keggFind",
    "keggGet",
    "kegg.charge",
    "regexpr",
    "regmatches"
  ))
}
sourceAll(file.path("backend", 
                    "scripts", 
                    "metaboanalyst"))

print("loaded global settings")

