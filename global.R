# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)


# --- source ---

library(shiny)
library(shinyBS)
library(shinyFiles)
library(MetaboAnalystR)
library(data.table)
library(gsubfn)
library(plotly)
library(colorRamps)
#library(randomForest)

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

mainmode <<- if(exists("dataSet")) dataSet$shinymode else("stat")

# --------------------------
options <- getOptions(".conf")

sourceDir("backend/scripts/joanna")
load("backend/umcfiles/adducts/AdductTableWKZ.RData")

# === LOAD LIBRARIES ===

data(isotopes, package = "enviPat")
session_cl <<- parallel::makeCluster(parallel::detectCores())

source("./Rsource/SwitchButton.R")

sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)

pos_adducts <<- wkz.adduct.confirmed[Ion_mode == "positive",
                                     c("Name")]
neg_adducts <<- wkz.adduct.confirmed[Ion_mode == "negative",
                                     c("Name")]
print("loaded global settings")
