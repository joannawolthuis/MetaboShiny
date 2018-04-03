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
library(enviPat)
library(stringr)
library(BatchCorrMetabolomics)

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
  #print(new_options)
  writeLines(opt_conn, text = unlist(new_options))
  close(opt_conn)
}

# --- beta stuff ---

mainmode <<- if(exists("dataSet")) dataSet$shinymode else("stat")

# --------------------------
options <- getOptions(".conf")

sourceDir("backend/scripts/joanna")
#fwrite(adducts_adj, file = "backend/umcfiles/adducts/AdductTable1.0.csv")
adducts <<- fread("backend/umcfiles/adducts/AdductTable1.0.csv", header = T)

# === LOAD LIBRARIES ===

data(isotopes, package = "enviPat")
session_cl <<- parallel::makeCluster(parallel::detectCores())

source("./Rsource/SwitchButton.R")

sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)

pos_adducts <<- adducts[Ion_mode == "positive",
                                     c("Name")]
neg_adducts <<- adducts[Ion_mode == "negative",
                                     c("Name")]
# interleave for sorting later ...
add_idx <- order(c(seq_along(pos_adducts$Name), seq_along(neg_adducts$Name)))
sort_order <<- unlist(c(pos_adducts$Name,neg_adducts$Name))[add_idx]

bar.css <<- nav.bar.css(options$col1, options$col2, options$col3, options$col4)
font.css <<- font.css(options$font1, options$font2, options$font3, options$font4)
plot.theme <<- ggplot2::theme_minimal

print("loaded global settings")
