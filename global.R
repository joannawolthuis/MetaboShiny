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

# --------------------------
options <- getOptions("user_options.txt")
adducts <- fread("backend/umcfiles/adducts/AdductTable1.0.csv", header = T)
home = normalizePath("~")

sourceDir("backend/scripts/joanna")

# === TEST GLOBAL VARIABLE SCOPE ===

global <- list(constants = list(ppm = 2,
                                nvars = 2,
                                packages = c("data.table", "DBI", "RSQLite", "ggplot2", "minval", "enviPat", 
                                             "plotly", "parallel", "shinyFiles", "curl", "httr", "pbapply", 
                                             "sqldf", "plyr", "ChemmineR", "gsubfn", "stringr", "heatmaply", 
                                             "reshape2", "XML", "xlsx", "colourpicker", "DT", "Rserve", "ellipse", 
                                             "scatterplot3d", "pls", "caret", "lattice", "compiler", "Cairo", 
                                             "randomForest", "e1071", "gplots", "som", "xtable", "RColorBrewer", 
                                             "impute", "pcaMethods", "siggenes", "globaltest", "GlobalAncova", 
                                             "Rgraphviz", "KEGGgraph", "preprocessCore", "genefilter", "pheatmap", 
                                             "igraph", "RJSONIO", "SSPA", "caTools", "ROCR", "pROC", "sva", 
                                             "rJava", "colorRamps", "grDevices", "KEGGREST", "manhattanly", 
                                             "BatchCorrMetabolomics", "R.utils", "rgl", "glmnet", "TSPred", 
                                             "VennDiagram", "rcdk", "SPARQL", "webchem", "WikidataQueryServiceR", 
                                             "openxlsx", "doParallel"),
                                images = list(list(name = 'cute_package', path = 'www/cat.png', dimensions = c(80, 80)),
                                                list(name = 'umc_logo_int', path = 'www/umcinternal.png', dimensions = c(120, 120)),
                                                list(name = 'umc_logo_noise', path = 'www/umcnoise.png', dimensions = c(120, 120)),
                                                list(name = 'hmdb_logo', path = 'www/hmdblogo.png', dimensions = c(150, 100)),
                                                list(name = 'metacyc_logo', path = 'www/metacyc.png', dimensions = c(300, 80)),
                                                list(name = 'chebi_logo', path = 'www/chebilogo.png', dimensions = c(140, 140)),
                                                list(name = 'wikipath_logo', path = 'www/wikipathways.png', dimensions = c(130, 150)),
                                                list(name = 'kegg_logo', path = 'www/kegglogo.gif', dimensions = c(200, 150)),
                                                list(name = 'pubchem_logo', path = 'www/pubchemlogo.png', dimensions = c(145, 90)),
                                                list(name = 'smpdb_logo', path = 'www/smpdb_logo_adj.png', dimensions = c(200, 160)),
                                                list(name = 'dimedb_logo', path = 'www/dimedb_logo.png', dimensions = c(310, 120)),
                                                list(name = 'wikidata_logo', path = 'www/wikidata.png', dimensions = c(250, 200)),
                                                list(name = 'vmh_logo', path = 'www/vmh.png', dimensions = c(250, 200)),
                                                list(name = 'pos_icon', path = 'www/handpos.png', dimensions = c(120, 120)),
                                                list(name = 'neg_icon', path = 'www/handneg.png', dimensions = c(120, 120)),
                                                list(name = 'excel_icon', path = 'www/excel.png', dimensions = c(120, 120)),
                                                list(name = 'excel_icon_2', path = 'www/excel.png', dimensions = c(120, 120)),
                                                list(name = 'db_icon', path = 'www/servers.png', dimensions = c(150, 150)),
                                                list(name = 'csv_icon', path = 'www/office.png', dimensions = c(100, 100)),
                                                list(name = 'dataset_icon', path = 'www/office.png', dimensions = c(100, 100))
                                ),
                                default.text = list(list(name='work_dir',text=options$work_dir),
                                                    list(name='curr_db_dir',text=options$db_dir),
                                                    list(name='ppm',text=options$ppm),
                                                    list(name='analUI',text="Please choose an analysis mode!"),
                                                    list(name='proj_name',text=options$proj_name),
                                                    list(name="curr_cpd", text="...")
                                )
                           ),
               functions = list(cf = rainbow,
                                color.function = rainbow),
               paths = list(patdb = file.path(options$work_dir, paste0(options$proj_name, ".db")),
                            csv_loc = file.path(options$work_dir, paste0(options$proj_name, ".csv")),
                            volumes =  c('MetaboShiny' = getwd(),
                                         'Home'=home,
                                         '~' = normalizePath("~"),
                                         'Downloads'=file.path(home, "Downloads"),
                                         'R Installation'=R.home(),
                                         'Desktop'=file.path(home, "Desktop"))
               ),
               tables = list(tbl = NA),
               vectors = list(db_list = c("internal", 
                                          "noise", 
                                          "hmdb", 
                                          "chebi", 
                                          "kegg",
                                          "metacyc",
                                          "wikipathways"
                                          ,"smpdb",
                                          "dimedb",
                                          "wikidata",
                                          "vmh"
                              ),
                              pos_adducts = adducts[Ion_mode == "positive",
                                                     c("Name")],
                              neg_adducts = adducts[Ion_mode == "negative",
                                                    c("Name")]
                              )
               )

# === LOAD LIBRARIES ===

data(isotopes, package = "enviPat")
session_cl <- parallel::makeCluster(max(c(1, parallel::detectCores()-1)))

source("./Rsource/SwitchButton.R")

sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)

# interleave for sorting later ...
add_idx <- order(c(seq_along(global$vectors$pos_adducts$Name), seq_along(global$vectors$neg_adducts$Name)))
sort_order <<- unlist(c(global$vectors$pos_adducts$Name, global$vectors$neg_adducts$Name))[add_idx]

bar.css <<- nav.bar.css(options$col1, options$col2, options$col3, options$col4)
font.css <<- font.css(options$font1, options$font2, options$font3, options$font4)
plot.theme <<- ggplot2::theme_minimal
taskbar_image <<- options$task_img

print("loaded global settings")
