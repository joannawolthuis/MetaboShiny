# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

### Anything loaded in this GLOBAL file will be avaliable for metaboShiny in general. ###

# --- source neccessary libraries ---

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

#' Sources all R scripts in a given directory. 
#' 
#' \code{sourceDir} searches the given directory for .R files and sources them into the current session.
#' 
#' @param path Path to search for scripts in.
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    #if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    #if(trace) cat("\n")
  }
  
}

#' Load user options saved in file.
#' 
#' \code{getOptions} returns all current user options defined in the given options file.
#' 
#' @param file.loc Path to user options file to read in.
#' @return R list with keys as option types and values as option values.
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

#' Changes an option in the given MetaboShiny user options file. 
#' 
#' \code{setOption} changes an option in the options file. Can also be used to add new options to the file.
#' 
#' @param file.loc Location of user options file. Usually .txt but any format is fine.
#' @param key Name of the new option / to change option
#' @param value Value of the option to change or add
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

# load in options
options <- getOptions("user_options.txt")

# load adduct table (if you add/remove any adducts, change them in the below file!)
#adducts <- fread("backend/umcfiles/adducts/AdductTable2.0.csv", header = T)
adducts <- fread("backend/umcfiles/adducts/AdductTable2.0.csv", header = T) # V2 has di/trimers

# set the home path
home = normalizePath("~")

# source all used functions (see shiny_plot.R, shiny_general.R, shiny_db.R etc.)
sourceDir("backend/scripts/joanna")

# === THE BELOW LIST CONTAINS ALL GLOBAL VARIABLES THAT METABOSHINY CALLS UPON LATER ===

global <- list(constants = list(ppm = 2, # TODO: re-add ppm as option for people importing their data through csv
                                timeseries = FALSE,
                                nvars = 2, # Default bivariate setting for statistics
                                max.cols = 8, # Maximum colours available to choose (need to change if anyone does ANOVA with >8 variables)
                                packages = unique(c(base.packs, "data.table", "DBI", "RSQLite", "ggplot2", "minval", "enviPat", 
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
                                             "openxlsx", "doParallel", "missForest", "InterpretMSSpectrum",
                                             "tm", "RISmed", "qdap", "extrafont", "sysfonts", "gmp", "shadowtext")
                                ), # these packages are listed in the first tab and should include all necessary packages
                                images = list(list(name = 'load_icon', path = 'www/cute.png', dimensions = c(100, 100)),
                                              list(name = 'cute_package', path = 'www/cat.png', dimensions = c(80, 80)),
                                              list(name = 'internal_logo', path = 'www/umcinternal.png', dimensions = c(120, 120)),
                                              list(name = 'noise_logo', path = 'www/umcnoise.png', dimensions = c(120, 120)),
                                              list(name = 'hmdb_logo', path = 'www/hmdblogo.png', dimensions = c(150, 100)),
                                              list(name = 'metacyc_logo', path = 'www/metacyc.png', dimensions = c(300, 80)),
                                              list(name = 'chebi_logo', path = 'www/chebilogo.png', dimensions = c(140, 140)),
                                              list(name = 'wikipathways_logo', path = 'www/wikipathways.png', dimensions = c(130, 150)),
                                              list(name = 'kegg_logo', path = 'www/kegglogo.gif', dimensions = c(200, 150)),
                                              list(name = 'pubchem_logo', path = 'www/pubchemlogo.png', dimensions = c(145, 90)),
                                              list(name = 'smpdb_logo', path = 'www/smpdb_logo_adj.png', dimensions = c(200, 160)),
                                              list(name = 'dimedb_logo', path = 'www/dimedb_logo.png', dimensions = c(310, 120)),
                                              list(name = 'wikidata_logo', path = 'www/wikidata.png', dimensions = c(250, 200)),
                                              list(name = 'respect_logo', path = 'www/respect_logo.png', dimensions = c(250, 100)),
                                              list(name = 'massbank_logo', path = 'www/massbank_logo.jpg', dimensions = c(250, 100)),
                                              list(name = 'metabolights_logo', path = 'www/metabolights_logo.png', dimensions = c(200, 200)),
                                              list(name = 'vmh_logo', path = 'www/vmh.png', dimensions = c(250, 200)),
                                              list(name = 'foodb_logo', path = 'www/foodb_logo.png', dimensions = c(200, 200)),
                                              list(name = 'pos_icon', path = 'www/handpos.png', dimensions = c(120, 120)),
                                              list(name = 'neg_icon', path = 'www/handneg.png', dimensions = c(120, 120)),
                                              list(name = 'excel_icon', path = 'www/excel.png', dimensions = c(120, 120)),
                                              list(name = 'excel_icon_2', path = 'www/excel.png', dimensions = c(120, 120)),
                                              list(name = 'db_icon', path = 'www/servers.png', dimensions = c(150, 150)),
                                              list(name = 'lipidmaps_logo', path = 'www/lipidmaps.png', dimensions = c(200, 150)),
                                              list(name = 'bloodexposome_logo', path = 'www/bloodexposome.png', dimensions = c(250, 200)),
                                              list(name = 'csv_icon', path = 'www/office.png', dimensions = c(100, 100)),
                                              list(name= 'magicball', path = 'www/magic-ball2.png', dimensions = c(200,200)),
                                              list(name = 'dataset_icon', path = 'www/office.png', dimensions = c(100, 100)),
                                              list(name = 'plus', path = 'www/add.png', dimensions = c(150, 150)),
                                              list(name= 'maconda_logo', path = 'www/maconda.png', dimensions = c(250,100)),
                                              list(name = 'sidebar_icon', path = 'www/detective.png', dimensions = c(60, 60))
                                              
                                ), # all image paths, if you add an image you can add it here
                                default.text = list(list(name='curr_exp_dir',text=options$work_dir),
                                                    list(name='curr_db_dir',text=options$db_dir),
                                                    list(name='curr_definition', text="No m/z selected"),
                                                    list(name='ppm',text=options$ppm),
                                                    list(name='proj_name',text=options$proj_name),
                                                    list(name="curr_cpd", text="..."),# default text options at startup
                                                    list(name="ml_train_ss", text="all"),
                                                    list(name="ml_test_ss", text="all")
                                ),
                                font.aes = list(font = options$font4,
                                                ax.num.size = 11,
                                                ax.txt.size = 15,
                                                ann.size = 20,
                                                title.size = 25),
                                db.build.info = list(
                                                    # internal = list(title = "UMC Internal", 
                                                    #                  description = "Internal commonly known metabolites",
                                                    #                  image_id = "internal_logo"),
                                                    #  noise = list(title = "UMC Noise", 
                                                    #               description = "Internal common pollutants found in DIMS using local method.",
                                                    #               image_id = "noise_logo"),
                                                     hmdb = list(title = "HMDB", 
                                                                 description = "Metabolites commonly found in human biological samples.",
                                                                 image_id = "hmdb_logo"),
                                                     metacyc = list(title = "MetaCyc", 
                                                                    description = "Large pathway database with over 10 000 available compounds. Spans several organisms!",
                                                                    image_id = "metacyc_logo"),
                                                     chebi = list(title = "ChEBI", 
                                                                  description = "A broad database with known chemicals of biological interest.",
                                                                  image_id = "chebi_logo"),
                                                     wikipathways = list(title = "WikiPathways", 
                                                                         description = "Open source biological pathway database. Currently only partially available. Requires CHEBI to be built.",
                                                                         image_id = "wikipathways_logo"),
                                                     kegg = list(title = "KEGG", 
                                                                 description = "Large pathway database with info on pathways in various organisms, involved enzymes, and connected disease phenotypes.",
                                                                 image_id = "kegg_logo"),
                                                     # pubchem = list(title = "PubChem", 
                                                     #                description = "Internal commonly known metabolites",
                                                     #                image_id = "internal_logo"),
                                                     smpdb = list(title = "SMPDB", 
                                                                  description = "Small molecule pathway database. Compounds overlap with HMDB.",
                                                                  image_id = "smpdb_logo"),
                                                     dimedb = list(title = "DimeDB", 
                                                                   description = "A direct infusion database of biologically relevant metabolite structures and annotations.",
                                                                   image_id = "dimedb_logo"),
                                                     wikidata = list(title = "Wikidata", 
                                                                     description = "Central storage for the data of its Wikimedia sister projects including Wikipedia, Wikivoyage, Wikisource, and others.",
                                                                     image_id = "wikidata_logo"),
                                                     respect = list(title = "ReSpect", 
                                                                    description = "RIKEN MSn spectral database for phytochemicals (ReSpect) is a collection of literature and in-house MSn spectra data for research on plant metabolomics.",
                                                                    image_id = "respect_logo"),
                                                     massbank = list(title = "MassBank", 
                                                                     description = "This site presents the database of comprehensive, high-resolution mass spectra of metabolites. Supported by the JST-BIRD project, it offers various query methods for standard spectra from Keio Univ., RIKEN PSC, and others. 
                                                                     In 2008, MassBank was authorized as the official mass spectral database of The Mass Spectrometry Society of Japan.",
                                                                     image_id = "massbank_logo"),
                                                     metabolights = list(title = "MetaboLights", 
                                                                         description = "MetaboLights is a database for Metabolomics experiments and derived information. The database is cross-species, cross-technique and covers metabolite structures and their reference spectra as well as their biological roles, locations and concentrations, and experimental data from metabolic experiments.",
                                                                         image_id = "metabolights_logo"),
                                                     vmh = list(title = "VMH",
                                                                description = "Virtual Metabolic Human (VMH) hosts ReconMap, an extensive network of human metabolism, and bacterial metabolites.",
                                                                image_id = "vmh_logo"),
                                                     foodb = list(title = "FooDB",
                                                                  description = "FooDB is the world’s largest and most comprehensive resource on food constituents, chemistry and biology. It provides information on both macronutrients and micronutrients, including many of the constituents that give foods their flavor, color, taste, texture and aroma.",
                                                                  image_id = "foodb_logo"),
                                                     bloodexposome = list(title = "Blood Exposome DB",
                                                                  description = " This new blood exposome database can be applied to prioritize literature-based chemical reviews, developing target assays in exposome research, identifying compounds in untargeted mass spectrometry and biological interpretation in metabolomics data.",
                                                                  image_id = "bloodexposome_logo"),
                                                     maconda = list(title = "MaConDa", 
                                                                    description = "MaConDa currently contains ca. 200 contaminant records detected across several MS platforms. The majority of records include theoretical as well as experimental MS data.",
                                                                    image_id = "maconda_logo"),
                                                     lipidmaps = list(title = "LIPID MAPS", 
                                                                    description = "The LIPID MAPS Structure Database (LMSD) is a relational database encompassing structures and annotations of biologically relevant lipids.",
                                                                    image_id = "lipidmaps_logo"),
                                                     # - - leave magicball last - - 
                                                     magicball = list(title = "MagicBall",
                                                                      description = "Algorithm to predict molecular formula from m/z value",
                                                                      image_id = "magicball"),
                                                     custom = list(title = "Custom",
                                                                   description = "Please select to start your own database creation!",
                                                                   image_id = "plus")
                                                     )
),
functions = list(# default color functions at startup, will be re-loaded from options
  cf = rainbow,
  color.function = rainbow,
  color.vec = rainbow,
  # available plot themes for ggplot. Can add more,also user-defined ones, 
  # but put them in shiny_general.R first so they are sourced properly.
  plot.themes = list(bw=ggplot2::theme_bw,
                     classic=ggplot2::theme_classic,
                     gray=ggplot2::theme_gray,
                     min=ggplot2::theme_minimal,
                     dark=ggplot2::theme_dark,
                     light=ggplot2::theme_light,
                     line=ggplot2::theme_linedraw),
  color.functions = { 
    # available colorbrewer themes to load into ggplot. These are the standard brew names used in their functions color.brewer etc.
    brew.cols <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", # - - sequential - -
                   "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                   "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", 
                   "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", # - - diverging - -
                   "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3" # - - qualitative - -
    )
    
    # generate direct functions from the brewer colours
    brew.opts <- lapply(brew.cols, function(opt) colorRampPalette(rev(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[opt,]$maxcolors, opt))))
    names(brew.opts) <- brew.cols
    
    # add the more general color scale functions as options
    base.opts <- list("rb"=rainbow,
                      "y2b"=ygobb,
                      "ml1"=matlab.like2,
                      "ml2"=matlab.like,
                      "m2g"=magenta2green,
                      "c2y"=cyan2yellow,
                      "b2y"=blue2yellow,
                      "g2r"=green2red,
                      "b2g"=blue2green,
                      "b2r"=blue2red,
                      "b2p"=cm.colors,
                      "bgy"=topo.colors,
                      "gyw"=terrain.colors,
                      "ryw"=heat.colors,
                      "bw"=blackwhite.colors)
    
    # add into a single list for use in interface
    append(base.opts, brew.opts)
  }),
# set default paths
paths = list(# path to .db file selected as default data source in the csv making pane
  patdb = file.path(options$work_dir, paste0(options$proj_name, ".db")),
  # path to .csv file selected as default source in the normalization pane
  csv_loc = file.path(options$work_dir, paste0(options$proj_name, ".csv")),
  # available paths when selecting a new file or folder
  volumes =  c('MetaboShiny' = getwd(),
               'Home'=home,
               '~' = normalizePath("~"),
               'Downloads'=file.path(home, "Downloads"),
               'R Installation'=R.home(),
               'Desktop'=file.path(home, "Desktop"))
),
# empty list to store result tables in at the statistics pane
tables = list(tbl = NA),
# default vectors to go through in metaboshiny
vectors = list(
  remove_match_cols = c("description","structure", "baseformula", "dppm", "source"),
  # default indices of chosen adducts
  pos_selected_add = c(1:3, nrow(adducts[Ion_mode == "positive",
                                         c("Name")])),
  neg_selected_add = c(1, 2, 14, 15, nrow(adducts[Ion_mode == "negative",
                                                  c("Name")])),
  # list of available databases!!
  db_list = c( # this dermines the show order of dbs in the app
              #"internal",
              #"noise", 
              "hmdb", 
              "chebi", 
              "kegg",
              "metacyc",
              "wikipathways",
              "smpdb",
              "dimedb",
              "wikidata",
              "vmh",
              "respect",
              "massbank",
              "metabolights",
              "foodb",
              "maconda",
              "bloodexposome",
              "lipidmaps",
              "magicball",
              "custom"
  ),
  # list of positive adducts
  pos_adducts = adducts[Ion_mode == "positive",
                        c("Name")],
  # list of negative adducts
  neg_adducts = adducts[Ion_mode == "negative",
                        c("Name")],
  wordcloud = list(top = 20)
)
)

global$vectors$wordcloud$skip <- unique(c( # manual curation(
                                  "on", "in", "and", "at", 
                                   "an", "by", "is", "it", "that", 
                                   "as", "be", "like", "can", "a", "of",
                                   "to", "but", "not", "mainly", "the",
                                   "", "which", "from",
                                   "found", "its", "two", "one", "if", "no",
                                   "yes", "any", "were", "observed", "also",
                                   "why", "other", "only", "known", "so", "do",
                                   "with", "resulting", "reaction", "via",
                                  "mtblc", "group", "groups",
                                   "metabolism", "mh", "ms", "position", "positions",
                                   "produced", "this", "ce", "mtblc",
                                   "has", "ko", "predicted", "are", "been", "isolated",
                                   "occurs", "form", "generated", "obtained", "et", "insource",
                                   "or", "nu", "substituted", "exhibits", "for", "ev",
                                   "attached", "constituent", "negative", "ph", "pmid",
                                   "compound", "residue", "unknown", "residues", "through",
                                   "lcesiitft", "lcesitof","lcesiqtof","have", "derived",
                                   "compounds", "having", "lcesiqq", "different", "more",
                                   "ec", "activity", "metabolite", "biotransformer","biotransformer¹",
                                  # pubmed words https://www.ncbi.nlm.nih.gov/books/NBK3827/table/pubmedhelp.T.stopwords/
                                  c("a", "about", "again", "all", "almost", "also", "although", 
                                    "always", "among", "an", "and", "another", "any", "are", "as", 
                                    "at", "be", "because", "been", "before", "being", "between", 
                                    "both", "but", "by", "can", "could", "did", "do", "does", "done", 
                                    "due", "during", "each", "either", "enough", "especially", "etc", 
                                    "for", "found", "from", "further", "had", "has", "have", "having", 
                                    "here", "how", "however", "i", "if", "in", "into", "is", "it", 
                                    "its", "itself", "just", "kg", "km", "made", "mainly", "make", 
                                    "may", "mg", "might", "ml", "mm", "most", "mostly", "must", "nearly", 
                                    "neither", "no", "nor", "obtained", "of", "often", "on", "our", 
                                    "overall", "perhaps", "pmid", "quite", "rather", "really", "regarding", 
                                    "seem", "seen", "several", "should", "show", "showed", "shown", 
                                    "shows", "significantly", "since", "so", "some", "such", "than", 
                                    "that", "the", "their", "theirs", "them", "then", "there", "therefore", 
                                    "these", "they", "this", "those", "through", "thus", "to", "upon", 
                                    "use", "used", "using", "various", "very", "was", "we", "were", 
                                    "what", "when", "which", "while", "with", "within", "without", 
                                    "would", LETTERS, letters, "acid", "cell", "cells", "human", "humans"),
                                  qdapDictionaries::Top200Words,
                                   global$vectors$db_list))

# load in custom databases
has.customs <- dir.exists(file.path(options$db_dir, "custom"))

if(has.customs){
  
  print("loading custom dbs...")
  
  customs = list.files(path = file.path(options$db_dir, "custom"), 
                       pattern = "\\.RData")
  
  dbnames = unique(tools::file_path_sans_ext(customs))
  
  for(db in dbnames){
    print(db)
    # add name to global
    dblist <- global$vectors$db_list
    dblist <- dblist[-which(dblist == "custom")]
    dblist <- c(dblist, db, "custom")
    global$vectors$db_list <- dblist
    
    metadata.path <- file.path(getOptions("user_options.txt")$db_dir, "custom", paste0(db, ".RData"))
    load(metadata.path)
    
    # add description to global
    global$constants$db.build.info[[db]] <- meta.dbpage
    
    # add image to global
    maxi = length(global$constants$images)
    global$constants$images[[maxi + 1]] <- meta.img 
  }
}

# split <- strsplit(str, split = ",")[[1]]
# split_words <- gsub(x = split, pattern = " ", replacement = "")
#' Gets the current used operating system. Important for parallel/multithreaded functions if using makeCluster("FORK")
 
#' \code{get_os} finds the name of the OS the user is running this function on.
#' 
#' @return osx, windows/win or linux
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
# === LOAD LIBRARIES ===

# load isotopes of all atoms, necessary for the build.extended.db function.
data(isotopes, package = "enviPat")

# create parallel workers, leaving 1 core for general use
# TODO: make this a user slider
session_cl <- parallel::makeCluster(max(c(1, parallel::detectCores()-1))) # leave 1 core for general use and 1 core for shiny session

# source the miniscript for the toggle buttons used in the interface (needs custom CSS)
source("./Rsource/SwitchButton.R")

#' Squishes HTML elements close together.
#' 
#' \code{sardine} used on consequtive html objects will make them sit next to each other neatly.
#'
#' @param content TagList or shiny function that generates HTML (such as an image..)
#' @return 'sardined' content that will sit close to its neighbors.
sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)

# interleave for sorting later ...
add_idx <- order(c(seq_along(global$vectors$pos_adducts$Name), seq_along(global$vectors$neg_adducts$Name)))
sort_order <<- unlist(c(global$vectors$pos_adducts$Name, global$vectors$neg_adducts$Name))[add_idx]

# generate CSS for the interface based on user settings for colours, fonts etc.
bar.css <<- nav.bar.css(options$col1, options$col2, options$col3, options$col4)
font.css <<- font.css(options$font1, options$font2, options$font3, options$font4,
                      options$size1, options$size2, options$size3, options$size4)

# google fonts

# name = the name of the font in Google's Library (https://fonts.google.com)
# family = how you want to refer to the font inside R
# regular.wt = the weight of font used in axis, labels etc.
# bolt.wt = the weight of font used in the title

# === GOOGLE FONT SUPPORT FOR GGPLOT2 ===

# Download a webfont
lapply(c(options[grepl(pattern = "font", names(options))]), function(font){
  sysfonts::font_add_google(name = font, family = font, regular.wt = 400, bold.wt = 700)
})

# Perhaps the only tricky bit is remembering to run the following function to enable webfonts
showtext::showtext_auto(enable = T)

# ======================================

# set default plot theme to minimal
plot.theme <<- ggplot2::theme_minimal

# set taskbar image as set in options
taskbar_image <<- options$task_img

# parse color options
global$vectors$mycols <- get.col.map("user_options.txt") # colours for discrete sets, like group A vs group B etc.
global$constants$spectrum <- options$gspec # gradient function for heatmaps, volcano plot etc.
global$vectors$project_names <- unique(tools::file_path_sans_ext(list.files(options$work_dir, pattern=".csv|.db"))) # the names listed in the 'choose project' tab of options.

# load existing file
fn <- paste0(tools::file_path_sans_ext(global$paths$patdb), ".metshi")
#

# if(exists("mSet")){
#   print("TESTING .. resetting analyses")
#   mSet$analSet <- NULL
# }

print("loaded scripts!")

if(file.exists(fn) & !exists("mSet")){
  print("loading existing mset ..")
  load(fn)
  print(".. done!")
}