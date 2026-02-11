read_images_csv <- function(path = system.file("extdata", "images.csv", package = "MetaboShiny")){
  if(is.null(path) || path == ""){
    path <- file.path("inst", "extdata", "images.csv")
  }
  img <- utils::read.csv(path, stringsAsFactors = FALSE)
  parse_dim <- function(x){
    if(grepl("^[0-9.]+$", x)) as.numeric(x) else x
  }
  lapply(seq_len(nrow(img)), function(i){
    list(
      name = img$name[i],
      path = img$path[i],
      dimensions = c(parse_dim(img$dim1[i]), parse_dim(img$dim2[i]))
    )
  })
}

read_db_build_info_csv <- function(path = system.file("extdata", "db_build_info.csv", package = "MetaboShiny")){
  if(is.null(path) || path == ""){
    path <- file.path("inst", "extdata", "db_build_info.csv")
  }
  db <- utils::read.csv(path, stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(db)), function(i){
    list(title = db$title[i],
         description = db$description[i],
         image_id = db$image_id[i])
  })
  names(res) <- db$key
  res
}

build_gbl_constants <- function(){
  list(
    ppm = 2, # TODO: re-add ppm as option for people importing their data through csv
    ml.twoonly = c("adaboost","logicBag","bartMachine","binda",
                   "ada","gamboost","glmboost","chaid",
                   "C5.0Cost","rpartCost","deepboost",
                   "dwdPoly","dwdRadial","glm","glmnet",
                   "glmStepAIC","glmnet_h2o","svmLinearWeights2",
                   "dwdLinear","svmLinearWeights","logreg","mlpKerasDropoutCost",
                   "mlpKerasDecayCost","ORFlog","ORFpls","ORFridge","ORFsvm",
                   "plsRglm","rotationForest","rotationForestCp",
                   "svmRadialWeights","nodeHarvest"),
    # get all caret models that can do classification and have some kind of importance metric
    ml.models = {
      caret.mdls <- caret::getModelInfo()
      fin = names(caret.mdls)[sapply(seq_along(caret.mdls), function(i){
        curr.mdl = caret.mdls[[i]]
        can.classify = if("Classification" %in% curr.mdl$type) TRUE else FALSE
        #has.importance = if("varImp" %in% names(curr.mdl)) TRUE else FALSE
        can.classify# & has.importance
      })]
      fin = c(fin, "glm (logistic)")
      caret.mdls <- NULL
      fin
    },
    max.cols = 30,
    images = read_images_csv(),
    default.text = list(list(name='curr_definition', text="No m/z selected"),
                        list(name="curr_cpd", text="..."),# default text options at startup
                        list(name="ml_train_ss", text="all"),
                        list(name="ml_test_ss", text="all")
    ),
    db.build.info = read_db_build_info_csv()
  )
}

build_gbl_functions <- function(){
  list(
    # default color functions at startup, will be re-loaded from options
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
      brew.opts <- lapply(brew.cols, function(opt) grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[opt,]$maxcolors, opt))))
      names(brew.opts) <- brew.cols
      
      # add the more general color scale functions as options
      base.opts <- list("rb"=rainbow,
                        "y2b"=colorRamps::ygobb,
                        "ml1"=colorRamps::matlab.like2,
                        "ml2"=colorRamps::matlab.like,
                        "m2g"=colorRamps::magenta2green,
                        "c2y"=colorRamps::cyan2yellow,
                        "b2y"=colorRamps::blue2yellow,
                        "g2r"=colorRamps::green2red,
                        "b2g"=colorRamps::blue2green,
                        "b2r"=colorRamps::blue2red,
                        "b2p"=grDevices::cm.colors,
                        "bgy"=grDevices::topo.colors,
                        "gyw"=grDevices::terrain.colors,
                        "ryw"=grDevices::heat.colors,
                        "bw"=MetaboShiny::blackwhite.colors)
      
      # add into a single list for use in interface
      fin = append(base.opts, brew.opts)
      brew.cols <- brew.opts <- base.opts <- NULL
      fin
    }
  )
}

build_gbl_paths <- function(){
  list(
    # available paths when selecting a new file or folder
    volumes = {
      vols = c("Recent",
              "Your Files",
               "Home",
               "Documents",
               "Downloads",
               "Desktop",
               "Examples",
               "System Root")
      home = path.expand('~')
      folders = lapply(vols, 
                       FUN = function(folder){
                         switch(folder,
                                `System Root` = "/",
                                Recent = system.file("examples",
                                                     package = "MetaboShiny"),
                                Home = home,
                                `Your Files` = file.path(home, "MetaboShiny"),
                                Documents = {
                                  loc = file.path(home, "Documents")
                                  if(dir.exists(loc)) loc else NULL
                                },
                                Downloads = {
                                  loc = file.path(home, "Downloads")
                                  if(dir.exists(loc)) loc else NULL
                                },
                                Desktop = {
                                  loc = file.path(home, "Desktop")
                                  if(dir.exists(loc)) loc else NULL
                                },
                                Examples = {
                                  system.file("examples",
                                              package = "MetaboShiny")
                                })
                       })
      names(folders) = vols

      for(diskletter in LETTERS){
        if (Sys.info()["sysname"] == "Windows") {
          dn <- paste0(diskletter, ":/")
          if(dir.exists(dn)){
            print(paste0("Adding disk: ", dn, " to file picker."))
            folders[dn] <- dn
          }
        }else{
          #("to do")
          NULL
        }
      }
      unlist(folders[!sapply(folders, is.null)])
    }
  )
}

build_gbl_vectors <- function(adducts){
  list(
    hide_match_cols = c("structure", "identifier","baseformula",
                        "isocat", "fullformula", "finalcharge", "query_mz"),
    # list of available databases!!
    db_no_build = c("cmmmediator",
                    "chemspider",
                    "magicball",
                    "knapsack",
                    "chemidplus",
                    'supernatural2',
                    "custom",
                    "pubchem"),
    db_categories = list(versatile = c("wikidata", "dimedb", "metacyc", "chebi", "massbank", "cmmediator"),
                         verbose = c("hmdb", "chebi", "t3db", "metabolights", "ymdb", "ecmdb", "pamdb", "metabolomicsworkbench", "markerdb"),
                         livestock = c("lmdb", "bmdb", "metacyc", "mcdb"),
                         human = c("hmdb", "metacyc", "expoexplorer", "t3db", "bloodexposome", "pharmgkb", "markerdb"),
                         microbial = c("ymdb", "ecmdb", "pamdb", "vmh", "mvoc", "npa"),
                         pathway = c("vmh", "smpdb", "kegg", "reactome"),
                         food = c("foodb", "phenolexplorer"),
                         plant = c("anpdb", "respect", "metacyc", "supernatural2", "knapsack"),
                         chemical = c("chebi", "massbank", "maconda", "stoff", "lipidmaps", "chemspider"),
                         massspec = c("massbank", "respect", "maconda", "dimedb"),
                         online = c("cmmmediator", "chemspider","pubchem","knapsack","chemidplus","supernatural2"),
                         studies = c("metabolights","metabolomicsworkbench"),
                         custom = c(),
                         predictive = c("magicball", "chemspider", "pubchem", "supernatural2", "knapsack","chemidplus")
    ),
    db_list = c( # this determines the show order of dbs in the app
      "hmdb",
      "chebi",
      "kegg",
      "metacyc",
      "mcdb",
      "pamdb",
      "mvoc",
      #"dsstox",
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
      "expoexplorer",
      "lipidmaps",
      't3db',
      'drugbank',
      'lmdb',
      'ymdb',
      'ecmdb',
      #'rmdb',
      'bmdb',
      'stoff',
      'anpdb',
      'pharmgkb',
      'reactome',
      'metabolomicsworkbench',
      'phenolexplorer',
      'npa',
      'markerdb',
      "magicball",
      'cmmmediator',
      'pubchem',
      'chemidplus',
      'supernatural2',
      'chemspider',
      'knapsack',
      "custom"
    ),
    # list of positive adducts
    pos_adducts = adducts[Ion_mode == "positive",
                          c("Name")],
    # list of negative adducts
    neg_adducts = adducts[Ion_mode == "negative",
                          c("Name")],
    wordcloud = list(top = 20),
    calc_adducts = c("M+H", "M-H")
  )
}
