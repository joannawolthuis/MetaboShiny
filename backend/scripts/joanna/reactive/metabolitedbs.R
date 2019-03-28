# everything below uses the dblist defined in global
# as well as the logos defined here
# if you add a db, both the name and associated logo need to be added

# create checkcmarks if database is present
lapply(global$vectors$db_list, FUN=function(db){
  # creates listener for if the 'check db' button is pressed
  observeEvent(input[[paste0("check_", db)]],{
    # see which db files are present in folder
    db_folder_files <- list.files(getOptions()$db_dir)
    is.present <- paste0(db, ".base.db") %in% db_folder_files
    check_pic <- if(is.present) "yes.png" else "no.png"
    # generate checkmark image objects
    output[[paste0(db,"_check")]] <- renderImage({
      filename <- normalizePath(file.path('www', check_pic))
      list(src = filename, width = 70,
           height = 70)
    }, deleteFile = FALSE)
  })
})

# these listeners trigger when build_'db' is clicked (loops through dblist in global)
lapply(global$vectors$db_list, FUN=function(db){
  observeEvent(input[[paste0("build_", db)]], {
    # ---------------------------
    library(RCurl)
    library(XML)
    library(SPARQL)
    # ---------------------------
    withProgress({

      # send necessary functions and libraries to parallel threads
      parallel::clusterEvalQ(session_cl, library(enviPat))
      parallel::clusterEvalQ(session_cl, library(KEGGREST))
      parallel::clusterEvalQ(session_cl, library(RCurl))
      parallel::clusterEvalQ(session_cl, library(SPARQL))
      parallel::clusterEvalQ(session_cl, library(XML))
      parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
        "isotopes",
        "subform.joanna",
        "mergeform.joanna",
        "multiform.joanna",
        "check.ded.joanna",
        "kegg.charge",
        "xmlParse",
        "getURL",
        "mape",
        "flattenlist"
      ))
      pkgs = c("data.table", "enviPat", "KEGGREST", "XML", "SPARQL", "RCurl")
      parallel::clusterCall(session_cl, function(pkgs) {
        for (req in pkgs) {
          require(req, character.only = TRUE)
        }
      }, pkgs = pkgs)
      shiny::setProgress(session = session, 0.1)

      # build base db (differs per db, parsers for downloaded data)
      build.base.db(db,
                    outfolder = getOptions()$db_dir,
                    cl = session_cl)
      shiny::setProgress(session = session, 0.5)

      if(!grepl(db, pattern = "maconda|noise")){
        # extend base db (identical per db, makes adduct and isotope variants of downloaded compounds)
        build.extended.db(db,
                          outfolder = getOptions()$db_dir,
                          adduct.table = adducts,
                          cl = F,#session_cl,
                          fetch.limit = 500) #TODO: figure out the optimal fetch limit...
      }
    })
  })
})
