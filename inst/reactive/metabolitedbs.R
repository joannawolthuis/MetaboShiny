# everything below uses the dblist defined in global
# as well as the logos defined here
# if you add a db, both the name and associated logo need to be added

# create checkcmarks if database is present
lapply(gbl$vectors$db_list, FUN=function(db){
  # creates listener for if the 'check db' button is pressed
  shiny::observeEvent(input[[paste0("check_", db)]],{
    # see which db files are present in folder
    db_folder_files <- list.files(lcl$paths$db_dir, full.names = T)
    dbname = paste0(db, ".db")
    is.present <- dbname %in% basename(db_folder_files)
    if(is.present){
      conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(file.path(lcl$paths$db_dir, dbname))) # change this to proper var later
      is.present <- RSQLite::dbExistsTable(conn, "base")
      RSQLite::dbDisconnect(conn)
    }
    check_pic <- if(is.present) "yes.png" else "no.png"
    # generate checkmark image objects
    output[[paste0(db,"_check")]] <- renderImage({
      filename <- normalizePath(file.path(getwd(), 'www', check_pic))
      list(src = filename, width = 70,
           height = 70)
    }, deleteFile = FALSE)
  })
})

shiny::observeEvent(input$build_custom_db, {

  cust_dir = file.path(lcl$paths$db_dir, paste0(input$my_db_short, "_source"))
  # make folder for this db
  if(dir.exists(cust_dir)) unlink(cust_dir)
  dir.create(cust_dir)
  
  # copy csv and imageto said folder
  img_path <- shinyFiles::parseFilePaths(gbl$paths$volumes, input$custom_db_img_path)$datapath
  file.copy(img_path, file.path(cust_dir, "logo.png"))
  
  csv_path <- shinyFiles::parseFilePaths(gbl$paths$volumes, input$custom_db)$datapath
  file.copy(csv_path, file.path(cust_dir, "base.csv"))

  dbinfo = list(title = input$my_db_name,
                description = input$my_db_description,
                image_id = paste0(input$my_db_short, "_logo"))
  
  save(dbinfo, file = file.path(cust_dir, "info.RData"))
  # print OK message and ask to restart
  print("Import OK! Please restart MetaboShiny to view and build your database :-)")
  
  shiny::removeModal()
  
})

# these listeners trigger when build_'db' is clicked (loops through dblist in global)
lapply(c(gbl$vectors$db_list), FUN=function(db){
  shiny::observeEvent(input[[paste0("build_", db)]], {
    withProgress({
      # send necessary functions and libraries to parallel threads
      parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
        "isotopes"
      )) 
      pkgs = c("data.table", "enviPat", 
               "KEGGREST", "XML", 
               "SPARQL", "RCurl", 
               "MetaDBparse")
      parallel::clusterCall(session_cl, function(pkgs) {
        for (req in pkgs) {
          library(req, character.only = TRUE)
        }
      }, pkgs = pkgs)
      
      shiny::setProgress(session = session, 0.1)

      if(input$db_build_mode %in% c("base", "both")){
        # check if custom
        if(!(db %in% gbl$vectors$db_list)){
          custom = TRUE
        }
        # - - - - - - - -
        MetaDBparse::buildBaseDB(dbname = db,
                                 outfolder = normalizePath(lcl$paths$db_dir), 
                                 cl = session_cl,
                                 custom_csv_path = if(!custom) NULL else file.path(lcl$paths$db_dir, paste0(db,"_source"), "base.csv"),
                                 silent = F)
      }
      
      # build base db (differs per db, parsers for downloaded data)
      shiny::setProgress(session = session, 0.5)

      if(input$db_build_mode %in% c("extended", "both")){

      if(!grepl(db, pattern = "maconda")){
        if(file.exists(file.path(lcl$paths$db_dir, paste0(db, ".db")))){
          my_range <- input$db_mz_range
          outfolder <- lcl$paths$db_dir
          MetaDBparse::buildExtDB(base.dbname = db,
                                  outfolder = outfolder,
                                  cl = session_cl,
                                  blocksize = 500,
                                  mzrange = my_range,
                                  adduct_table = adducts,
                                  adduct_rules = adduct_rules, 
                                  silent = T,
                                  ext.dbname = "extended") #TODO: figure out the optimal fetch limit... seems 200 for now
        }else{
          print("Please build base DB first! > _<")
        }
      }
        } 
    })
  })
})
