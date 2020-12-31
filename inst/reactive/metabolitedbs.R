# everything below uses the dblist defined in global
# as well as the logos defined here
# if you add a db, both the name and associated logo need to be added

# ==== DATABASES UI ====

shiny::observe({
  if(dbmanager$build[1] != "none"){
    for(db in dbmanager$build){
      #shiny::withProgress({
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
        
        if(input$db_build_mode %in% c("base", "both")){
          # check if custom
          custom_csv = file.path(lcl$paths$db_dir, paste0(db,"_source"), "base.csv")
          custom = file.exists(custom_csv)
          # - - - - - - - -
          MetaDBparse::buildBaseDB(dbname = db,
                                   outfolder = normalizePath(lcl$paths$db_dir), 
                                   cl = 0,#session_cl,
                                   custom_csv_path = if(!custom) NULL else custom_csv,
                                   silent = F)
        }
        
        if(input$db_build_mode %in% c("extended", "both")){
          if(!grepl(db, pattern = "maconda")){
            if(file.exists(file.path(lcl$paths$db_dir, paste0(db, ".db")))){
              my_range <- input$db_mz_range
              outfolder <- lcl$paths$db_dir
              all.isos <- input$db_all_iso
              count.isos <- input$db_count_iso
              MetaDBparse::buildExtDB(base.dbname = db,
                                      outfolder = outfolder,
                                      cl = session_cl,
                                      blocksize = 500,
                                      mzrange = my_range,
                                      adduct_table = adducts,
                                      adduct_rules = adduct_rules, 
                                      silent = T,
                                      all.isos = all.isos,
                                      count.isos = count.isos,
                                      ext.dbname = "extended") #TODO: figure out the optimal fetch limit... seems 200 for now
            }else{
              MetaboShiny::metshiAlert("Please build base DB first! (can be changed in settings)")
            }
          }
        } 
      #})
    }
    dbmanager$build <- "none"
  }  
})

shiny::observeEvent(input$db_build_sel_all, {
  for(db in gbl$vectors$db_list){
    input.id = paste0("build_queue_", db)
    if(input.id %in% names(input)){
      shinyWidgets::updatePrettyToggle(session, 
                                       inputId = input.id,
                                       value = if(input[[input.id]]) F else T)    
    }
  }
})

shiny::observeEvent(input$db_build_multi_all, {
  build.me <- unlist(sapply(gbl$vectors$db_list, function(db){
    input.id = paste0("build_queue_", db)
    if(input.id %in% names(input)){
      if(input[[input.id]]) db else NULL
    }else{
      NULL
    }
  }))
  if(length(build.me) > 0){
    dbmanager$build <- build.me
  }
})

shiny::observe({
  if(db_section$load){
    shiny::showNotification("Loading database screen...")
    # - - - load version numbers - - - 
    db.paths = list.files(lcl$paths$db_dir, pattern = "\\.db$",full.names = T)
    versions = lapply(db.paths,
           function(path){
             ver = "unknown"
             date = "unknown"
             try({
               conn <- RSQLite::dbConnect(RSQLite::SQLite(), path) # change this to proper var later
               meta = RSQLite::dbGetQuery(conn, "SELECT * FROM metadata")
               ver = meta$version[[1]]
               suppressWarnings({
                 numver = as.numeric(ver)
               })
               if(!is.na(numver)){
                 if(numver > 18000){
                   ver = as.character(as.Date(ver, origin = "1970-01-01"))
                 }  
               }
               date = meta$date[[1]]
               date = as.character(as.Date(date, origin = "1970-01-01"))
               RSQLite::dbDisconnect(conn)  
             },silent = T)
             if(is.na(ver)) ver <- "unknown"
             dbname = gsub(basename(path), 
                           pattern = "\\.db$", 
                           replacement="")
             if(ver != "unknown" & date != "unknown"){
               if(ver != date){
                 ver = as.character(gsubfn::fn$paste("Version $ver downloaded on $date")) 
               }else{
                 ver = as.character(gsubfn::fn$paste("Downloaded on $date")) 
               }
               output[[paste0(dbname, "_version")]] <- renderText({ver})
             }else{
               output[[paste0(dbname, "_version")]] <- renderText({""})
              }
             ver
           })
    names(versions) <- gsub(basename(db.paths), 
                            pattern = "\\.db$", 
                            replacement="")
    
    lcl$vectors$db.version <<- versions
                    
    lapply(c("db", "db_prematch"), function(midfix){
      shiny::observeEvent(input[[paste0("select_", midfix, "_all")]], {
        
        if(length(lcl$vectors$built_dbs) == 0){
          MetaboShiny::metshiAlert("Please create at least one database to use this feature!")
          NULL
        }else{
          dbs <- lcl$vectors$built_dbs[-which(lcl$vectors$built_dbs %in% gbl$vectors$db_no_build)]
          
          currently.on <- sapply(dbs, function(db){
            input[[paste0(switch(midfix, 
                                 "db" = "search_",
                                 "db_prematch" = "prematch_"), db)]]
          })
          
          if(any(unlist(currently.on))){
            set.to = F
          }else{
            set.to = T
          }
          
          for(db in dbs){
            shiny::updateCheckboxInput(session, paste0(switch(midfix, 
                                                              "db" = "search_",
                                                              "db_prematch" = "prematch_"), db), value = set.to)
          } 
        }
      })      
    })
    
    # check for the magicball-requiring databases...
    
    # render the database download area
    output$db_build_ui <- renderUI({
      dbs_per_line = 4
      max_col_width = 12
      rows = ceiling(length(gbl$vectors$db_list) / dbs_per_line)
      database_layout = lapply(1:rows, function(i){
        min_i = (dbs_per_line * i) - (dbs_per_line - 1)
        max_i = (dbs_per_line * i)
        if(max_i > length(gbl$vectors$db_list)) max_i <- length(gbl$vectors$db_list)
        # create 3 fluidrows followed by a break
        list(
          # row 1: name
          shiny::fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
            shiny::column(width=3,align="center", h2(gbl$constants$db.build.info[[db]]$title))
          })),
          # row 2: description
          shiny::fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
            shiny::column(width=3,align="center", shiny::helpText(gbl$constants$db.build.info[[db]]$description))
          })),
          # row 3: image
          shiny::fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
            if(db != "custom"){
              shiny::column(width=3,align="center", 
                            shiny::imageOutput(gbl$constants$db.build.info[[db]]$image_id, inline=T),
                            br(),br(),
                            shiny::div(shiny::tags$i(shiny::textOutput(paste0(db, "_version"))),style='font-size:70%; color: grey')
                            ,br()
                            )
            }else{
              shiny::column(width=3,align="center", shinyWidgets::circleButton("make_custom_db",
                                                                               size = "lg",
                                                                               icon = shiny::icon("plus",class = "fa-lg"),
                                                                               style = "stretch",
                                                                               color = "default"))
            }
          })),
          # row 4: button
          shiny::fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
            shiny::column(width=3, align="center",
                          if(!(db %in% gbl$vectors$db_no_build)){
                            list(
                              shinyBS::tipify(shiny::actionLink(paste0("check_", db),
                                                label = "",
                                                icon = icon("check")),
                                              title = "is base database built?"),
                              MetaboShiny::sardine(shiny::conditionalPanel("input.db_build_multi == false", 
                                                      shinyBS::tipify(shiny::actionLink(paste0("build_", db),
                                                                                        label = "",
                                                                                        icon = shiny::icon("wrench")),
                                                                      title = "build this database")
                                                      )),
                              MetaboShiny::sardine(shiny::conditionalPanel("input.db_build_multi == true", 
                                                      shinyWidgets::prettyToggle(
                                                        status_off = "default", 
                                                        status_on = "success",
                                                        inline=T,bigger=F,
                                                        animation="pulse",
                                                        inputId = paste0("build_queue_", db),
                                                        label_on = "", 
                                                        label_off = "",
                                                        outline = TRUE,
                                                        plain = TRUE,
                                                        value = db %in% gbl$vectors$db_categories$favorite,
                                                        icon_on = shiny::icon("wrench",lib ="glyphicon"), 
                                                        icon_off = shiny::icon("unchecked",lib ="glyphicon")
                                                      ))),
                              shinyBS::tipify(shinyWidgets::prettyToggle(
                                status_off = "default", 
                                status_on = "danger",
                                inline=T,bigger=F,
                                animation="pulse",
                                inputId = paste0("favorite_", db),
                                label_on = "", 
                                label_off = "",
                                outline = TRUE,
                                plain = TRUE,
                                value = db %in% gbl$vectors$db_categories$favorite,
                                icon_on = icon("heart",lib ="glyphicon"), 
                                icon_off = icon("heart-empty",lib ="glyphicon")
                              ), title = "add this database to favorites category"),
                              shiny::br(),shiny::br(),
                              shiny::imageOutput(paste0(db, "_check"),inline = T))
                          }else{
                            list()
                          }
            )
          })),
          shiny::br())
      })
      # return
      database_layout
    })
    
    db_button_prefixes = c("search", "prematch")
    
    # generate all the fadebuttons for the database selection
    lapply(db_button_prefixes, function(prefix){
      
      output[[paste0("db_", prefix, "_select")]] <- renderUI({
        db.paths = list.files(lcl$paths$db_dir, pattern = "\\.db$",full.names = T)
        built.dbs <- c(gsub(x = basename(db.paths), 
                            pattern = "\\.db", replacement = ""), 
                       gbl$vectors$db_no_build)
        really.built.dbs <- sapply(db.paths, function(path) {
          conn <- RSQLite::dbConnect(RSQLite::SQLite(), path) # change this to proper var later
          exists = RSQLite::dbExistsTable(conn, "base")
          if(exists) exists = RSQLite::dbGetQuery(conn, "select count(*) from base")[1,] > 0 
          RSQLite::dbDisconnect(conn)
          exists
        })
        really.built.dbs <- db.paths[really.built.dbs]
        really.built.dbs <- gsub(x = basename(really.built.dbs), 
                                 pattern = "\\.db", replacement = "")
        
        no.need.build = c("cmmmediator", "pubchem","chemspider","supernatural2","knapsack","chemidplus", "magicball")
        if(length(really.built.dbs) > 0){
          built.dbs <- unique(c(no.need.build,
                                intersect(really.built.dbs,
                                          gbl$vectors$db_list)))
        }else{
          built.dbs <- list(no.need.build)
        }
        
        lcl$vectors$built_dbs <<- built.dbs
        
        if(length(lcl$vectors$built_dbs) == 0){
          MetaboShiny::metshiAlert("Please create at least one database to use this feature!")
          shiny::fluidRow(align="center", 
                          br(),
                          helpText("No databases built..."),
                          br())
        }else{
          iconPicks = list(
            all = "cart-plus",
            versatile = "map-signs",
            verbose = "book",
            livestock = "piggy-bank",
            human = "male",
            microbial = "splotch",
            pathway = "road",
            food = "utensil-spoon",
            plant = "seedling",
            massspec = "fingerprint",
            chemical = "flask",
            online = "globe",
            study = "scroll",
            predictive = "magic",
            custom = "cart-plus",
            favorite = "heart")  
          
          iconWrap <- sapply(iconPicks, function(ic){
            gsubfn::fn$paste("<i class='fa fa-$ic'></i>")
          })
          
          choices = names(iconPicks)
          names(choices) <- iconWrap
          
          tooltips = lapply(as.character(choices), function(choice){
            radioTooltip(id = paste0(prefix, "_db_categories"),
                         choice = choice,
                         title = paste(choice, "databases"),
                         trigger = "hover",
                         placement="right")
          })
          
          list(
            shiny::fluidRow(align="center",
                            shinyWidgets::checkboxGroupButtons(
                              inputId = paste0(prefix, "_db_categories"),
                              label = "", 
                              choices = choices, selected = "all",
                              justified = TRUE,size = "sm"
                            )
                            
            ),
            tooltips,
            shiny::wellPanel(id = "def",
                             style = "overflow-y:scroll; max-height: 250px; border:1px dashed #e3e3e3; background-color: #ffffff;",
                             shiny::uiOutput(paste0(prefix,"_db_categ")))
          )
        }
      })
    })
    
    lapply(db_button_prefixes, function(prefix){
      shiny::observeEvent(input[[paste0(prefix,"_db_categories")]], {
        output[[paste0(prefix,"_db_categ")]] <- shiny::renderUI({
          
          considered_all = gbl$vectors$db_list[which(gbl$vectors$db_list != "custom" & gbl$vectors$db_list %in% lcl$vectors$built_dbs)]
          
          lapply(considered_all, function(db){
            tag = paste0(prefix, "_", db)
            shinyjs::runjs('Shiny.onInputChange("$tag, null)')
          })
          
          dbs_categ <- intersect(considered_all, unlist(gbl$vectors$db_categories[input[[paste0(prefix,"_db_categories")]]]))
          display = intersect(dbs_categ, considered_all)
          shiny::fluidRow(
            lapply(display, function(db){
              which_idx = grep(sapply(gbl$constants$images, function(x) x$name), pattern = db) # find the matching image (NAME MUST HAVE DB NAME IN IT COMPLETELY)
              shinyBS::tipify(shiny::div(style="display: inline-block;vertical-align:top;",
                                         MetaboShiny::fadeImageButton(inputId = paste0(prefix, "_", db), 
                                                                      img.path = gbl$constants$images[[which_idx]]$path)),
                              title = gbl$constants$db.build.info[[db]]$title) # generate fitting html
            })
          )
        })
      })
    })
    
    # check if these buttons are selected or notr
    lapply(db_button_prefixes, function(prefix){
      shiny::observe({
        # ---------------------------------
        db_path_list <- lapply(gbl$vectors$db_list, # go through the dbs defined in db_lists
                               FUN = function(db){
                                 button_id = input[[paste0(prefix, "_", db)]]
                                 if(is.null(button_id)){
                                   NA
                                 }else{
                                   if(!button_id){
                                     c(db)# add path to list of dbpaths
                                   }
                                   else{NA}
                                 }
                               }
        )
        # save the selected database paths to global
        lcl$vectors[[paste0("db_", prefix, "_list")]] <<- db_path_list[!is.na(db_path_list)]
      })
    })
    
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
          filename <- normalizePath(file.path('www', check_pic))
          list(src = filename, width = 70,
               height = 70)
        }, deleteFile = FALSE)
      })
    })
    # these listeners trigger when build_'db' is clicked (loops through dblist in global)
    lapply(c(gbl$vectors$db_list), FUN=function(db){
      shiny::observeEvent(input[[paste0("favorite_", db)]], {
        favorites = names(which(unlist(sapply(gbl$vectors$db_list, function(db) input[[paste0("favorite_", db)]]))))
        if(!is.null(favorites)){
          if(length(favorites) > 0 ){
            MetaboShiny::setOption(lcl$paths$opt.loc, "dbfavs", paste0(favorites, collapse=","))
            gbl$vectors$db_categories$favorite <<- favorites     
          }
        }
      })
    })
    
    # these listeners trigger when build_'db' is clicked (loops through dblist in global)
    lapply(c(gbl$vectors$db_list), FUN=function(db){
      shiny::observeEvent(input[[paste0("build_", db)]], {
        dbmanager$build <- db
      })
    })
    modifyStyle("body", background = "white")
    shiny::removeModal()
  }
})

shiny::observeEvent(input$build_custom_db, {
  
  cust_dir = file.path(lcl$paths$db_dir, paste0(input$my_db_short, 
                                                "_source"))
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
  shiny::showNotification("Import OK! Please restart MetaboShiny to view and build your database.")
  
  shiny::removeModal()
  
  #Sys.sleep(5)
  #setHeartLoader(40)
})



output$db_example <- DT::renderDataTable({
  DT::datatable(data = data.table::data.table(
    compoundname = c("1-Methylhistidine", "1,3-Diaminopropane", "2-Ketobutyric acid"),
    description = c("One-methylhistidine (1-MHis) is derived mainly from the anserine of dietary flesh sources, especially poultry.",
                    "1,3-Diaminopropane is a stable, flammable and highly hydroscopic fluid. It is a polyamine that is normally quite toxic if swallowed, inhaled or absorbed through the skin.",
                    "2-Ketobutyric acid is a substance that is involved in the metabolism of many amino acids (glycine, methionine, valine, leucine, serine, threonine, isoleucine) as well as propanoate metabolism and C-5 branched dibasic acid metabolism. "),
    baseformula = c("C7H11N3O2", "C3H10N2", "C4H6O3"),
    identifier = c("HMDB1", "HMDB2", "HMDB3"),
    charge = c(0, 0, 0),
    structure = c("CN1C=NC(C[C@H](N)C(O)=O)=C1", "NCCCN", "CCC(=O)C(O)=O")
  ),
  options = list(searching = FALSE,
                 paging = FALSE,
                 info = FALSE))
})

# observeEvent
shiny::observeEvent(input$make_custom_db, {
  
  # get window
  showModal(modalDialog(
    shiny::fluidRow(align="center",
                    shiny::textInput("my_db_name", label = "Database full name", value = "MyDb"),
                    shiny::textInput("my_db_short", label = "Database shorthand name", value = "mydb"),
                    shiny::textInput("my_db_description", label = "Database description", value = "Custom database for MetaboShiny."),
                    shiny::hr(),
                    shiny::helpText("Please input a CSV file with at these columns (example below):"),
                    shiny::helpText("'baseformula', 'charge', 'compoundname', 'identifier', 'description', 'structure' (in SMILES!)"),
                    shiny::div(DT::dataTableOutput("db_example"), style="font-size:60%"),
                    shiny::br(),
                    shinyFiles::shinyFilesButton("custom_db", title = "Please pick a .csv file", multiple = F, label = "Select"),
                    shiny::hr(),
                    shiny::helpText("Please upload a database logo"),
                    shinyFiles::shinyFilesButton("custom_db_img_path",
                                                 'Select image',
                                                 'Please select a .png file',
                                                 FALSE),shiny::br(),
                    shiny::imageOutput("custom_db_img", inline=T),shiny::br(),
                    shiny::hr(),
                    shinyWidgets::circleButton("build_custom_db",icon = shiny::icon("arrow-right", "fa-lg"), size = "lg")
    ),
    title = "Custom database creation",
    easyClose = TRUE,
    size = "l",
    footer = NULL
  )
  )
})
