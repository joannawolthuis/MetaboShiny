# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

shinyServer(function(input, output, session) {

  # ================================= DEFAULTS ===================================

  source('./backend/scripts/joanna/shiny_general.R')

  # set progress bar style to 'old' (otherwise it's not movable with CSS)
  shinyOptions(progress.style="old")

  options(digits=22,
          spinner.size = 0.5,
          spinner.type = 6,
          spinner.color = "black",
          spinner.color.background = "white")

  # loading screen
  loadModal <- function(failed = FALSE) {
    modalDialog(
      fluidRow(align="center",
               helpText("Initializing MetaboShiny"),
               shinycssloaders::withSpinner(helpText(NULL))
      )
    )
  }

  showModal(loadModal())

  # send specific functions/packages to other threads
  parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
    "mape",
    "flattenlist"))

  parallel::clusterEvalQ(session_cl, library(data.table))

  # create default text objects in UI
  lapply(global$constants$default.text, FUN=function(default){
    output[[default$name]] = renderText(default$text)
  })

  library(showtext)

  # import google fonts
  for(font in unlist(options[grep(names(options), pattern = "font")])){
    if(font %in% sysfonts::font.families()){
      NULL
    }else{
      sysfonts::font_add_google(font,db_cache = T)
    }
  }

  showtext::showtext_auto() ## Automatically use showtext to render text for future devices

  # default match table fill
  output$match_tab <-DT::renderDataTable({
    if(is.null(global$tables$last_matches)){
      DT::datatable(data.table("no m/z chosen" = "Please choose m/z value from results ٩(｡•́‿•̀｡)۶	"),
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15),
                                   pageLength = 5))
    }
  })

  # create image objects in UI
  lapply(global$constants$images, FUN=function(image){
    output[[image$name]] <- renderImage({
      filename <- normalizePath(image$path)
      # Return a list containing the filename and alt text
      list(src = filename,
           width = image$dimensions[1],
           height = image$dimensions[2])
    }, deleteFile = FALSE)
  })

  # create color pickers based on amount of colours allowed in global
  output$colorPickers <- renderUI({
    lapply(c(1:global$constants$max.cols), function(i) {
      colourpicker::colourInput(inputId = paste("col", i, sep="_"),
                                label = paste("Choose colour", i),
                                value = global$vectors$mycols[i],
                                allowTransparent = F)
    })
  })

  # create color1, color2 etc variables to use in plotting functions
  # and update when colours picked change
  observe({
    values <- unlist(lapply(c(1:global$constants$max.cols), function(i) {
      input[[paste("col", i, sep="_")]]
    }))
    if(!any(is.null(values))){
      set.col.map(file.path(optfolder, 'user_options.txt'), values)
      global$vectors$mycols <<- get.col.map(file.path(optfolder, paste0('user_options_', runmode, ".txt")))
    }
  })

  shown_matches <- reactiveValues(table = global$tables$last_matches)

  output$match_tab <- DT::renderDataTable({
      remove_cols = global$vectors$remove_match_cols
      remove_idx <- which(colnames(global$tables$last_matches) %in% remove_cols)
      # don't show some columns but keep them in the original table, so they can be used
      # for showing molecule descriptions, structure
      DT::datatable(shown_matches$table,
      selection = 'single',
      autoHideNavigation = T,
      options = list(lengthMenu = c(5, 10, 15),
                     pageLength = 5,
                     columnDefs = list(list(visible=FALSE, targets=remove_idx)))
      )

  })
  # ===== UI SWITCHER ====

  # create interface mode storage object.
  interface <- reactiveValues()

  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  observe({

    # hide all tabs by default, easier to hide them and then make visible selectively

    hide.tabs <- list(
      list("statistics", "inf"),
      list("dimred", "pca"),
      list("dimred", "plsda"),
      list("permz", "asca"),
      list("permz", "meba"),
      list("permz", "aov"),
      list("statistics", "ml"),
      list("overview", "volc"),
      list("overview", "venn"),
      list("overview", "heatmap"),
      list("permz", "tt"),
      list("permz", "fc"),
      list("dimred", "tsne")
    )
    # check mode of interface (depends on timeseries /yes/no and bivariate/multivariate)
    # then show the relevent tabs
    # TODO: enable multivariate time series analysis
    if(is.null(interface$mode)) {
      show.tabs <- hide.tabs[1]
    }else if(interface$mode == 'multivar'){
      show.tabs <- hide.tabs[c(1,2,3,6,7,9,10,13)]
      heatbutton$status <- NULL
      #show.tabs <- c("inf","pca", "aov", "heatmap", "enrich", "venn")
    }else if(interface$mode == 'bivar'){
      show.tabs <- hide.tabs[c(1,2,3,7,8,9,10,11,12,13)]
      heatbutton$status <- "ttfc"
      #show.tabs <- c("inf","pca", "plsda", "tt", "fc", "volc", "heatmap", "ml", "enrich", "venn")
    }else if(interface$mode == 'time'){
      show.tabs <- hide.tabs[c(1,2,4,5,6,7,9,10,13)]
      heatbutton$status <- "asmb"
      #show.tabs <- c("inf", "pca", "aov", "asca", "meba", "heatmap", "ml", "venn")
    }else{
      show.tabs <- hide.tabs[1]
      #show.tabs <- c("inf") # 'info' tab that loads when no data is loaded currently
    }

    # hide all the tabs to begin with
    for(tab in hide.tabs){
      hideTab(inputId = unlist(tab)[1], 
              unlist(tab)[2], 
              session = session)
    }

    i=1
    # show the relevant tabs
    for(tab in show.tabs){
      showTab(inputId = unlist(tab)[1], unlist(tab)[2], session = session, select = ifelse(i==1, TRUE, FALSE))
      #showTab(inputId = "statistics", tab, select = ifelse(i==1, TRUE, FALSE), session = session)
      i = i + 1
    }

  })

  # -----------------
  observeEvent(input$undo_match_filt, {
    shown_matches$table <- global$tables$last_matches
  })

  # change ppm accuracy, ONLY USEFUL if loading in from CSV
  # TODO: make it possible to change this and re-make user database (mzranges table specifically)
  observeEvent(input$set_ppm, {
    ppm <<- input$ppm
    # show ppm amount in UI
    output$ppm <- renderText(ppm)
    # change in options file
    setOption(key="ppm", value=ppm)
  })

  # triggers when probnorm or compnorm is selected
  # let user pick a reference condition
  ref.selector <- reactive({
    # -------------
    if(input$norm_type == "ProbNorm" | input$norm_type == "CompNorm"){
      fluidRow(
        hr(),
        selectInput('ref_var',
                    'What is your reference condition?',
                    choices = c("")),
        actionButton("check_csv",
                     "Get options",
                     icon=icon("search")),
        hr()
      )
    }
  })

  # triggers when check_csv is clicked - get factors usable for normalization
  observeEvent(input$check_csv, {
    req(global$paths$csv_loc)
    switch(input$norm_type,
           ProbNorm=updateSelectInput(session, "ref_var",
                                      choices = get_ref_vars(fac = "label") # please add options for different times later, not difficult
           ),
           CompNorm=updateSelectInput(session, "ref_var",
                                      choices = get_ref_cpds() # please add options for different times later, not difficult
           ))
  })

  # render the created UI
  output$ref_select <- renderUI({ref.selector()})

  # init all observers
  for(fp in list.files("./backend/scripts/joanna/reactive", full.names = T)){
    #print(paste("Loading:", fp))
    source(fp, local = T)
  }

  # triggered when user enters the statistics tab

  observeEvent(input$dimred, {
    # check if an mset is present, otherwise abort
    if(!exists("mSet")) return(NULL)
    # depending on the present tab, perform analyses accordingly
    if(input$dimred %not in% names(mSet$analSet)){
      statsmanager$calculate <- input$dimred
    }
    datamanager$reload <- input$dimred
  })

  observeEvent(input$permz, {
    # check if an mset is present, otherwise abort
    if(!exists("mSet")) return(NULL)
    # depending on the present tab, perform analyses accordingly
    if(input$permz %not in% names(mSet$analSet)){
      statsmanager$calculate <- input$permz
    }
    datamanager$reload <- input$permz
  })

  observeEvent(input$overview, {
    # check if an mset is present, otherwise abort
    if(!exists("mSet")) return(NULL)
    # depending on the present tab, perform analyses accordingly
    if(input$overview %not in% names(mSet$analSet) | input$overview == "venn"){
      statsmanager$calculate <- input$overview
    }
    print(input$overview)
    datamanager$reload <- input$overview
  })

  observe({
    mdl = caret::getModelInfo()[[input$ml_method]]
    params <- mdl$parameters
    output$ml_params <- renderUI({
      list(
        helpText(mdl$label),
        hr(),
        h2("Tuning settings"),
        lapply(1:nrow(params), function(i){
          row = params[i,]
          list(
            textInput(inputId = paste0("ml_", row$parameter),
                      label = row$parameter,
                      value=if(input$ml_method=="glmnet"){
              switch(row$parameter,
                     alpha = 1,
                     lambda = "0:1:0.01")
            }),
            helpText(paste0(row$label, " (", row$class, ")."))
          )
        })
      )
    })
  })

  observe({
    if(exists("mSet")){
      if("analSet" %in% names(mSet)){
        if("ml" %in% names(mSet$analSet)){
          # update choice menu
          choices = c()
          for(method in c("rf", "ls")){
            if(method %in% names(mSet$analSet$ml)){
              model.names = names(mSet$analSet$ml[[method]])
              choices <- c(choices, paste0(method, " - ", paste0(model.names)))
            }
          }
          updateSelectInput(session, "show_which_ml", choices = choices, selected = paste0(mSet$analSet$ml$last$method, " - ", mSet$analSet$ml$last$name))
        }
      }
    }
  })

  observeEvent(input$show_which_ml,{
    split.name = strsplit(input$show_which_ml, split = " - ")[[1]]
    mSet$analSet$ml$last$method <<- split.name[[1]]
    mSet$analSet$ml$last$name <<- split.name[[2]]
    datamanager$reload <- "ml"

  },ignoreNULL = T, ignoreInit = T)

  observeEvent(input$select_db_all, {

    dbs <- global$vectors$db_list[-which(global$vectors$db_list == "custom")]

    currently.on <- sapply(dbs, function(db){
      input[[paste0("search_", db)]]
    })

    if(any(currently.on)){
      set.to = F
    }else{
      set.to = T
    }

    for(db in dbs){
      updateCheckboxInput(session, paste0("search_", db), value = set.to)
    }

  })

  # render the database download area
  output$db_build_ui <- renderUI({
    dbs_per_line = 4
    max_col_width = 12
    rows = ceiling(length(global$vectors$db_list) / dbs_per_line)
    database_layout = lapply(1:rows, function(i){
      min_i = (dbs_per_line * i) - (dbs_per_line - 1)
      max_i = (dbs_per_line * i)
      if(max_i > length(global$vectors$db_list)) max_i <- length(global$vectors$db_list)
      # create 3 fluidrows followed by a break
      list(
        # row 1: name
        fluidRow(lapply(global$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center", h2(global$constants$db.build.info[[db]]$title))
        })),
        # row 2: description
        fluidRow(lapply(global$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center", helpText(global$constants$db.build.info[[db]]$description))
        })),
        # row 3: image
        fluidRow(lapply(global$vectors$db_list[min_i:max_i], function(db){
          if(db != "custom"){
            column(width=3,align="center", imageOutput(global$constants$db.build.info[[db]]$image_id, inline=T))
          }else{
            column(width=3,align="center", shinyWidgets::circleButton("make_custom_db",
                                                                    size = "lg",
                                                                    icon = icon("plus",class = "fa-lg"),
                                                                    style = "stretch",
                                                                    color = "default")
                   )
          }
        })),
        # row 4: button
        fluidRow(lapply(global$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center",
            if(!(db %in% c("magicball", "custom"))){
              list(actionButton(paste0("check_", db), "Check", icon = icon("check")),
              actionButton(paste0("build_", db), "Build", icon = icon("wrench")),
              br(),
              imageOutput(paste0(db, "_check"),inline = T))
            }else{
              helpText("")
            }
          )
        })),
        br(),br()
      )
    })
    # return
    database_layout
  })

  db_button_prefixes = c("search", "add", "enrich")

  # generate all the fadebuttons for the database selection
  lapply(db_button_prefixes, function(prefix){
    output[[paste0("db_", prefix, "_select")]] <- renderUI({
      fluidRow(
        lapply(global$vectors$db_list[-which(global$vectors$db_list == "custom")], function(db){
          which_idx = grep(sapply(global$constants$images, function(x) x$name), pattern = db) # find the matching image (NAME MUST HAVE DB NAME IN IT COMPLETELY)
          sardine(fadeImageButton(inputId = paste0(prefix, "_", db), img.path = basename(global$constants$images[[which_idx]]$path))) # generate fitting html
        })
      )
    })
  })

  # check if these buttons are selected or not
  lapply(db_button_prefixes, function(prefix){
    observe({
      # ---------------------------------
      db_path_list <- lapply(global$vectors$db_list[-which(global$vectors$db_list == "custom")], # go through the dbs defined in db_lists
                             FUN = function(db){
                               button_id = input[[paste0(prefix, "_", db)]]
                               if(is.null(button_id)){
                                 NA
                               }else{
                                 if(!button_id){
                                   c(file.path(getOptions()$db_dir, paste0(db, ".base.db"))) # add path to list of dbpaths
                                 }
                                 else{NA}
                               }
                             }
      )
      # save the selected database paths to global
      global$vectors[[paste0("db_", prefix, "_list")]] <<- db_path_list[!is.na(db_path_list)]
    })
  })

  # render icon for search bar
  output$find_mol_icon <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('www/search.png'))
    # Return a list containing the filename and alt text
    list(src = filename,
         width=70,
         height=70)
  }, deleteFile = FALSE)


  observeEvent(input$save_mset, {
    # save mset
    withProgress({
      fn <- paste0(tools::file_path_sans_ext(global$paths$patdb), ".metshi")
      save(mSet, file = fn)
    })
  })
  
  observeEvent(input$load_mset, {
    # load mset
    withProgress({
      fn <- paste0(tools::file_path_sans_ext(global$paths$patdb), ".metshi")
      load(fn)
    },env = .GlobalEnv)
    print("loading...")
    datamanager$reload <- "general"
  })
  
  observeEvent(input$debug, {
    input <<- isolate(as.list(input))
    dput(isolate(as.list(input)))
  })
    
  observeEvent(input$ml_train_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    global$vectors$ml_train <<- c(input$subset_var, input$subset_group)
    print(global$vectors$ml_train)
    output$ml_train_ss <- renderText(subset.name)
  })

  observeEvent(input$ml_test_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    global$vectors$ml_test <<- c(input$subset_var, input$subset_group)
    print(global$vectors$ml_test)
    output$ml_test_ss <- renderText(subset.name)
  })

  # check if a dataset is already loaded in
  # change mode according to how many levels the experimental variable has
  # change interface based on that
  observe({
    if(exists("mSet")){
      if(is.null(mSet$timeseries)) mSet$timeseries <<- FALSE
      datamanager$reload <- "general"
    }else{
      # hide time series button
      timebutton$status <- "off"
      heatbutton$status <- "asmb"
    }
  })

  onStop(function() {
    print("closing metaboShiny ~ヾ(＾∇＾)")
    if(exists("mSet")){
      mSet$storage[[mSet$dataSet$cls.name]] <<- list()
      mSet$storage[[mSet$dataSet$cls.name]]$analysis <<- mSet$analSet
    }
    # remove metaboshiny csv files
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if(all(file.remove(rmv))) NULL
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
  observeEvent(input$make_custom_db, {

    # get window
    showModal(modalDialog(
      fluidRow(align="center",
               textInput("my_db_name", label = "Database full name", value = "MyDb"),
               textInput("my_db_short", label = "Database shorthand name", value = "mydb"),
               textInput("my_db_description", label = "Database description", value = "Custom database for MetaboShiny."),
               hr(),
               helpText("Please input a CSV file with at these columns (example below):"),
               helpText("'baseformula', 'charge', 'compoundname', 'identifier', 'description', 'structure' (in SMILES!)"),
               div(DT::dataTableOutput("db_example"), style="font-size:60%"),
               br(),
               shinyFiles::shinyFilesButton("custom_db", title = "Please pick a .csv file", multiple = F, label = "Select"),
               hr(),
               helpText("Please upload a database logo"),
               shinyFilesButton("custom_db_img_path",
                                'Select image',
                                'Please select a png file',
                                FALSE),br(),
               imageOutput("custom_db_img", inline=T),br(),
               hr(),
               shinyWidgets::circleButton("build_custom_db",icon = icon("arrow-right", "fa-lg"), size = "lg")
      ),
      title = "Custom database creation",
      easyClose = TRUE,
      size = "l",
      footer = NULL
    )
    )

  })

  observeEvent(input$build_custom_db, {

    csv_path <- parseFilePaths(global$paths$volumes, input$custom_db)$datapath

    # build base db
    db.build.custom(
      db.name = input$my_db_name,
      db.short = input$my_db_short,
      db.description = input$my_db_description,
      db.icon = global$paths$custom.db.path,
      csv = csv_path
    )

    # build extended db
    build.extended.db(tolower(input$my_db_short),
                      outfolder = file.path(getOptions()$db_dir, "custom"),
                      adduct.table = adducts,cl = 0)

  })

  removeModal()

})
