# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

shinyServer(function(input, output, session) {

  mSet <- NULL
  bgcol <- "black"
  font.css <- ""
  bar.css <- ""
  opts <- list()

  logged <- reactiveValues(status = "notlogged",
                           text = "please log in! ( •́ .̫  •̀ )")

  lcl = list(
    curr_mz = "nothing selected",
    patdb = "",
    csv_loc = "",
    proj_name ="",
    last_mset="",
    tables=list(last_matches=data.table::data.table(query_mz = "none")),
    aes = list(font = list(),
               mycols = c(),
               spectrum = "rb",
               theme = "min"),
    vectors = list(proj_names = c()),
    paths = list(opt.loc = "",
                 patdb = "",
                 work_dir="")
    )

  ####### !!!!!!!!!!! #########
  metshi_mode <<- "one_user" #" multi_user" # for server mode
  
  observe({
    if(exists("lcl")){
      if(metshi_mode == "one_user"){
        print("Single-user mode activated~")
        # - - - 
        userfolder = "~/MetaboShiny/saves/admin"
        dbdir = "~/MetaboShiny/databases"
        # - - - 
        if(!dir.exists(userfolder)) dir.create(userfolder,recursive = T)
        if(!dir.exists(dbdir)) dir.create(dbdir,recursive = T)
        lcl$paths$opt.loc <<- file.path(userfolder, "options.txt")
        lcl$paths$work_dir <<- userfolder
        lcl$paths$db_dir <<- dbdir

        if(!file.exists(lcl$paths$opt.loc)){
          print("welp re-making options...")
          contents = gsubfn::fn$paste('db_dir = $dbdir
work_dir = $userfolder
proj_name = MY_METSHI
ppm = 2
packages_installed = Y
font1 = Pacifico
font2 = Pacifico
font3 = Open Sans
font4 = Open Sans
col1 = #000000
col2 = #DBDBDB
col3 = #FFFFFF
col4 = #FFFFFF
size1 = 40
size2 = 20
size3 = 15
size4 = 11
taskbar_image = gemmy_rainbow.png
gtheme = classic
gcols = #FF0004&#38A9FF&#FFC914&#2E282A&#8A00ED&#00E0C2&#95C200&#FF6BE4
gspec = RdBu
mode = complete')
          writeLines(contents, lcl$paths$opt.loc)
        }
        opts = getOptions(lcl$paths$opt.loc)
        logged$status <<- switch(opts$mode, dbonly = "db_only", complete= "logged")
      }
    }
  })
  
  output$login_status <- renderText({
    logged$text
  })

  # init all observers
  for(fp in list.files("./backend/scripts/joanna/reactive", full.names = T)){
    source(fp, local = T)
  }

  # ================================== LOGIN =====================================

  # default

  userdb = normalizePath("./users.db")

  # if logged in, check if exists
  observeEvent(input$login,{
    if(input$username != "" & input$password != ""){
      # get user role
      role = get_user_role(input$username, input$password)
      if(is.null(role)){
        logged$text <<- "wrong username/password (｡•́︿•̀｡)"
      }else{

        runmode <- if(file.exists(".dockerenv")) 'docker' else 'local'

        work_dir <- switch(runmode,
                           docker = "/userfiles/saves",
                           local = normalizePath("~/MetaboShiny/saves"))

        dbdir <- switch(runmode,
                           docker = "/userfiles/databases",
                           local = normalizePath("~/MetaboShiny/databases"))

        # check if user folder exists, otherwise make it
        userfolder = file.path(work_dir, input$username)

        if(!dir.exists(userfolder)){
          logged$text <<- "creating new user..."
          dir.create(userfolder)
        }

        logged$text <<- "logging in..."

        username = input$username

        # check if opts file exists, otherwise make it with the proper files
        lcl$paths$opt.loc <<- file.path(userfolder, "options.txt")
        lcl$paths$work_dir <<- userfolder
        lcl$paths$db_dir <<- dbdir

        if(!file.exists(lcl$paths$opt.loc)){
          contents = gsubfn::fn$paste('db_dir = $dbdir
work_dir = $userfolder
proj_name = MY_METSHI
ppm = 2
packages_installed = Y
font1 = Pacifico
font2 = Pacifico
font3 = Open Sans
font4 = Open Sans
col1 = #000000
col2 = #DBDBDB
col3 = #FFFFFF
col4 = #FFFFFF
size1 = 40
size2 = 20
size3 = 15
size4 = 11
taskbar_image = gemmy_rainbow.png
gtheme = classic
gcols = #FF0004&#38A9FF&#FFC914&#2E282A&#8A00ED&#00E0C2&#95C200&#FF6BE4
gspec = RdBu')
          writeLines(contents, lcl$paths$opt.loc)
        }

        logged$status <- "logged"
      }
    }else{
      logged$text <- "try again (｡•́︿•̀｡)"
    }
  })

  
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
  lapply(gbl$constants$default.text, FUN=function(default){
    output[[default$name]] = renderText(default$text)
  })

  showtext::showtext_auto() ## Automatically use showtext to render text for future devices

  # default match table fill
  output$match_tab <-DT::renderDataTable({
    if(is.null(lcl$tables$last_matches)){
      DT::datatable(data.table("no m/z chosen" = "Please choose m/z value from results ٩(｡•́‿•̀｡)۶	"),
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15),
                                   pageLength = 5))
    }
  })

  # create image objects in UI
  lapply(gbl$constants$images, FUN=function(image){
    output[[image$name]] <- renderImage({
      filename <- normalizePath(image$path)
      # Return a list containing the filename and alt text
      list(src = filename,
           width = image$dimensions[1],
           height = image$dimensions[2])
    }, deleteFile = FALSE)
  })

  shown_matches <- reactiveValues(table = data.table())

  output$match_tab <- DT::renderDataTable({
    remove_cols = gbl$vectors$remove_match_cols
    remove_idx <- which(colnames(shown_matches$table) %in% remove_cols)
    # don't show some columns but keep them in the original table, so they can be used
    # for showing molecule descriptions, structure
    DT::datatable(shown_matches$table,
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 15),
                                 pageLength = 5,
                                columnDefs = list(list(visible=FALSE, targets=remove_idx))
                  )
                  
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
    if(mSet$metshiParams$prematched){
      shown_matches$table <- get_prematches(mz = lcl$curr_mz,
                                            patdb = lcl$paths$patdb)
    }else{
      shown_matches$table <- lcl$tables$last_matches  
    }
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
    req(lcl$paths$csv_loc)
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

  # triggered when user enters the statistics tab
  
  observeEvent(input$statistics, {
    if(!is.null(mSet)){
      # check if an mset is present, otherwise abort
      switch(input$statistics,
             dimred = {
               if(!is.null(input$dimred)){
                 if(input$dimred %not in% names(mSet$analSet)){
                   statsmanager$calculate <- input$dimred
                 }
                 datamanager$reload <- input$dimred
               }
             },permz = {
               if(!is.null(input$permz)){
                 if(input$permz %not in% names(mSet$analSet)){
                   statsmanager$calculate <- input$permz
                 }
                 datamanager$reload <- input$permz
               }
             },overview = {
               if(!is.null(input$overview)){
                 if(input$overview %not in% names(mSet$analSet) | input$overview %in% c("heatmap", "venn")){
                    print(input$overview)
                    statsmanager$calculate <- input$overview
                 }
                 datamanager$reload <- input$overview
               }
             }, ml = {
               if(!is.null(input$ml)){
                 datamanager$reload <- "ml"
               }
             })
    }
  })
  
  
  observeEvent(input$tab_iden_4, {
    if(!is.null(mSet)){
      switch(input$tab_iden_4,
             pie_db = {
               statsmanager$calculate <- "match_pie"
               datamanager$reload <- "match_pie"
             },
             pie_add = {
               statsmanager$calculate <- "match_pie"
               datamanager$reload <- "match_pie"
             },
             word_cloud = {
               statsmanager$calculate <- "match_wordcloud"
               datamanager$reload <- "match_wordcloud"
             }) 
    }
  })
  
  observeEvent(input$dimred, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$dimred %not in% names(mSet$analSet)){
        statsmanager$calculate <- input$dimred
      }
      datamanager$reload <- input$dimred
    }
  })

  observeEvent(input$permz, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$permz %not in% names(mSet$analSet)){
        statsmanager$calculate <- input$permz
      }
      datamanager$reload <- input$permz
    }
  })

  observeEvent(input$overview, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$overview %not in% names(mSet$analSet) | input$overview == "venn"){
        statsmanager$calculate <- input$overview
      }
      datamanager$reload <- input$overview
    }
  })

  observeEvent(input$ml_top_x, {
    if(!is.null(mSet)){
      datamanager$reload <- "ml"
    }
  },ignoreInit = T, ignoreNULL = T)
  
  observe({
    if(!is.null(input$ml_method)){
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
    }
  })

  observeEvent(input$show_which_ml,{
    if(!is.null(mSet)){
      split.name = strsplit(input$show_which_ml, split = " - ")[[1]]
      mSet$analSet$ml$last$method <<- split.name[[1]]
      mSet$analSet$ml$last$name <<- split.name[[2]]
      datamanager$reload <- "ml"
    }
  },ignoreNULL = T, ignoreInit = T)

  observeEvent(input$select_db_all, {

    dbs <- gbl$vectors$db_list[-which(gbl$vectors$db_list %in% c("custom", "magicball"))]

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
    rows = ceiling(length(gbl$vectors$db_list) / dbs_per_line)
    database_layout = lapply(1:rows, function(i){
      min_i = (dbs_per_line * i) - (dbs_per_line - 1)
      max_i = (dbs_per_line * i)
      if(max_i > length(gbl$vectors$db_list)) max_i <- length(gbl$vectors$db_list)
      # create 3 fluidrows followed by a break
      list(
        # row 1: name
        fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center", h2(gbl$constants$db.build.info[[db]]$title))
        })),
        # row 2: description
        fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center", helpText(gbl$constants$db.build.info[[db]]$description))
        })),
        # row 3: image
        fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
          if(db != "custom"){
            column(width=3,align="center", imageOutput(gbl$constants$db.build.info[[db]]$image_id, inline=T))
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
        fluidRow(lapply(gbl$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center",
                 if(!(db %in% c("magicball", "custom"))){
                   list(actionButton(paste0("check_", db), "Check", icon = icon("check")),
                        actionButton(paste0("build_", db), "Build", icon = icon("wrench")),
                        br(),br(),
                        imageOutput(paste0(db, "_check"),inline = T))
                 }else{
                   helpText("")
                 }
          )
        })),
        br())
    })
    # return
    database_layout
  })

  db_button_prefixes = c("search", "add", "enrich", "prematch")

  # generate all the fadebuttons for the database selection
  lapply(db_button_prefixes, function(prefix){
    output[[paste0("db_", prefix, "_select")]] <- renderUI({
      built.dbs <- c(gsub(x = list.files(lcl$paths$db_dir, pattern = "\\.db"), 
                        pattern = "\\.db", replacement = ""), "custom")
      fluidRow(
        lapply(gbl$vectors$db_list[-which(gbl$vectors$db_list == "custom" | !(gbl$vectors$db_list %in% built.dbs))], function(db){
          which_idx = grep(sapply(gbl$constants$images, function(x) x$name), pattern = db) # find the matching image (NAME MUST HAVE DB NAME IN IT COMPLETELY)
          sardine(fadeImageButton(inputId = paste0(prefix, "_", db), img.path = basename(gbl$constants$images[[which_idx]]$path))) # generate fitting html
        })
      )
    })
  })

  # check if these buttons are selected or not
  lapply(db_button_prefixes, function(prefix){
    observe({
      # ---------------------------------
      db_path_list <- lapply(gbl$vectors$db_list[-which(gbl$vectors$db_list == "custom")], # go through the dbs defined in db_lists
                             FUN = function(db){
                               button_id = input[[paste0(prefix, "_", db)]]
                               if(is.null(button_id)){
                                 NA
                               }else{
                                 if(!button_id){
                                   c(file.path(lcl$paths$db_dir, paste0(db, ".db"))) # add path to list of dbpaths
                                 }
                                 else{NA}
                               }
                             }
      )
      # save the selected database paths to global
      lcl$vectors[[paste0("db_", prefix, "_list")]] <<- db_path_list[!is.na(db_path_list)]
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
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$patdb), ".metshi")
      if(exists("mSet")){
        save(mSet, file = fn)
      }
    })
  })

  observeEvent(input$load_mset, {
    # load mset
    withProgress({
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$patdb), ".metshi")
      if(file.exists(fn)){
        load(fn)
        mSet <<- mSet
        opts <- getOptions(lcl$paths$opt.loc)
        lcl$proj_name <<- opts$proj_name
        lcl$paths$patdb <<- file.path(opts$work_dir, paste0(opts$proj_name, ".db"))
        lcl$paths$csv_loc <<- file.path(opts$work_dir, paste0(opts$proj_name, ".csv"))
      }
      datamanager$reload <- "general"
    })
    # reload current plot
    updateNavbarPage(session, "statistics", selected = "inf")
  })

  observeEvent(input$debug, {
    debug_input <<- isolate(reactiveValuesToList(input))
    debug_lcl <<- lcl
    debug_mSet <<- mSet
  })

  observeEvent(input$ml_train_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    gbl$vectors$ml_train <<- c(input$subset_var,
                                  input$subset_group)
    output$ml_train_ss <- renderText(subset.name)
  })

  observeEvent(input$ml_test_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    gbl$vectors$ml_test <<- c(input$subset_var, input$subset_group)
    output$ml_test_ss <- renderText(subset.name)
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
  
  observeEvent(input$export_plot,{
    export_plotly(p = plotly::last_plot(), 
                  file=paste0(basename(tempfile()), input$export_format),
                  port = 9091)
  })

  observeEvent(input$build_custom_db, {

    csv_path <- parseFilePaths(gbl$paths$volumes, input$custom_db)$datapath

    # build base db
    db.build.custom(
      db.name = input$my_db_name,
      db.short = input$my_db_short,
      db.description = input$my_db_description,
      db.icon = gbl$paths$custom.db.path,
      csv = csv_path
    )

    # build extended db
    build.extended.db(tolower(input$my_db_short),
                      outfolder = file.path(lcl$paths$db_dir, "custom"),
                      adduct.table = adducts,cl = 0)

  })
  
  removeModal()

  # ==== ON EXIT ====
  
  onStop(function() {
    print("closing metaboShiny ~ヾ(＾∇＾)")
    debug_input <<- isolate(reactiveValuesToList(input))
    debug_lcl <<- lcl
    debug_mSet <<- mSet
    # remove metaboshiny csv files
    switch(runmode,
           local = {
             stop_orca(orca_serv_id)
           },
           docker = {
             print("how to stop orca server here??")
           })
    parallel::stopCluster(session_cl)
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if(all(file.remove(rmv))) NULL
  })
})
