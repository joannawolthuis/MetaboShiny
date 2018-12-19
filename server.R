# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

shinyServer(function(input, output, session) {
  
  # ================================= DEFAULTS ===================================
  
  source('./backend/scripts/joanna/shiny_general.R')
  
  # set progress bar style to 'old' (otherwise it's not movable with CSS)
  shinyOptions(progress.style="old")
  
  options(digits=22,
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
    DT::datatable(data.table("no m/z chosen" = "Please choose m/z value from results ٩(｡•́‿•̀｡)۶	"),
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 15), 
                                 pageLength = 5))
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
      set.col.map("user_options.txt", values)
      global$vectors$mycols <<- get.col.map("user_options.txt")
    }
  })
  
  # ===== UI SWITCHER ====
  
  # create interface mode storage object.
  interface <- reactiveValues()
  
  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  observe({
    
    # hide all tabs by default, easier to hide them and then make visible selectively
    hide.tabs <- c("inf", "pca", "plsda", 
                   "tt", "fc", "aov", 
                   "meba", "asca", "ml", 
                   "volc", "heatmap", "enrich", 
                   "venn")
    
    # check mode of interface (depends on timeseries /yes/no and bivariate/multivariate)
    # then show the relevent tabs
    # TODO: enable multivariate time series analysis
    if(is.null(interface$mode)) {
      show.tabs <- c("inf")
    }else if(interface$mode == 'multivar'){ 
      show.tabs <- c("pca", "aov", "heatmap", "enrich", "venn")
    }else if(interface$mode == 'bivar'){  
      show.tabs <- c("pca", "plsda", "tt", "fc", "volc", "heatmap", "ml", "enrich", "venn")
    }else if(interface$mode == 'time'){
      show.tabs <- c("pca", "aov", "asca", "meba", "heatmap", "ml", "venn")
    }else{
      show.tabs <- c("inf") # 'info' tab that loads when no data is loaded currently
    }
    
    # hide all the tabs to begin with
    for(tab in hide.tabs){
      hideTab(inputId = "statistics", tab, session = session)
    }
    i=1
    # show the relevant tabs
    for(tab in show.tabs){
      showTab(inputId = "statistics", tab, select = ifelse(i==1, TRUE, FALSE), session = session)
      i = i + 1
    }
  })
  
  # -----------------
  
  # change ppm accuracy, ONLY USEFUL if loading in from CSV
  # TODO: make it possible to change this and re-make user database (mzranges table specifically)
  observeEvent(input$set_ppm, {
    ppm <<- input$ppm
    # show ppm amount in UI
    output$ppm <- renderText(ppm)
    # change in options file
    setOption("user_options.txt", "ppm", ppm)
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
    print(paste("Loading:", fp))
    source(fp, local = T)
  }
  
  # triggered when user enters the statistics tab
  observeEvent(input$statistics, {
    
    # check if an mset is present, otherwise abort
    if(!exists("mSet")) return(NULL)
    
    # depending on the present tab, perform analyses accordingly
    if(input$statistics %not in% names(mSet$analSet)){
      statsmanager$calculate <- input$statistics
    }
    datamanager$reload <- input$statistics 
    
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
          column(width=3,align="center",imageOutput(global$constants$db.build.info[[db]]$image_id, inline=T))
        })),
        # row 4: button
        fluidRow(lapply(global$vectors$db_list[min_i:max_i], function(db){
          column(width=3,align="center", list(
            actionButton(paste0("check_", db), "Check", icon = icon("check")),
            actionButton(paste0("build_", db), "Build", icon = icon("wrench")), 
            br(),
            imageOutput(paste0(db, "_check"),inline = T)
          ))
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
        lapply(global$vectors$db_list, function(db){
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
      db_path_list <- lapply(global$vectors$db_list, # go through the dbs defined in db_lists
                             FUN = function(db){
                               button_id = input[[paste0(prefix, "_", db)]]
                               if(is.null(button_id)){
                                 NA
                               }else{
                                 if(!button_id){
                                   c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db"))) # add path to list of dbpaths
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
  
  # this SHOULD trigger on closing the app
  observe({
    if (input$nav_general == "stop"){ # if on the 'close' tab...
      # interrupt and ask if you're sure 
      # TODO: save prompt here too
      shinyalert::shinyalert(title = "Question", 
                             text = "Do you want to close metaboShiny?", 
                             type = "warning",
                             #imageUrl = "www/question.png", 
                             showCancelButton = T, 
                             cancelButtonText = "No",
                             showConfirmButton = T, 
                             confirmButtonText = "Yes",
                             callbackR = function(x) {if(x == TRUE){
                               print("closing metaboShiny ~ヾ(＾∇＾)")
                               if(any(!is.na(session_cl))) parallel::stopCluster(session_cl) # close parallel threads
                               R.utils::gcDLLs() # flush dlls 
                               stopApp() 
                             }else{
                               NULL # if not, do nothing
                             }})
    }
  })
  
  removeModal()
  
})
