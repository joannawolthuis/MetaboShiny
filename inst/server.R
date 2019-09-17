# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

function(input, output, session) {
  
  shiny::showModal(MetaboShiny::loadModal())
  
  print("...server...")
  #detach("package:MetaboShiny", unload=T)
  # used to be in startshiny.R
  options("download.file.method" = "libcurl")
  # make metaboshiny_storage dir in home first..
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
  # with autorun
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R
  # current instructions
  #Rshiny app to analyse untargeted metabolomics data! BASH INSTRUCTIONS: STEP 1: mydir=~"/MetaboShiny" #or another of your choice | STEP 2: mkdir $mydir | STEP 3: docker run -p 8080:8080 -v $mydir:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /start.sh
  
  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'
  
  runmode <- if(file.exists(".dockerenv")) 'docker' else 'local'
  
  options('unzip.unzip' = getOption("unzip"),
          'download.file.extra' = switch(runmode, 
                                         docker="--insecure",
                                         local=""),  # bad but only way to have internet in docker...
          'download.file.method' = 'curl',
          width = 1200, height=800)
  
  # - - - - - - - -
  
  mSet <- NULL
  opts <- list()
  showtext::showtext_auto(enable = T)
  
  logged <- shiny::reactiveValues(status = "notlogged",
                                  text = "please log in! ( •́ .̫  •̀ )")
  
  lcl = list(
    proj_name ="",
    last_mset="",
    load_ui = F,
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
  
  # == REACTIVE VALUES ==
  
  interface <- shiny::reactiveValues()
  
  shown_matches <- shiny::reactiveValues(forward_full = data.table::data.table(),
                                         forward_unique=data.table::data.table(),
                                         reverse = data.table::data.table())
  
  browse_content <- shiny::reactiveValues(table = data.table::data.table())
  
  result_filters <- shiny::reactiveValues(add = c(), 
                                          db = c(), 
                                          iso = c())
  
  my_selection = shiny::reactiveValues(mz = "",
                                       struct = "",
                                       revstruct = "",
                                       name = "")
  
  scanmode = shiny::reactiveValues(positive = FALSE,
                                   negative = FALSE)
  
  pieinfo = shiny::reactiveValues(add = c(),
                                  db = c(),
                                  iso = c())
  
  search_button = shiny::reactiveValues(on=TRUE)
  
  search = shiny::reactiveValues(go = F)
  
  ####### !!!!!!!!!!! #########
  shiny::observe({
    print("...loading user settings...")
    if(exists("lcl")){
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
      
      # = = = = = = =
      
      print("...loading ui...")
      opts <- MetaboShiny::getOptions(lcl$paths$opt.loc)
      
      # === GOOGLE FONT SUPPORT FOR GGPLOT2 ===
      online = MetaboShiny::internetWorks()
      
      # Download a webfont
      if(online){
        lapply(c(opts[grepl(pattern = "font", names(opts))]), function(font){
          try({
            sysfonts::font_add_google(name = font,
                                      family = font,
                                      regular.wt = 400,
                                      bold.wt = 700)
          })
        })
      }
      
      # parse color opts
      lcl$aes$mycols <<- MetaboShiny::get.col.map(lcl$paths$opt.loc) # colours for discrete sets, like group A vs group B etc.
      lcl$aes$theme <<- opts$gtheme # gradient function for heatmaps, volcano plot etc.
      lcl$aes$spectrum <<- opts$gspec # gradient function for heatmaps, volcano plot etc.
      
      # = = update css.. = =
      
      # generate CSS for the interface based on user settings for colours, fonts etc.
      bar.css <- MetaboShiny::nav.bar.css(opts$col1, opts$col2, opts$col3, opts$col4)
      font.css <- MetaboShiny::app.font.css(opts$font1, opts$font2, opts$font3, opts$font4,
                                            opts$size1, opts$size2, opts$size3,
                                            opts$size4, online=online, font.col=opts$col3)
      foot.css <- MetaboShiny::footer.css(opts$col1, font.col = opts$col3)
      
      # set taskbar image as set in options
      taskbar_image <- opts$task_img
      
      # $("head").append('<style type="text/css"></style>');
      jq = paste0('$("head")',
                  ".append('", '<style type="text/css">', 
                  bar.css, font.css, foot.css,
                  "</style>');")
      
      shinyjs::runjs(jq)
      
      shinyjs::removeClass(class="hidden",id = "metshi")
      
      # - - load custom dbs - -
      
      # load in custom databases
      has.customs <- dir.exists(file.path(lcl$paths$db_dir, "custom"))
      
      if(has.customs){
        
        customs = list.files(path = file.path(lcl$paths$db_dir, "custom"),
                             pattern = "\\.RData")
        
        dbnames = unique(tools::file_path_sans_ext(customs))
        
        for(db in dbnames){
          # add name to global
          dblist <- gbl$vectors$db_list
          dblist <- dblist[-which(dblist == "custom")]
          if(!(db %in% dblist)){
            dblist <- c(dblist, db, "custom")
            gbl$vectors$db_list <- dblist
          }
          metadata.path <- file.path(lcl$paths$db_dir, "custom", paste0(db, ".RData"))
          load(metadata.path)
          
          # add description to global
          gbl$constants$db.build.info[[db]] <- meta.dbpage
          
          # add image to global
          maxi = length(gbl$constants$images)
          gbl$constants$images[[maxi + 1]] <- meta.img
        }
      }
      
      # init stuff that depends on opts file
      
      lcl$proj_name <<- opts$proj_name
      lcl$paths$patdb <<- file.path(opts$work_dir,
                                    paste0(opts$proj_name, ".db"))
      lcl$paths$csv_loc <<- file.path(opts$work_dir,
                                      paste0(opts$proj_name, ".csv"))
      lcl$texts <<- list(
        list(name='curr_exp_dir', text=lcl$paths$work_dir),
        list(name='curr_db_dir', text=lcl$paths$db_dir),
        list(name='ppm', text=opts$ppm),
        list(name='proj_name', text=opts$proj_name)
      )
      
      lcl$vectors$project_names <<- unique(gsub(list.files(opts$work_dir,pattern = "\\.csv"),
                                                pattern = "(_no_out\\.csv)|(\\.csv)",
                                                replacement=""))
      
      lapply(1:4, function(i){
        colourpicker::updateColourInput(session=session,
                                        inputId = paste0("bar.col.", i),
                                        value = opts[[paste0("col", i)]])
      })
      
      lapply(1:4, function(i){
        shiny::updateTextInput(session=session,
                               inputId = paste0("font.", i),
                               value = opts[[paste0("font", i)]])
      })
      
      lapply(1:4, function(i){
        shiny::updateSliderInput(session=session,
                                 inputId = paste0("size.", i),
                                 value = opts[[paste0("size", i)]])
      })
      
      shiny::updateSelectInput(session, "ggplot_theme", selected = opts$gtheme)
      
      shiny::updateSelectInput(session, "color_ramp", selected = opts$gspec)
      
      #TODO: set these
      #shiny::div(class="plus", img(class="imagetop", src=opts$taskbar_image, width="100px", height="100px")),
      #shiny::div(class="minus", img(class="imagebottom", src=opts$taskbar_image, width="100px", height="100px"))
      
      shiny::updateCheckboxInput(session, "db_only", 
                                 value=switch(opts$mode, dbonly=T, complete=F))
      shiny::updateSelectizeInput(session,
                                  "proj_name",
                                  choices = lcl$vectors$project_names,
                                  selected = opts$proj_name)
      shiny::updateTextInput(session,"proj_name_new",
                             value = opts$proj_name)
      
      lapply(lcl$texts, FUN=function(default){
        output[[default$name]] = shiny::renderText(default$text)
      })
      
      lcl$aes$font <- list(family = opts$font4,
                           ax.num.size = 11,
                           ax.txt.size = 15,
                           ann.size = 20,
                           title.size = 25)
      
      # create color pickers based on amount of colours allowed in global
      output$colorPickers <- shiny::renderUI({
        lapply(c(1:gbl$constants$max.cols), function(i) {
          colourpicker::colourInput(inputId = paste("col", i, sep="_"),
                                    label = paste("Choose colour", i),
                                    value = lcl$aes$mycols[i],
                                    allowTransparent = F)
        })
      })
      
      # create color1, color2 etc variables to use in plotting functions
      # and update when colours picked change
      shiny::observe({
        values <- unlist(lapply(c(1:gbl$constants$max.cols), function(i) {
          input[[paste("col", i, sep="_")]]
        }))
        if(!any(is.null(values))){
          if(lcl$paths$opt.loc != ""){
            MetaboShiny::set.col.map(optionfile = lcl$paths$opt.loc, values)
            lcl$aes$mycols <<- values
          }
        }
      })
      
      shiny::updateSelectInput(session, "ggplot_theme", selected = opts$gtheme)
      shiny::updateSelectInput(session, "color_ramp", selected = opts$gspec)
      
    }
    
    shiny::removeModal()
    
  })
  
  
  
  # ================================== LOGIN =====================================
  
  output$manual_search <- shiny::renderUI({
    if(search_button$on){
      shiny::tags$button(
        id = "search_mz",
        class = "btn btn-default action-button",
        img(src = "detective.png",
            height = "50px")
      )
    }else{
      shiny::fluidRow(align="center", shiny::icon("paw","fa-s fa-rotate-90"), shiny::br(),
                      shiny::tags$i("pre-matched"))
    }
  })
  
  # ================================= DEFAULTS ===================================
  
  # set progress bar style to 'old' (otherwise it's not movable with CSS)
  shiny::shinyOptions(progress.style="old")
  
  options(digits=22,
          spinner.size = 0.5,
          spinner.type = 6,
          spinner.color = "black",
          spinner.color.background = "white")
  
  
  # send specific functions/packages to other threads
  parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
    "mape",
    "flattenlist"))
  
  parallel::clusterEvalQ(session_cl, library(data.table))
  
  # create default text objects in UI
  lapply(gbl$constants$default.text, FUN=function(default){
    output[[default$name]] = shiny::renderText(default$text)
  })
  
  showtext::showtext_auto() ## Automatically use showtext to render text for future devices
  
  # create image objects in UI
  lapply(gbl$constants$images, FUN=function(image){
    output[[image$name]] <- shiny::renderImage({
      filename <- normalizePath(image$path)
      # Return a list containing the filename and alt text
      list(src = filename,
           width = image$dimensions[1],
           height = image$dimensions[2])
    }, deleteFile = FALSE)
  })
  
  shiny::observe({
    # - - filters - -
    if(search$go){
      shiny::withProgress({

        matches = data.table::as.data.table(MetaboShiny::get_prematches(who = my_selection$mz,
                                                                        what = "query_mz",
                                                                        patdb = lcl$paths$patdb,
                                                                        showadd = result_filters$add,
                                                                        showdb = result_filters$db,
                                                                        showiso = result_filters$iso))
        
        shiny::setProgress(0.2)
        
        uniques = data.table::as.data.table(unique(data.table::as.data.table(matches)[,-c("source", "description"),with=F]))
        if(nrow(uniques)>1){
          uniques = uniques[, .(structure = list(structure)), 
                                    by = setdiff(names(uniques), "structure")]
          
          structs <- lapply(1:nrow(uniques), function(i){
            row = uniques[i,]
            opts = row$structure[[1]]
            if(length(opts)>1){
              # check for square brackets
              if(any(grepl(pattern = "\\[", x = opts))){
                opts[!grepl(pattern="\\[",x = opts)]
              }else{
                # just grab the first
                opts[1]
              }
            }else{
              opts
            }
          })
          uniques$structure <- structs
          
          }
          shiny::setProgress(0.4)
        
        shiny::setProgress(0.6)
        
        shown_matches$forward_unique <- uniques
        
        shown_matches$forward_full <- matches[,c("name", "source", "description"),with=F]
        
        
        if(nrow(shown_matches$forward_full)>0){
          pieinfo$add <- reshape::melt(table(shown_matches$forward_unique$adduct))
          pieinfo$db <- reshape::melt(table(shown_matches$forward_full$source))
          pieinfo$iso <- reshape::melt(table(shown_matches$forward_unique$isocat))
        }

        my_selection$name <- ""
        my_selection$struct <- ""
        
        shiny::setProgress(0.8)
        
      })
      search$go <- FALSE #reset self
    }
  })
  
  shiny::observe({
    if(my_selection$mz != ""){
      if(mSet$metshiParams$prematched){
        search$go <- TRUE
        my_selection$name <- ""
      }
    }
  })
  
  shiny::observe({
    if(my_selection$name != ""){
      
      subsec = data.table::as.data.table(shown_matches$forward_full)[name == my_selection$name]
      
      # render descriptions seperately
      output$desc_ui <- shiny::renderUI({lapply(1:nrow(subsec), function(i){
        row = subsec[i,]
        # icon(s) in one row
        dbs = strsplit(row$source, split=",")[[1]]
        img = sapply(dbs, function(db){
          id = gbl$constants$db.build.info[[db]]$image_id
          address = unlist(sapply(gbl$constants$images, function(item) if(item$name == id) item$path else NULL))
        })
        
        # text cloud underneath https://codepen.io/rikschennink/pen/mjywQb
        desc = row$description
        lapply(1:length(img), function(i){
          # ICONS IN A ROW
          name = names(img)[i]
          address = img[i]
          output[[paste0("desc_icon_", name)]] <- shiny::renderImage({
            # When input$n is 3, filename is ./images/image3.jpeg
            list(src = address,
                 #width=30,
                 height=30)
          }, deleteFile = FALSE)
        })
        
        desc_id = paste0("curr_desc_", dbs, collapse="_")
        output[[desc_id]] <- shiny::renderText({trimws(desc)})
        
        # ui
        list(shiny::fluidRow(align="center", 
                             lapply(dbs, function(db){
                               MetaboShiny::sardine(shiny::imageOutput(paste0("desc_icon_",db),inline = T))}),
                             shiny::br(),
                             textOutput(desc_id)),
             shiny::hr()
        )
      })
      }) 
    }
  })
  
  shiny::observe({
    if(my_selection$struct != ""){
      width = shiny::reactiveValuesToList(session$clientData)$output_empty_width
      if(width > 300) width = 300
      
      output$curr_struct <- shiny::renderPlot({plot.mol(my_selection$struct,
                                                        style = "cow")},
                                              width=width, 
                                              height=width) # plot molecular structure WITH CHEMMINER
      output$curr_formula <- shiny::renderText({unlist(shown_matches$forward_unique[curr_row,'baseformula'])}) # render text of current formula
    }
  })
  
  shiny::observe({
    if(my_selection$revstruct != ""){
      if(!mSet$metshiParams$prematched){
        print("Please perform pre-matching first to enable this feature!")
        return(NULL)
      }else{
        shown_matches$reverse <- unique(MetaboShiny::get_prematches(who = my_selection$revstruct,
                                                                    what = "map.structure", #map.mz as alternative
                                                                    patdb = lcl$paths$patdb)[,c("query_mz", "adduct", "%iso", "dppm")])
      } 
    }
  })
  
  # pie charts
  lapply(c("add", "db", "iso"), function(which_pie){
    output[[paste0("match_pie_", which_pie)]] <- plotly::renderPlotly({
      print("pies reloading...")
      pievec = pieinfo[[which_pie]]
      print(pievec)
      if(length(pievec)>0){
        plotly::plot_ly(pievec, labels = ~Var.1, 
                        values = ~value, size=~value*10, type = 'pie',
                        textposition = 'inside',
                        textinfo = 'label+percent',
                        insidetextfont = list(color = '#FFFFFF'),
                        hoverinfo = 'text',
                        text = ~paste0(Var.1, ": ", value, ' matches'),
                        marker = list(colors = colors,
                                      line = list(color = '#FFFFFF', width = 1)),
                        #The 'pull' attribute can also be used to create space between the sectors
                        showlegend = FALSE) %>%
          layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))  
      }
    })
  })
  
  # ===== UI SWITCHER ====
  
  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  shiny::observe({
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
      shiny::hideTab(inputId = unlist(tab)[1],
                     unlist(tab)[2],
                     session = session)
    }
    
    i=1
    # show the relevant tabs
    for(tab in show.tabs){
      shiny::showTab(inputId = unlist(tab)[1], 
                     unlist(tab)[2], session = session, 
                     select = ifelse(i==1, TRUE, FALSE))
      #showTab(inputId = "statistics", tab, select = ifelse(i==1, TRUE, FALSE), session = session)
      i = i + 1
    }
    
  })
  
  # print current compound in sidebar
  output$curr_mz <- shiny::renderText(my_selection$mz)
  
  # make miniplot for sidebar with current compound
  output$curr_plot <- plotly::renderPlotly({
    if(my_selection$mz != ""){
      MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, 
                                 cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                                 styles = input$ggplot_sum_style,
                                 add_stats = input$ggplot_sum_stats, 
                                 col.fac = input$col_var,
                                 txt.fac = input$txt_var,
                                 plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                 font = lcl$aes$font)
    }
  })
  
  # -----------------
  shiny::observeEvent(input$undo_match_filt, {
    result_filters$add <- c()
    result_filters$db <- c()
    result_filters$iso <- c()
    search$go <- T
  })
  
  shiny::observe({
    if(my_selection$mz != ""){
      scanmode$positive <- F
      scanmode$negative <- F
      ion_mode <- getIonMode(my_selection$mz, lcl$paths$patdb)
      for(mode in ion_mode){
        scanmode[[mode]] <- TRUE
      }
    }
  })
  
  # triggers when probnorm or compnorm is selected
  # let user pick a reference condition
  ref.selector <- reactive({
    # -------------
    if(input$norm_type == "ProbNorm" | input$norm_type == "CompNorm"){
      shiny::fluidRow(
        shiny::hr(),
        selectInput('ref_var',
                    'What is your reference condition?',
                    choices = c("")),
        actionButton("check_csv",
                     "Get options",
                     icon=shiny::icon("search")),
        shiny::hr()
      )
    }
  })
  
  # triggers when check_csv is clicked - get factors usable for normalization
  shiny::observeEvent(input$check_csv, {
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
  
  shiny::observeEvent(input$statistics, {
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
  
  shiny::observeEvent(input$dimred, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$dimred %not in% names(mSet$analSet)){
        statsmanager$calculate <- input$dimred
      }
      datamanager$reload <- input$dimred
    }
  })
  
  shiny::observeEvent(input$permz, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$permz %not in% names(mSet$analSet)){
        statsmanager$calculate <- input$permz
      }
      datamanager$reload <- input$permz
    }
  })
  
  shiny::observeEvent(input$overview, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$overview %not in% names(mSet$analSet) | input$overview == "venn"){
        statsmanager$calculate <- input$overview
      }
      datamanager$reload <- input$overview
    }
  })
  
  shiny::observeEvent(input$ml_top_x, {
    if(!is.null(mSet)){
      datamanager$reload <- "ml"
    }
  },ignoreInit = T, ignoreNULL = T)
  
  shiny::observe({
    if(!is.null(input$ml_method)){
      mdl = caret::getModelInfo()[[input$ml_method]]
      params <- mdl$parameters
      output$ml_params <- renderUI({
        list(
          shiny::helpText(mdl$label),
          shiny::hr(),
          h2("Tuning settings"),
          lapply(1:nrow(params), function(i){
            row = params[i,]
            list(
              shiny::textInput(inputId = paste0("ml_", row$parameter),
                               label = row$parameter,
                               value=if(input$ml_method=="glmnet"){
                                 switch(row$parameter,
                                        alpha = 1,
                                        lambda = "0:1:0.01")
                               }),
              shiny::helpText(paste0(row$label, " (", row$class, ")."))
            )
          })
        )
      })
    }
  })
  
  shiny::observeEvent(input$show_which_ml,{
    if(!is.null(mSet)){
      split.name = strsplit(input$show_which_ml, split = " - ")[[1]]
      mSet$analSet$ml$last$method <<- split.name[[1]]
      mSet$analSet$ml$last$name <<- split.name[[2]]
      datamanager$reload <- "ml"
    }
  },ignoreNULL = T, ignoreInit = T)
  
  shiny::observeEvent(input$select_db_all, {
    
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
  
  shiny::observeEvent(input$select_db_prematch_all, {
    dbs <- gbl$vectors$db_list[-which(gbl$vectors$db_list %in% c("custom", "magicball"))]
    currently.on <- sapply(dbs, function(db){
      input[[paste0("prematch_", db)]]
    })
    if(any(currently.on)){
      set.to = F
    }else{
      set.to = T
    }
    for(db in dbs){
      updateCheckboxInput(session, paste0("prematch_", db), value = set.to)
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
            shiny::column(width=3,align="center", shiny::imageOutput(gbl$constants$db.build.info[[db]]$image_id, inline=T))
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
          shiny::column(width=3,align="center",
                        if(!(db %in% c("magicball", "custom"))){
                          list(actionButton(paste0("check_", db), "Check", icon = shiny::icon("check")),
                               actionButton(paste0("build_", db), "Build", icon = shiny::icon("wrench")),
                               shiny::br(),shiny::br(),
                               shiny::imageOutput(paste0(db, "_check"),inline = T))
                        }else{
                          shiny::helpText("")
                        }
          )
        })),
        shiny::br())
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
      shiny::fluidRow(
        lapply(gbl$vectors$db_list[-which(gbl$vectors$db_list == "custom" | !(gbl$vectors$db_list %in% built.dbs))], function(db){
          which_idx = grep(sapply(gbl$constants$images, function(x) x$name), pattern = db) # find the matching image (NAME MUST HAVE DB NAME IN IT COMPLETELY)
          MetaboShiny::sardine(fadeImageButton(inputId = paste0(prefix, "_", db), img.path = basename(gbl$constants$images[[which_idx]]$path))) # generate fitting html
        })
      )
    })
  })
  
  # check if these buttons are selected or not
  lapply(db_button_prefixes, function(prefix){
    shiny::observe({
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
  
  
  shiny::observeEvent(input$save_mset, {
    # save mset
    shiny::withProgress({
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$patdb), ".metshi")
      if(exists("mSet")){
        save(mSet, file = fn)
      }
    })
  })
  
  shiny::observeEvent(input$load_mset, {
    # load mset
    shiny::withProgress({
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$patdb), ".metshi")
      print(fn)
      if(file.exists(fn)){
        load(fn)
        mSet <<- mSet
        opts <- MetaboShiny::getOptions(lcl$paths$opt.loc)
        if("ml" %in% names(mSet$analSet)){
          datamanager$reload <- "ml"
        }
        lcl$proj_name <<- opts$proj_name
        lcl$paths$patdb <<- file.path(opts$work_dir, paste0(opts$proj_name, ".db"))
        lcl$paths$csv_loc <<- file.path(opts$work_dir, paste0(opts$proj_name, ".csv"))
      }
      datamanager$reload <- "general"
    })
    # reload current plot
    updateNavbarPage(session, "statistics", selected = "inf")
  })
  
  shiny::observeEvent(input$debug, {
    debug_input <<- shiny::isolate(shiny::reactiveValuesToList(input))
    debug_lcl <<- lcl
    debug_mSet <<- mSet
  })
  
  shiny::observeEvent(input$ml_train_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$vectors$ml_train <<- c(input$subset_var,
                               input$subset_group)
    output$ml_train_ss <- renderText(subset.name)
  })
  
  shiny::observeEvent(input$ml_test_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$vectors$ml_test <<- c(input$subset_var, input$subset_group)
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
                                                   'Please select a png file',
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
  
  shiny::observeEvent(input$export_plot,{
    switch(runmode,
           docker = plotly::orca(p = plotly::last_plot(), 
                                 file=file.path(lcl$paths$work_dir,paste0(lcl$proj_name, "_", basename(tempfile()), input$export_format))),
           local = export_plotly(p = plotly::last_plot(), 
                                 file=file.path(lcl$paths$work_dir,paste0(lcl$proj_name, "_", basename(tempfile()), input$export_format)),
                                 port = 9091))
  })
  
  # shiny::observeEvent(input$build_custom_db, {
  #   
  #   csv_path <- parseFilePaths(gbl$paths$volumes, input$custom_db)$datapath
  #   
  #   # build base db
  #   db.build.custom(
  #     db.name = input$my_db_name,
  #     db.short = input$my_db_short,
  #     db.description = input$my_db_description,
  #     db.icon = gbl$paths$custom.db.path,
  #     csv = csv_path
  #   )
  #   
  #   # build extended db
  #   build.extended.db(tolower(input$my_db_short),
  #                     outfolder = file.path(lcl$paths$db_dir, "custom"),
  #                     adduct.table = adducts,cl = 0)
  #   
  # })
  
  # ==== LOAD LOGIN UI ====
  
  # init all observer
  for(fp in list.files("reactive", full.names = T)){
    source(fp, local = T)
  }  
  
  # ==== ON EXIT ====
  
  onStop(function() {
    print("closing metaboShiny ~ヾ(＾∇＾)")
    debug_input <<- shiny::isolate(shiny::reactiveValuesToList(input))
    debug_lcl <<- lcl
    debug_mSet <<- mSet
    debug_matches <<- shown_matches
    debug_selection <<- my_selection
    
    # remove metaboshiny csv files
    switch(runmode,
           local = {
             stop_orca(orca_serv_id)
           },
           docker = {
             NULL
           })
    parallel::stopCluster(session_cl)
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if(all(file.remove(rmv))) NULL
  })
}
