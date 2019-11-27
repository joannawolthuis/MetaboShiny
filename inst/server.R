# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

function(input, output, session) {
  
  shiny::showModal(MetaboShiny::loadModal())
  
  shiny::showNotification("Starting server process...")
  
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
                                  text = "please log in!")
  
  lcl = list(
    proj_name ="",
    last_mset="",
    load_ui = F,
    tables=list(last_matches=data.table::data.table(query_mz = "none")),
    aes = list(font = list(),
               mycols = c(),
               spectrum = "rb",
               theme = "min"),
    vectors = list(proj_names = c(),
                   built_dbs = gbl$vectors$db_list,
                   custom_db_list = c()),
    paths = list(opt.loc = "",
                 patdb = "",
                 work_dir="")
  )
  
  reinstall_restart <- function(){
    devtools::install(reload = T,upgrade = F)
    require(MetaboShiny)
    start.metshi(inBrowser = T)
  }
  
  # == REACTIVE VALUES ==
  
  interface <- shiny::reactiveValues()

  mSetter <- shiny::reactiveValues(do = NULL)
  
  shown_matches <- shiny::reactiveValues(forward_full = data.table::data.table(),
                                         forward_unique = data.table::data.table(),
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
  
  db_section = shiny::reactiveValues(load = FALSE)
  
  search_button = shiny::reactiveValues(on=TRUE)
  
  search = shiny::reactiveValues(go = F)
  
  ####### !!!!!!!!!!! #########
  shiny::observe({
    shiny::showNotification("Loading user settings.")
    
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
      
      if("adducts.csv" %in% basename(list.files(lcl$paths$work_dir))){
        adducts <<- fread(file.path(lcl$paths$work_dir, "adducts.csv"))
        }
      
      if("adduct_rules.csv" %in% basename(list.files(lcl$paths$work_dir))){
        adduct_rules <<- fread(file.path(lcl$paths$work_dir, "adduct_rules.csv"))
        }
      
      adducts[adducts == ''|adducts == ' '] <<- NA
      
      addResourcePath('www', system.file('www', package = 'MetaboShiny'))
      
      # - - - - check for custom databases - - - - 
      
      last_db = length(gbl$vectors$db_list) - 2
      
      files_db_folder = list.files(lcl$paths$db_dir)
      all_dbs = unique(gsub(files_db_folder, pattern = "_source|\\.db", replacement=""))
      which.custom = which(!((all_dbs %in% gbl$vectors$db_list)|grepl(all_dbs,pattern="^ext.*$")))
      
      gbl$vectors$db_list <<- append(gbl$vectors$db_list, 
                                     values = all_dbs[which.custom], 
                                     after = last_db)
      
      for(db in all_dbs[which.custom]){
        dbfolder = paste0(db, "_source")
        shiny::addResourcePath(prefix = db,
                               directoryPath = normalizePath(file.path(lcl$paths$db_dir, dbfolder)))
        
        # add image to gbl images
        gbl$constants$images <<- append(gbl$constants$images, values = 
                                          list(list(name = paste0(db, '_logo'), 
                                                    path = file.path(db, "logo.png"),
                                                    dimensions = c(150, 150))))
        
        # add db info
        load(file = file.path(lcl$paths$db_dir, dbfolder, "info.RData"))
        gbl$constants$db.build.info[[db]] <<- dbinfo
      }
      
      db_section$load <- TRUE
      
      # look for existing source folder that DOESN'T MATCH the files
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
gcols = #1C1400&#FFE552&#D49C1A&#EBC173&#8A00ED&#00E0C2&#95C200&#FF6BE4&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF
gspec = RdBu
mode = complete
cores = 1
apikey = ')
        writeLines(contents, lcl$paths$opt.loc)
      }
      
      # = = = = = = =
      
      shiny::showNotification("Loading interface...")
      opts <- MetaboShiny::getOptions(lcl$paths$opt.loc)
      
      shiny::updateSliderInput(session, "ncores", value = as.numeric(opts$cores))
      # send specific functions/packages to other threads

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
      
      if(opts$apikey != ""){
        output$api_set <- shiny::renderText("key saved!")
      }
      
      # parse color opts
      lcl$aes$mycols <<- MetaboShiny::get.col.map(lcl$paths$opt.loc) # colours for discrete sets, like group A vs group B etc.
      lcl$aes$theme <<- opts$gtheme # gradient function for heatmaps, volcano plot etc.
      lcl$aes$spectrum <<- opts$gspec # gradient function for heatmaps, volcano plot etc.
      lcl$apikey <<- opts$apikey
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
      
      shinyjs::removeClass(class="hidden",
                           id = "metshi")
      
      # load in custom databases
      has.customs <- dir.exists(file.path(lcl$paths$db_dir, 
                                          "custom"))
      
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
          metadata.path <- file.path(lcl$paths$db_dir, 
                                     "custom", 
                                     paste0(db, ".RData"))
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
        list(name='curr_exp_dir', 
             text=lcl$paths$work_dir),
        list(name='curr_db_dir', 
             text=lcl$paths$db_dir),
        list(name='ppm', 
             text=opts$ppm),
        list(name='proj_name', 
             text=opts$proj_name)
      )
      
      lcl$vectors$project_names <<- unique(gsub(list.files(opts$work_dir,
                                                           pattern = "\\.csv"),
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
      
      shiny::updateSelectInput(session, 
                               "ggplot_theme", 
                               selected = opts$gtheme)
      
      shiny::updateSelectInput(session, 
                               "color_ramp", 
                               selected = opts$gspec)
      
      #TODO: set these
      #shiny::div(class="plus", img(class="imagetop", src=opts$taskbar_image, width="100px", height="100px")),
      #shiny::div(class="minus", img(class="imagebottom", src=opts$taskbar_image, width="100px", height="100px"))
      
      shiny::updateCheckboxInput(session, "db_only", 
                                 value=switch(opts$mode, 
                                              dbonly=T, 
                                              complete=F))
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
            MetaboShiny::set.col.map(optionfile = lcl$paths$opt.loc, 
                                     values)
            lcl$aes$mycols <<- values
          }
        }
      })
      
      shiny::updateSelectInput(session, "ggplot_theme", 
                               selected = opts$gtheme)
      shiny::updateSelectInput(session, "color_ramp", 
                               selected = opts$gspec)
      

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
      shiny::fluidRow(align="center", 
                      shiny::icon("paw","fa-s fa-rotate-90"), 
                      shiny::br(),
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
  
  # create default text objects in UI
  lapply(gbl$constants$default.text, 
         FUN=function(default){
    output[[default$name]] = shiny::renderText(default$text)
  })
  
  showtext::showtext_auto() ## Automatically use showtext to render text for future devices
  
  shiny::observe({
    # - - filters - -
    if(search$go){
      shiny::withProgress({

        matches = data.table::as.data.table(MetaboShiny::get_prematches(who = my_selection$mz,
                                                                        what = "map.query_mz",
                                                                        patdb = lcl$paths$patdb,
                                                                        showadd = result_filters$add,
                                                                        showdb = result_filters$db,
                                                                        showiso = result_filters$iso))
        shiny::setProgress(0.2)
        
        uniques = data.table::as.data.table(unique(data.table::as.data.table(matches)[,-c("source", "description"),
                                                                                      with=F]))
        if(nrow(uniques)>1){
          uniques = uniques[, .(structure = list(structure)),
                                    by = setdiff(names(uniques),
                                                 "structure")]
          
          structs <- lapply(1:nrow(uniques), function(i){
            row = uniques[i,]
            opts = row$structure[[1]]
            if(length(opts)>1){
              if(any(grepl(opts, pattern = "\\[.*]-?\\d",perl = T))){
                res = opts[!grepl(opts, pattern = "\\[.*]-?\\d", perl=T)][1]
              }else{
                res = opts[1]
              }
            }else{
              res = opts
            }
            res
          })
          uniques$structure <- unlist(structs)
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
        my_selection$form <- ""
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
      keep <- which(trimws(subsec$description) %not in% c("","Unknown","unknown", " ",
                                                  "Involved in pathways: . More specifically: . Also associated with compound classes:"))
      subsec <- subsec[keep,]
      if(nrow(subsec)>0){
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
      }else{
        output$desc_ui <- shiny::renderUI({
          helpText("No additional info available!")
        })
      }
     
    }
  })
  
  shiny::observe({
    if(!MetaDBparse::is.empty(my_selection$struct)){
      width = shiny::reactiveValuesToList(session$clientData)$output_empty_width
      if(width > 300) width = 300
      
      output$curr_struct <- shiny::renderPlot({plot.mol(my_selection$struct,
                                                        style = "cow")},
                                              width=width, 
                                              height=width) # plot molecular structure WITH CHEMMINER
      output$curr_formula <- shiny::renderUI({
        tags$div(
          HTML(gsub(my_selection$form,pattern="(\\d+)", replacement="<sub>\\1</sub>",perl = T))
          )
        }) # render text of current formula
    }
  })
  
  shiny::observe({
    if(my_selection$revstruct != ""){
      if(!mSet$metshiParams$prematched){
        MetaboShiny::metshiAlert("Please perform pre-matching first to enable this feature!")
        return(NULL)
      }else{
        rev_matches = MetaboShiny::get_prematches(who = my_selection$revstruct,
                                     what = "map.structure",
                                     patdb = lcl$paths$patdb)
        if(nrow(rev_matches)>0){
          shown_matches$reverse <- unique(rev_matches[,c("query_mz", "adduct", "%iso", "dppm")])
        }else{
          shown_matches$reverse <- data.table()
        }
      } 
    }
  })
  
  # pie charts
  lapply(c("add", "db", "iso"), function(which_pie){
    output[[paste0("match_pie_", which_pie)]] <- plotly::renderPlotly({
      pievec = pieinfo[[which_pie]]
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
  
  shiny::observeEvent(input$heattable,{
    if(!is.null(input$overview)){
      if(input$overview == "heatmap"){
        statsmanager$calculate <- "heatmap"
        datamanager$reload <- "heatmap"
      }  
    }
  })
  
  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  shiny::observe({
    # hide all tabs by default, easier to hide them and then make visible selectively
    hide.tabs <- list(
      list("statistics", "inf"),#1
      list("dimred", "pca"),#2
      list("dimred", "plsda"),#3
      list("permz", "asca"),#4
      list("permz", "meba"),#5
      list("permz", "aov"),#6
      list("statistics", "ml"),#7
      list("overview", "volc"),#8
      list("overview", "venn"),#9
      list("overview", "heatmap"),#10
      list("permz", "tt"),#11
      list("permz", "fc"),#12
      list("dimred", "tsne")#13
    )
    
    # check mode of interface (depends on timeseries /yes/no and bivariate/multivariate)
    # then show the relevent tabs
    # TODO: enable multivariate time series analysis
      if(is.null(interface$mode)) {
        show.tabs <- hide.tabs[1]
      }else if(interface$mode == '1fb'){
        show.tabs <- hide.tabs[c(1,2,3,7,8,9,10,11,12,13)]
        shiny::updateSelectInput(session, "ml_method",
                                 selected = "rf",
                                 choices = as.list(gbl$constants$ml.models))
      }else if(interface$mode == '1fm'){
        show.tabs <- hide.tabs[c(1,2,3,6,7,9,10,13)]
        shiny::updateSelectInput(session, "ml_method",
                                 selected = "rf",
                                 choices = as.list(setdiff(gbl$constants$ml.models, gbl$constants$ml.twoonly)))
      }else if(interface$mode == '2f'){
        show.tabs <- hide.tabs[c(1,2,4,6,9,10,13)]
      }else if(interface$mode == 't1f'){
        show.tabs = hide.tabs[c(1,2,4,5,6,9,10,13)]
      }else if(interface$mode == 't'){
        show.tabs = hide.tabs[c(1,2,5,6,7,9,10,13)]
      }else{
        show.tabs <- hide.tabs[1]
      }
      
      # hide all the tabs to begin with
      for(tab in hide.tabs){
        shiny::hideTab(inputId = unlist(tab)[1],
                       unlist(tab)[2],
                       session = session)
      }
      
      i = 1
      # show the relevant tabs
      for(tab in show.tabs){
        shiny::showTab(inputId = unlist(tab)[1], 
                       unlist(tab)[2], session = session, 
                       select = ifelse(i==1, TRUE, FALSE))
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
  
  output$curr_add <- shiny::renderText({
    paste0(result_filters$add, collapse=", ")
  })
  output$curr_iso <- shiny::renderText({
    paste0(result_filters$iso, collapse=", ")
  })
  output$curr_db <- shiny::renderText({
    paste0(result_filters$db, collapse=", ")
  })
  
  shiny::observe({
    if(my_selection$mz != ""){
      scanmode$positive <- F
      scanmode$negative <- F
      ion_mode <- MetaboShiny::getIonMode(my_selection$mz, lcl$paths$patdb)
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
           ProbNorm=shiny::updateSelectInput(session, "ref_var",
                                      choices = get_ref_vars(fac = "label") # please add options for different times later, not difficult
           ),
           CompNorm=shiny::updateSelectInput(session, "ref_var",
                                      choices = get_ref_cpds() # please add options for different times later, not difficult
           ))
  })
  
  # render the created UI
  output$ref_select <- shiny::renderUI({ref.selector()})
  
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
                 if(input$overview %not in% names(mSet$analSet) | input$overview ==  "venn"){
                   statsmanager$calculate <- input$overview
                 }
                 datamanager$reload <- input$overview
               }
             }, ml = {
               datamanager$reload <- "ml"
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
  
  shiny::observeEvent(input$ncores, {
    if(!is.null(session_cl)){
      
      shiny::showNotification("Stopping threads...")
      parallel::stopCluster(session_cl)
    }
    shiny::showNotification("Starting new threads...")
    
    session_cl <<- parallel::makeCluster(input$ncores)#,outfile="") # leave 1 core for general use and 1 core for shiny session
    # send specific functions/packages to other threads
    parallel::clusterEvalQ(session_cl, library(data.table))
    MetaboShiny::setOption(lcl$paths$opt.loc, "cores", input$ncores)
  })
  
  shiny::observeEvent(input$set_api, {
    MetaboShiny::setOption(lcl$paths$opt.loc, "apikey", input$apikey)
    lcl$apikey <<- input$apikey
    output$api_set <- shiny::renderText("key saved!")
  })
  
  shiny::observeEvent(input$ml, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      datamanager$reload <- input$overview
    }
  })
  
  shiny::observeEvent(input$tab_iden_4, {
    # check if an mset is present, otherwise abort
    if(!is.null(mSet)){
      # depending on the present tab, perform analyses accordingly
      if(input$tab_iden_4 == "word_cloud"){
        if(nrow(shown_matches$forward_full) > 0){
          statsmanager$calculate <- "match_wordcloud"
          datamanager$reload <- "match_wordcloud"
        }
      }
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
  
  shiny::observeEvent(input$save_mset, {
    # save mset
    if(!is.null(mSet)){
      shiny::withProgress({
        fn <- paste0(tools::file_path_sans_ext(lcl$paths$patdb), ".metshi")
        if(exists("mSet")){
          save(mSet, file = fn)
        }
      })  
    }
  })
  shiny::observeEvent(input$stats_type,{
    if(!is.null(mSet)){
      datamanager$reload <- "statspicker"
    }
    })
  
  shiny::observeEvent(input$load_mset, {
    # load mset
    shiny::withProgress({
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$patdb), ".metshi")
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
        shiny::updateCheckboxInput(session,
                                   "paired",
                                   value = mSet$dataSet$paired)
      }
      datamanager$reload <- c("general","statspicker")
    })
    # reload current plot
    shiny::updateNavbarPage(session, "statistics", selected = "inf")
  })
  
  shiny::observeEvent(input$debug, {
    debug_input <<- shiny::isolate(shiny::reactiveValuesToList(input))
    debug_lcl <<- lcl
    debug_mSet <<- mSet
    debug_matches <<- shown_matches
    debug_selection <<- my_selection
  })
  
  shiny::observeEvent(input$ml_train_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$vectors$ml_train <<- c(input$subset_var,
                               input$subset_group)
    output$ml_train_ss <- shiny::renderText(subset.name)
  })
  
  shiny::observeEvent(input$ml_test_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$vectors$ml_test <<- c(input$subset_var, input$subset_group)
    output$ml_test_ss <- shiny::renderText(subset.name)
  })
  
  shiny::observeEvent(input$subset_var, {
    lvls = levels(as.factor(mSet$dataSet$covars[[input$subset_var]]))
    shiny::updateSelectizeInput(session, "subset_group", choices = lvls)
  },ignoreInit = T)
  
  shiny::observeEvent(input$export_plot,{
    success=F
    try({
      orca_server$export(plotly::last_plot(), file.path(lcl$paths$work_dir,paste0(lcl$proj_name, "_",
                                                                                  basename(tempfile()), 
                                                                                  input$export_format)))
      success=T
    })
    if(!success) MetaboShiny::metshiAlert("Orca isn't working, please check your installation. If on Mac, please try starting Rstudio from the command line with the command 'open -a Rstudio'", session=session)
  })
  
  # ==== LOAD LOGIN UI ====
  
  # init all observer
  for(fp in list.files("reactive", full.names = T)){
    source(fp, local = T)
  }  
  
  # ==== ON EXIT ====
  
  onStop(function() {
    print("Closing metaboShiny ~ヾ(＾∇＾)")
    debug_input <<- shiny::isolate(shiny::reactiveValuesToList(input))
    debug_lcl <<- lcl
    debug_mSet <<- mSet
    debug_matches <<- shown_matches
    debug_selection <<- my_selection
    if(!is.null(session_cl)){
      parallel::stopCluster(session_cl)
    }
    session_cl <<- NULL
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if(all(file.remove(rmv))) NULL
  })
}
