# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

function(input, output, session) {
  
  library(data.table)
  library(plotly)
  library(ggplot2)
  library(SPARQL)
  library(MetaboShiny)
  library(MetaDBparse)
  library(plyr)
  library(dplyr)
  
  options(shiny.maxRequestSize=50000*1024^2)
  setTimeLimit(cpu = Inf)
  
  # EMPTY DEBUG
  debug_browse_content <- debug_result_filters <- debug_pieinfo <- debug_input <- debug_lcl <- debug_mSet <- debug_matches <- debug_enrich <- debug_selection <- debug_venn_yes <- debug_report_yes <- list()
  
  # FIX FOR NORMALIZATION
  OFFtoJSON <- function(obj, ...){
    print("Disabled in MetaboShiny!")
  }
  
  rlang::env_unlock(env = asNamespace('RJSONIO'))
  rlang::env_binding_unlock(env = asNamespace('RJSONIO'))
  assign('toJSON', OFFtoJSON, envir = asNamespace('RJSONIO'))
  rlang::env_binding_lock(env = asNamespace('RJSONIO'))
  rlang::env_lock(asNamespace('RJSONIO'))
  
  AddErrMsg <- function(msg){
    print(msg)
    try({
      shiny::showNotification(msg)
    })
  }
  
  assignInNamespace("AddErrMsg", AddErrMsg, ns="MetaboAnalystR", 
                    envir=as.environment("package:MetaboAnalystR"))
  
  shiny::showNotification("Starting server process...")
  
  # detach("package:MetaboShiny", unload=T)
  # used to be in startshiny.R
  options("download.file.method" = "libcurl")
  options(expressions = 5e5)
  online = MetaboShiny::internetWorks()
  
  # make MetaboShiny_storage dir in home first..
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
  # with autorun
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R
  # current instructions
  # Rshiny app to analyse untargeted metabolomics data! BASH INSTRUCTIONS: STEP 1: mydir=~"/MetaboShiny" #or another of your choice | STEP 2: mkdir $mydir | STEP 3: docker run -p 8080:8080 -v $mydir:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /start.sh
  
  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'
  
  runmode <- if(file.exists(".dockerenv")) 'docker' else 'local'
  
  mSet <- NULL
  opts <- list()
  showtext::showtext_auto(enable = T)
  
  lcl = list(
    proj_name = "",
    last_mset = "",
    hasChanged = FALSE,
    load_ui = FALSE,
    last_dir = c(),
    lists = list(),
    prev_mz = "",
    prev_struct = "",
    beep = F,
    tables = list(last_matches = data.table::data.table(query_mz = "none"),
                  prev_pie = data.table::data.table()),
    functions = list(),
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
    #require(MetaboShiny)
    MetaboShiny::start_metshi(inBrowser = T)
  }
  
  shiny::showModal(loadModal())
  
  # == REACTIVE VALUES ==
  
  interface <- shiny::reactiveValues()
  interface$mode <- NULL
  
  mSetter <- shiny::reactiveValues(do = NULL)
  
  filemanager <- shiny::reactiveValues(do = "nothing")
  
  ggplotly <- shiny::reactiveValues(fancy = T)
  
  shown_matches <- shiny::reactiveValues(forward_full = data.table::data.table(),
                                         forward_unique = data.table::data.table(),
                                         reverse = data.table::data.table())
  
  ml_queue <- shiny::reactiveValues(jobs = list())
  
  enrich <- shiny::reactiveValues(overview =data.table::data.table("not run" = "run enrichment in sidebar :)"),
                                  current = data.table::data.table("nothing selected" = "please select a pathway to see selected m/z!"))
  
  browse_content <- shiny::reactiveValues(table = data.table::data.table())
  
  result_filters <- shiny::reactiveValues(add = list(positive = c(),
                                                     negative = c()),
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
  
  dbmanager <- shiny::reactiveValues(build = c("none"))
  
  search_button = shiny::reactiveValues(on=TRUE)
  
  search = shiny::reactiveValues(go = F)
  
  shiny::observe({
    shiny::showNotification("Loading user settings.")
    
    if(exists("lcl")){
      # - - - 
      home = path.expand('~')
      
      # create dir for options
      basedir=gsubfn::fn$paste("$home/MetaboShiny")
      if(!dir.exists(basedir)) dir.create(basedir,recursive = T)
      
      lcl$paths$opt.loc <<- file.path(basedir, "options.txt")
      old.wd = getwd()
      
      # == LOAD OPTIONS ==
      # look for existing source folder that DOESN'T MATCH the files
      if(!file.exists(lcl$paths$opt.loc)){
        shiny::showNotification("Welcome! Creating new user options file...")
        contents = gsubfn::fn$paste('db_dir = $home/MetaboShiny/databases
work_dir = $home/MetaboShiny/saves/admin
proj_name = MY_METSHI
ppm = 2
packages_installed = Y
font1 = Dosis
font2 = Dosis
font3 = Open Sans
font4 = Open Sans
col1 = #1961AB
col2 = #FFFFFF
col3 = #FFFFFF
col4 = #000000
size1 = 40
size2 = 20
size3 = 15
size4 = 11
taskbar_image = metshi_logo.png
gtheme = classic
gcols = #1C1400&#FFE552&#D49C1A&#EBC173&#8A00ED&#00E0C2&#95C200&#FF6BE4&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF
gspec = RdBu
gfont = 15
mode = complete
cores = 1
apikey =  
dbfavs =  
omit_unknown = yes
beep = no')
        writeLines(contents, lcl$paths$opt.loc)
      }
      
      opts <- MetaboShiny::getOptions(lcl$paths$opt.loc)
      
      # ==================
      
      if(opts$work_dir == ""){
        userfolder = file.path(home, "/MetaboShiny/saves")
      }else{
        userfolder = opts$work_dir
      }
      
      if(opts$db_dir == ""){
        dbdir = file.path(home, "/MetaboShiny/saves")
      }else{
        dbdir = opts$db_dir
      }
      
      if(!is.null(opts$add_paths)){
        print("Trying to add custom user paths to file selection menu")
        try({
          for(path in stringr::str_split(opts$add_paths, pattern = ";")[[1]]){
            print(path)
            
            spl_path = stringr::str_split(path, pattern = "\\:")[[1]]
            gbl$paths$volumes[spl_path[1]] <<-normalizePath(spl_path[2], mustWork=T)
          }  
        })
      }
      
      lcl$paths$work_dir <<- userfolder
      lcl$paths$db_dir <<- dbdir
      
      # - - - 
      if(!grepl("^//root", lcl$paths$work_dir)){
        print("Can't find 'root' directory -  activating docker > R autofix...")
        for(path in names(lcl$paths)){
          lcl$paths[[path]] <<- gsub("/root", "~", lcl$paths[[path]])
        }
      }else{
        if(!dir.exists(lcl$paths$work_dir)) dir.create(lcl$paths$work_dir,recursive = T)
        if(!dir.exists(lcl$paths$db_dir)) dir.create(lcl$paths$db_dir,recursive = T)
      }
      
      if("adducts.csv" %in% basename(list.files(lcl$paths$work_dir))){
        print("!")
        adducts <<- data.table::fread(file.path(lcl$paths$work_dir, "adducts.csv"))
      }
      
      if("adduct_rules.csv" %in% basename(list.files(lcl$paths$work_dir))){
        adduct_rules <<- data.table::fread(file.path(lcl$paths$work_dir, "adduct_rules.csv"))
      }
      
      adducts[adducts == ''|adducts == ' '] <<- NA
      
      lapply(c("mummi_adducts", 
               "score_adducts",
               "filter_adducts"), function(id){
                 cats = setdiff(colnames(adducts), colnames(MetaDBparse::adducts))
                 pickerID = paste0(id, "_cat_picker")
                 output[[paste0(id, "_cats")]] <- shiny::renderUI(if(length(cats) > 0){
                   shinyWidgets::checkboxGroupButtons(
                     inputId = pickerID,
                     label = "Adduct categories:",
                     choices = cats,
                     justified = TRUE,
                     checkIcon = list(
                       yes = icon("plus", 
                                  lib = "glyphicon"))
                   )
                 }else{
                   list()
                 })
                 # observers
                 shiny::observeEvent(input[[pickerID]],{
                   if(length(input[[pickerID]]) > 0){
                     sel_adducts = lapply(input[[pickerID]], function(colu){
                       col = adducts[[colu]] == "v"
                       col[is.na(col)] <- FALSE
                       col
                     })
                     if(length(sel_adducts) == 1){
                       sel_adducts = sel_adducts[[1]]
                     }else{
                       sel_adducts = lapply(sel_adducts, as.list)
                       sel_adduct_table = data.table::rbindlist(sel_adducts)
                       sel_adduct_sums = colSums(sel_adduct_table)
                       sel_adducts = sel_adduct_sums > 0
                     }
                     if(!grepl("filter", pickerID)){
                       shinyWidgets::updatePickerInput(session, id, selected = adducts$Name[sel_adducts])
                     }else{
                       # also do pie chart filter?
                       if(!is.null(pieinfo)){
                         adds_in_search = as.character(pieinfo$add$Var.1)
                         for(mzMode in c("positive", "negative")){
                           sel_adducts_mode = sel_adducts & adducts$Ion_mode == mzMode
                           keep_sel_adducts <- intersect(adducts$Name[sel_adducts_mode], adds_in_search)
                           result_filters$add[[mzMode]] <- keep_sel_adducts
                           search$go <- T
                         }  
                       }
                     }  
                   }else{
                     if(grepl("filter", pickerID)){
                       if(!is.null(pieinfo) & my_selection$mz != ""){
                         result_filters$add$positive <- result_filters$add$negative <- character(0)
                         search$go <- T 
                       }  
                     }else{
                       shinyWidgets::updatePickerInput(session, id, selected = character(0))
                     }
                   }
                 },ignoreNULL = FALSE)
               })
      
      addResourcePath('www', system.file('www', package = 'MetaboShiny'))
      tmp = tempdir()
      addResourcePath('tmp', tmp)
      
      shiny::observeEvent(input$nav_general, {
        if(input$nav_general == "help"){
          if(online){
            shiny::showNotification("Loading help file...")
            url <- "https://github.com/joannawolthuis/MetaboShiny/blob/master/README.md"
            source <- readLines(url, encoding = "UTF-8")
            parsed_doc <- XML::htmlParse(source, encoding = "UTF-8")
            helpFile = tempfile()
            XML::xpathSApply(parsed_doc, path = '//*[@id="readme"]/article', function(x) XML::saveXML(x,helpFile))
            
            html = suppressWarnings(paste0(readLines(helpFile),collapse="<br>"))
            html = gsub("([\\w|\\/]+)www", "www", html, perl=T)
            html = gsub("(img src=.+?) ", "\\1 class=\"readme_image\" ", html, perl = T)
            html = gsub("[ |]\\(Figure.*?\\)", "", html)
            html = gsub("id=\"user-content-", "id=\"", html)
            html = gsub("<svg.*?\\/svg>", "", html)
            html = gsub("\\/joannawolthuis", "https://github.com/joannawolthuis", html)
            html = gsub("Issues page here,", "Issues page on <a href=\"https://github.com/joannawolthuis/MetaboShiny/issues\">GitHub<\\/a>", html)
            
            output$help <- shiny::renderUI(shiny::HTML(html))  
          }else{
            output$help <- shiny::renderUI(shiny::h3("You are offline!"))
          }
        }
      })
      
      # - - - - check for custom databases - - - - 
      
      last_db = length(gbl$vectors$db_list) - 2
      
      files_db_folder = list.files(lcl$paths$db_dir,pattern = "_source|\\.db$")
      all_dbs = unique(gsub(files_db_folder, pattern = "_source|\\.db$", replacement=""))
      which.custom = which(!((all_dbs %in% gbl$vectors$db_list)|grepl(all_dbs,pattern="^ext.*$")))
      
      gbl$vectors$db_list <<- append(gbl$vectors$db_list, 
                                     values = all_dbs[which.custom], 
                                     after = last_db)
      
      for(db in all_dbs[which.custom]){
        
        dbfolder = paste0(db, "_source")
        imgpath = normalizePath(file.path(lcl$paths$db_dir, dbfolder, "logo.png"))
        
        # add image to gbl images
        imgname = paste0(paste0(db, '_logo'), ".png")
        file.copy(imgpath, file.path(tmp, imgname))
        gbl$constants$images <<- append(gbl$constants$images, values = 
                                          list(list(name = paste0(db, '_logo'), 
                                                    path = file.path("tmp", imgname),
                                                    dimensions = c(150, 150))))
        
        # add db info
        load(file = file.path(lcl$paths$db_dir, dbfolder, "info.RData"))
        gbl$constants$db.build.info[[db]] <<- dbinfo
      }
      
      gbl$vectors$db_categories$custom <<- all_dbs[which.custom]
      gbl$vectors$db_categories$all <<- gbl$vectors$db_list
      
      # create image objects in UI
      lapply(gbl$constants$images, FUN=function(image){
        output[[image$name]] <- shiny::renderImage({
          image$path <- if(grepl("tmp", image$path)) gsub("tmp", tmp, image$path) else image$path
          filename <- normalizePath(image$path)
          # Return a list containing the filename and alt text
          list(src = filename)
        }, deleteFile = FALSE)
      })
      
      db_section$load <- TRUE
      
      if(!("beep" %in% names(opts))){
        opts$beep <- "none"
        MetaboShiny::setOption(lcl$paths$opt.loc, "beep", "none")
      }
      
      seeb = 1337
      if("seed" %in% names(opts)){
        shiny::showNotification("Setting user-defined seed...")
        seeb <- as.numeric(opts$seed)
      }
      
      set.seed(seeb)
      lcl$seed <<- seeb
      shiny::updateNumericInput(session, "seed", value = seeb)
      
      options('unzip.unzip' = getOption("unzip"),
              'download.file.extra' = switch(runmode, 
                                             docker="--insecure",
                                             local=""),  # bad but only way to have internet in docker...
              'download.file.method' = 'curl',
              width = 1200, height=800)
      
      maxcores = if(Sys.getenv("SLURM_CPUS_ON_NODE") != ""){
        as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
      }else parallel::detectCores()
      
      shiny::updateSliderInput(session, "ncores",
                               value = min(maxcores-1, as.numeric(opts$cores)), 
                               min = 1, 
                               max = maxcores)
      
      if("adducts" %in% names(opts)){
        fav_adducts <- opts$adducts
        fav_adducts <- stringr::str_split(fav_adducts, pattern = "&")[[1]]  
      }else{
        fav_adducts <- adducts$Name
      }
      
      for(id in c("mummi_adducts", 
                  "score_adducts",
                  "fav_adducts")){
        shinyWidgets::updatePickerInput(session, id, selected = intersect(fav_adducts,adducts$Name))
      }
      
      # === GOOGLE FONT SUPPORT FOR GGPLOT2 ===
      
      # Download a webfont
      if(online){
        try({
          lapply(c(opts[grepl(pattern = "^font", names(opts))]), function(font){
            showtextdb::font_install(showtextdb::google_fonts(font))
          })  
        })
      }
      
      shiny::updateCheckboxInput(session, "omit_unknown", switch(opts$omit_unknown, yes=TRUE, no=FALSE))
      
      #if(opts$apikey != " "){
      #  output$api_set <- shiny::renderText("key saved!")
      #}
      
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
                  ".append('", 
                  '<style type="text/css">', 
                  font.css, bar.css, foot.css,
                  "</style>",
                  "');")
      
      shinyjs::runjs(jq)
      
      shinyjs::removeClass(class = "hidden",
                           id = "metshi")
      
      # init stuff that depends on opts file
      lcl$proj_name <<- opts$proj_name
      lcl$paths$proj_dir <<- file.path(lcl$paths$work_dir, 
                                       lcl$proj_name)
      lcl$paths$patdb <<- file.path(lcl$paths$proj_dir,
                                    paste0(opts$proj_name, ".db"))
      lcl$paths$csv_loc <<- file.path(lcl$paths$proj_dir,
                                      paste0(opts$proj_name, ".csv"))
      
      old.dir = getwd()
      on.exit({
        setwd(old.dir)
      })
      
      if(dir.exists(lcl$paths$work_dir)){
        setwd(lcl$paths$work_dir)
      }else{
        warning(paste("Please go to settings and re-specify work dir. Current does not exist:", lcl$paths$work_dir))
      }
      
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
      
      lcl$vectors$project_names <<- setdiff(basename(list.dirs(lcl$paths$work_dir)),"admin")
      
      lapply(1:4, function(i){
        colourpicker::updateColourInput(session=session,
                                        inputId = paste0("bar.col.", i),
                                        value = opts[[paste0("col", i)]])
        shiny::updateTextInput(session=session,
                               inputId = paste0("font.", i),
                               value = opts[[paste0("font", i)]])
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
      
      # fix for new font size in plots setting not being in options file yet
      if(!("gfont" %in% names(opts))){
        opts$gfont="15"
      }
      
      lcl$aes$font <<- list(family = opts$font4,
                            ax.num.size = as.numeric(opts$size4),
                            ax.txt.size = as.numeric(opts$size3),
                            ann.size = as.numeric(opts$size4),
                            title.size = as.numeric(opts$size2),
                            plot.font.size = as.numeric(opts$gfont))
      
      # create color pickers based on amount of colors allowed in global
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
      shiny::updateSelectInput(session, "beep", 
                               selected = opts$beep)
      shiny::updateSelectInput(session, "ggplot_theme", 
                               selected = opts$gtheme)
      shiny::updateSelectInput(session, "color_ramp", 
                               selected = opts$gspec)
      
      
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
      list("overview", "volcano"),#8
      list("overview", "venn"),#9
      list("overview", "heatmap"),#10
      list("overview", "power"),#11
      list("permz", "tt"),#12
      list("permz", "fc"),#13
      list("dimred", "tsne"),#14
      list("permz", "corr"),#15
      list("overview", "enrich"),#16
      list("dimred", "umap"),#17
      list("dimred", "ica"),#18
      list("overview", "featsel"),#19,
      list("permz", "cliffd")#20
    )
    # check mode of interface (depends on timeseries /yes/no and bivariate/multivariate)
    # then show the relevent tabs
    # TODO: enable multivariate time series analysis
    if(is.null(interface$mode)){
      show.tabs <- hide.tabs[1]
    }else if(interface$mode == '1fb'){
      show.tabs <- hide.tabs[c(1,2,3,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
      shiny::updateSelectInput(session, "ml_method",
                               selected = "rf",
                               choices = as.list(gbl$constants$ml.models))
    }else if(interface$mode == '1fm'){
      show.tabs <- hide.tabs[c(1,2,3,6,7,9,10,11,14,15,16,17,18,19)]
      shiny::updateSelectInput(session, "ml_method",
                               selected = "rf",
                               choices = as.list(setdiff(gbl$constants$ml.models,
                                                         gbl$constants$ml.twoonly)))
    }else if(interface$mode == '2f'){
      show.tabs <- hide.tabs[c(1,2,4,6,9,10,11,14,16,17,18,19)]
    }else if(interface$mode == 't1f'){
      show.tabs = hide.tabs[c(1,2,4,5,6,9,10,11,14,16,17,18,19)]
    }else if(interface$mode == 't'){
      show.tabs = hide.tabs[c(1,2,5,6,7,9,10,11,14,15,16,17,18,19)]
    }else{
      show.tabs <- hide.tabs[1]
    }
    
    # hide all the tabs to begin with
    if(length(show.tabs) > 1){
      for(bigtab in c("dimred", "permz", "overview", "ml")){
        shiny::showTab("statistics", bigtab)  
      }
      for(bigtab in c("search", "plot_aes", "switch", "subset", "metadata")){
        shiny::showTab("anal_sidebar", bigtab)  
      }
      shiny::hideTab("anal_sidebar", "start")
      shiny::updateTabsetPanel(session, "anal_sidebar", "switch")
    }else{
      for(bigtab in c("dimred", "permz", "overview", "ml")){
        shiny::hideTab("statistics", bigtab)  
      }
      for(bigtab in c("search", "plot_aes", "switch", "subset", "metadata")){
        shiny::hideTab("anal_sidebar", bigtab)
      }
    }
    
    for(tab in hide.tabs){
      shiny::hideTab(inputId = "statistics",#unlist(tab)[1],
                     unlist(tab)[2],
                     session = session)
    }
    
    i = 1
    # show the relevant tabs
    for(tab in show.tabs){
      shiny::showTab(inputId = "statistics",#unlist(tab)[1], 
                     unlist(tab)[2], session = session, 
                     select = ifelse(i==1, TRUE, FALSE))
      i = i + 1
    }
    
  })
  
  shiny::observe({
    if(my_selection$mz != ""){
      scanmode$positive <- F
      scanmode$negative <- F
      mzMode = if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
      for(mode in mzMode){
        scanmode[[mode]] <- TRUE
      }
    }
  })
  
  shiny::observeEvent(input$ncores, {
    try({
      shiny::showNotification("Stopping threads...")
      parallel::stopCluster(session_cl)
    })
    net_cores = input$ncores# - 1
    if(net_cores > 0){
      shiny::showNotification("Starting new threads...")
      logfile <<- tempfile()
      print(paste("Log file:", logfile))
      #if(file.exists(logfile)) file.remove(logfile)
      session_cl <<- parallel::makeCluster(net_cores,
                                           outfile=logfile)#,setup_strategy = "sequential") # leave 1 core for general use and 1 core for shiny session
      # send specific functions/packages to other threads
      parallel::clusterEvalQ(session_cl, {
        library(data.table)
        library(iterators)
        library(MetaboShiny)
        library(MetaDBparse)
      })  
    }else{
      session_cl <<- 0
    }
    
    MetaboShiny::setOption(lcl$paths$opt.loc, "cores", input$ncores)
  })
  
  shiny::observeEvent(input$set_api, {
    MetaboShiny::setOption(lcl$paths$opt.loc, "apikey", input$apikey)
    lcl$apikey <<- input$apikey
    output$api_set <- shiny::renderText("key saved!")
  })
  
  shiny::observeEvent(input$set_seed, {
    shiny::showNotification("Setting RNG seed...")
    MetaboShiny::setOption(lcl$paths$opt.loc, "seed", input$seed)
    lcl$seed <<- input$seed
    set.seed(as.numeric(lcl$seed))
    output$curr_seed <- shiny::renderText(as.character(lcl$seed))
  })
  
  shiny::observeEvent(input$beep, {
    if(input$nav_general == "options"){
      MetaboShiny::setOption(lcl$paths$opt.loc, "beep", input$beep)
      beepr::beep(sound = input$beep)
      lcl$aes$which_beep <<- input$beep
      lcl$beep <<- if(input$beep != "none") T else F      
    }
  },ignoreInit = T)
  
  shiny::observeEvent(input$omit_unknown, {
    new.val = if(input$omit_unknown) "yes" else "no"
    MetaboShiny::setOption(lcl$paths$opt.loc, "omit_unknown", new.val)
  })
  
  shiny::observeEvent(input$save_mset, {
    # save mset
    if(!is.null(mSet)){
      shiny::withProgress({
        filemanager$do <- "save"
      })  
    }
  })
  
  shiny::observeEvent(input$load_mset, {
    filemanager$do <- "load"
  })
  
  shiny::observeEvent(input$stats_type,{
    if(!is.null(mSet)){
      uimanager$refresh <- "statspicker"
    }
  })
  
  output$combi_anal1_picker <- shiny::renderUI({
    if(!is.null(input$combi_anal1)){
      anal = mSet$analSet[[input$combi_anal1]]
      target.mat = grep("\\.mat", names(anal),value = T)
      choices = colnames(anal[[target.mat]])
      shiny::selectInput("combi_anal1_var",label = "A result column:", choices=choices,selected = 1)
    }else{list()}
  })
  
  output$combi_anal2_picker <- shiny::renderUI({
    if(!is.null(input$combi_anal2)){
      anal = mSet$analSet[[input$combi_anal2]]
      target.mat = grep("\\.mat", names(anal),value = T)
      choices = colnames(anal[[target.mat]])
      shiny::selectInput("combi_anal2_var",label = "B result column:", choices=choices,selected = 1)
    }else{list()}
  })
  
  shiny::observeEvent(input$debug_metshi, {
    assign("lcl", lcl, envir = .GlobalEnv)
    assign("mSet", mSet, envir = .GlobalEnv)
    #assign("clientData", shiny::isolate(shiny::reactiveValuesToList(session$clientData)), envir = .GlobalEnv)
    assign("input", shiny::isolate(shiny::reactiveValuesToList(input)), envir = .GlobalEnv)
    assign("enrich", shiny::isolate(shiny::reactiveValuesToList(enrich)), envir = .GlobalEnv)
    assign("shown_matches", shiny::isolate(shiny::reactiveValuesToList(shown_matches)), envir = .GlobalEnv)
    assign("my_selection", shiny::isolate(shiny::reactiveValuesToList(my_selection)), envir = .GlobalEnv)
    assign("browse_content",  shiny::isolate(shiny::reactiveValuesToList(browse_content)), envir = .GlobalEnv)
    assign("pieinfo",  shiny::isolate(shiny::reactiveValuesToList(pieinfo)), envir = .GlobalEnv)
    assign("result_filters",  shiny::isolate(shiny::reactiveValuesToList(result_filters)), envir = .GlobalEnv)
    assign("report_yes",  shiny::isolate(shiny::reactiveValuesToList(report_yes)), envir = .GlobalEnv)
    assign("venn_yes",  shiny::isolate(shiny::reactiveValuesToList(venn_yes)), envir = .GlobalEnv)
    assign("ml_queue",  shiny::isolate(shiny::reactiveValuesToList(ml_queue)), envir = .GlobalEnv)
  })
  
  shiny::observeEvent(input$export_plot,{
    success=F
    try({
      orca_server$export(plotly::last_plot(), file.path(lcl$paths$proj_dir,paste0(basename(tempfile()), 
                                                                                  input$export_format)))
      success=T
    })
    if(!success) MetaboShiny::metshiAlert("Orca isn't working, please check your installation. If on Mac, please try starting Rstudio from the command line with the command 'open -a Rstudio'", session=session)
  })
  
  # triggered when user enters the statistics tab
  shinyjs::runjs('$("#mainPanel").resizable({
                                              handles: "e",
                                              resize: function() {
                                                $("#sidePanel").outerWidth($("#panelContainer").innerWidth() - $("#mainPanel").outerWidth());
                                              }
                                            });')
  
  observeEvent(input$statistics, { 
    if(!is.null(mSet)){
      if(!is.null(input$statistics)){
        
        uimanager$refresh <- input$statistics
        
        if(input$statistics %in% c("venn", "enrich", "heatmap", "network", "ml")){
          print("vennrich")
          statsmanager$calculate <- "vennrich"
          tablemanager$make <- "vennrich"
          uimanager$refresh <- c(input$statistics, "vennrich")
        }
        
        checkMe = input$statistics
        if(checkMe == "meba") checkMe = "MB"
        if(checkMe == "aov" & "aov2" %in% names(mSet$analSet)) checkMe = "aov2"
        
        if(checkMe %in% names(mSet$analSet)){
          shinyjs::show(selector = paste0("div.panel[value=collapse_", input$statistics, "_tables]"))
          tablemanager$make <- c(tablemanager$make, input$statistics)
          shinyjs::show(selector = paste0("div.panel[value=collapse_", input$statistics, "_plots]"))
          plotmanager$make <- c(plotmanager$make, input$statistics)
          shinyBS::updateCollapse(session, paste0("collapse_",input$statistics),open = paste0("collapse_", 
                                                                                              input$statistics, 
                                                                                              c("_tables","_plots")))
        }else{
          shinyBS::updateCollapse(session, paste0("collapse_",input$statistics),open = paste0("collapse_", 
                                                                                              input$statistics, 
                                                                                              "_settings")) 
          shinyjs::hide(selector = paste0("div.panel[value=collapse_", input$statistics, "_plots]"))
          shinyjs::hide(selector = paste0("div.panel[value=collapse_", input$statistics, "_tables]"))
        }
      }  
    }
  })
  
  analyses <- c("ml", "wordcloud", "plsda", 
                "pca", "tsne", "tt", "aov","combi",
                "fc", "volcano", "heatmap", 
                "meba", "asca", "corr", 
                "enrich", "network", "power",
                "umap", "ica", "featsel",
                "cliffd")
  
  lapply(analyses, function(an){
    shiny::observeEvent(input[[paste0("do_", an)]], {
      try({
        statsmanager$calculate <- c(statsmanager$calculate, an)
        shinyjs::show(selector = paste0("div.panel[value=collapse_", an, "_tables]"))
        tablemanager$make <- c(tablemanager$make, an)
        shinyjs::show(selector = paste0("div.panel[value=collapse_", an, "_plots]"))
        if(an != "featsel"){
          plotmanager$make <- c(plotmanager$make, an)
        }
        uimanager$refresh <- c(uimanager$refresh, an)
      })
    })    
  })
  
  # ==== LOAD LOGIN UI ====
  
  # init all observer
  for(fp in list.files("reactive", full.names = T)){
    source(fp, local = T)
  }  
  
  observeEvent(input$quit_metshi, {
    if(!is.null(mSet) & isolate(save_info$has_changed)){
      shinyWidgets::confirmSweetAlert(
        session = session,
        inputId = "save_exit",
        text = tags$div(
          tags$b("Click upper right ", icon("times"), " button to cancel."),br(),
          br(),
          shiny::img(class = "imagetop", 
                     src = "www/metshi_heart_bezel.png", 
                     height = "70px"),
          br()
        ),
        btn_labels = c("No", "Yes"),
        title = "Save before exiting?",
        #showCloseButton = T,
        html = TRUE
      )
    }else{
      shiny::stopApp()
      shinyjs::js$closeWindow()
    }
  })
  
  observeEvent(input$save_exit,{
    if(isTRUE(input$save_exit)){
      filemanager$do <- "save"
    }
    shiny::stopApp()
    shinyjs::js$closeWindow()
  },ignoreNULL = T)
  
  shinyjs::runjs("delay=1000;
                  setTimeoutConst = setTimeout(function(){
                                        $('#loading-page').fadeIn(1000);
                                      }, delay);
                  setInterval(function(){
                    if ($('html').attr('class')=='shiny-busy') {
                      setTimeoutConst = setTimeout(function(){
                        $('#loading-page').fadeIn(1000);
                      }, delay);
                    } else {
                      clearTimeout(setTimeoutConst);
                      $('#loading-page').hide();
                    }
                  },1000)")
  # new version check for either github or docker
  try({
    if(online){
      pkgs = c("MetaboShiny", "MetaDBparse")
      need.new = lapply(pkgs, function(pkg){
        remote = remotes:::github_remote(paste0("joannawolthuis/", pkg),
                                         host = "api.github.com",
                                         repos = getOption("repos"),
                                         type = getOption("pkgType"))
        package_name <- remotes:::remote_package_name(remote)
        local_sha <- remotes:::local_sha(package_name)
        remote_sha <- remotes:::remote_sha(remote, local_sha)
        need = local_sha != remote_sha
        msg = NA
        if(need){
          sha = remote_sha
          url = gsubfn::fn$paste("https://api.github.com/repos/joannawolthuis/$pkg/git/commits/$sha")
          msg = ""
          try({
            response = httr::GET(
              url
            )
            resParsed = httr::content(response, as = "parsed")
            msg <- resParsed$message
          }, silent=T)
        }
        list(need=need, msg=msg, pkg=pkg)
      })
      
      needs.new = sapply(need.new, function(x) x[[1]])
      if(any(needs.new)){
        alertContent = lapply(need.new, function(pkgRes){
          if(pkgRes$need){
            textID=tolower(paste0(pkgRes$pkg, "_update"))
            output[[textID]] <<- shiny::renderText(pkgRes$msg)
            if(pkgRes$msg != ""){
              shiny::tagList(shiny::h2(pkgRes$pkg),
                             shiny::wellPanel(shiny::helpText(pkgRes$msg))
              )  
            }else{
              list()
            }
          }else{
            list()
          }
        })
        if(".dockerenv" %in% list.files("/", all.files=T)){
          headerMsg = paste0("New ", paste(pkgs[needs.new], collapse = " and "), " version available! Please pull the latest docker image.")
        }else{
          headerMsg = paste0("New ", paste(pkgs[needs.new], collapse = " and "), " version available! Please install the latest GitHub version", if(sum(needs.new) == 2) "s." else ".")
        }
        metshiAlert(shiny::tagList(h3(headerMsg),
                                   hr(),
                                   alertContent), 
                    session = session,
                    title = "Notification")
      }
      
    }  
  },silent = T)
  
  observeEvent(input$fancy, {
    beFancy = !input$fancy
    if(beFancy){
      shinyjs::runjs("$(document).ready(startCursor())")#paste0(readLines(readPath),collapse="\n"))
      shinyjs::addClass(id = "fancy_pic", class = "imagetop")
      shinyjs::addClass(id = "appHeader",class = "rainbow")#"color-text-flow text")
    }else{
      removeUI(selector = "[id^=dots]", multiple=T)
      removeUI(selector = ".smallsparkle", multiple=T)
      shinyjs::removeClass(id="appHeader",
                           class="rainbow")
      shinyjs::removeClass(id="fancy_pic",
                           class="imagetop")
    }
  })
  
  onStop(function() {
    print("Closing MetaboShiny...")
    if(!is.null(session_cl)){
      parallel::stopCluster(session_cl)
    }
    session_cl <<- NULL
    gc()
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if(all(file.remove(rmv))) NULL
  })
}