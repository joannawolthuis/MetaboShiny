# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

function(input, output, session) {
  
  library(data.table)
  library(plotly)
  library(ggplot2)
  library(SPARQL)
  library(MetaboShiny)
  
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
  
  # ====
  shiny::showModal(MetaboShiny::loadModal())
  
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
  
  # == REACTIVE VALUES ==
  
  interface <- shiny::reactiveValues()
  interface$mode <- NULL
  
  mSetter <- shiny::reactiveValues(do = NULL)
  
  ggplotly <- shiny::reactiveValues(fancy = T)
  
  shown_matches <- shiny::reactiveValues(forward_full = data.table::data.table(),
                                         forward_unique = data.table::data.table(),
                                         reverse = data.table::data.table())
  
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
  
  search_button = shiny::reactiveValues(on=TRUE)
  
  search = shiny::reactiveValues(go = F)
  
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
      
      old.wd = getwd()
      setwd(lcl$paths$work_dir)
      on.exit({
        setwd(old.wd)
      })
      
      if("adducts.csv" %in% basename(list.files(lcl$paths$work_dir))){
        adducts <<- data.table::fread(file.path(lcl$paths$work_dir, "adducts.csv"))
      }
      
      if("adduct_rules.csv" %in% basename(list.files(lcl$paths$work_dir))){
        adduct_rules <<- data.table::fread(file.path(lcl$paths$work_dir, "adduct_rules.csv"))
      }
      
      adducts[adducts == ''|adducts == ' '] <<- NA
      
      addResourcePath('www', system.file('www', package = 'MetaboShiny'))
      tmp = tempdir()
      addResourcePath('tmp', tmp)
      
      # - - - - - - - - - - - - - - - - - - - - - -
      
      
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
      
      # look for existing source folder that DOESN'T MATCH the files
      if(!file.exists(lcl$paths$opt.loc)){
        shiny::showNotification("Welcome! Creating new user options file...")
        contents = gsubfn::fn$paste('db_dir = $dbdir
work_dir = $userfolder
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
mode = complete
cores = 1
apikey =  
dbfavs =  
omit_unknown = yes')
        writeLines(contents, lcl$paths$opt.loc)
      }
      
      # = = = = = = =
      
      shiny::showNotification("Loading interface...")
      opts <- MetaboShiny::getOptions(lcl$paths$opt.loc)
      
      options('unzip.unzip' = getOption("unzip"),
              'download.file.extra' = switch(runmode, 
                                             docker="--insecure",
                                             local=""),  # bad but only way to have internet in docker...
              'download.file.method' = 'curl',
              width = 1200, height=800)
      
      shiny::updateSliderInput(session, "ncores", value = as.numeric(opts$cores))
      # send specific functions/packages to other threads
      
      
      fav_adducts <- opts$adducts
      fav_adducts <- stringr::str_split(fav_adducts, pattern = "&")[[1]]
      for(id in c("mummi_adducts", 
                  "score_adducts",
                  "fav_adducts")){
        shinyWidgets::updatePickerInput(session, id, selected = intersect(fav_adducts,adducts$Name))
      }
      
      # === GOOGLE FONT SUPPORT FOR GGPLOT2 ===
      
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
      lcl$paths$proj_dir <<- file.path(lcl$paths$work_dir, lcl$proj_name)
      lcl$paths$patdb <<- file.path(lcl$paths$proj_dir, paste0(opts$proj_name, ".db"))
      lcl$paths$csv_loc <<- file.path(lcl$paths$proj_dir, paste0(opts$proj_name, ".csv"))
      
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
      
      lcl$aes$font <<- list(family = opts$font4,
                            ax.num.size = as.numeric(opts$size4),
                            ax.txt.size = as.numeric(opts$size3),
                            ann.size = as.numeric(opts$size4),
                            title.size = as.numeric(opts$size2))
      
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
  
  observe({
    if(my_selection$mz != ""){
      for(pie in c("add", "iso","db")){
        if(pie == "add"){
          mzMode =if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
          result_filters$add[[mzMode]] <- result_filters$add[[mzMode]][!is.na(result_filters$add[[mzMode]])]
        }else{
          result_filters[[pie]] <- result_filters[[pie]][!is.na(result_filters[[pie]])]
        }
      }   
    }
  })
  
  shiny::observe({
    # - - filters - -
    if(search$go){
      shiny::withProgress({
        if(input$tab_iden_2 == "mzmol"){
          if(lcl$prev_mz != my_selection$mz & !identical(lcl$vectors$prev_dbs, lcl$vectors$db_search_list)){
            matches = data.table::as.data.table(MetaboShiny::get_prematches(who = gsub(my_selection$mz, 
                                                                                       pattern="/.*$|RT.*$", 
                                                                                       replacement=""),
                                                                            what = "map.query_mz",
                                                                            patdb = lcl$paths$patdb,
                                                                            showadd = c(),
                                                                            showdb = c(),
                                                                            showiso = c()))
            
            if(nrow(matches) == 0){
              shown_matches$forward_unique <- data.table::data.table()
              shown_matches$forward_full <- data.table::data.table()
              return(NULL)
            }
            lcl$prev_mz <<- my_selection$mz
            lcl$vectors$prev_dbs <<- lcl$vectors$db_search_list
            pieinfo$db <- reshape::melt(table(matches$source))
            pieinfo$add <- reshape::melt(table(matches$adduct))
            pieinfo$iso <- reshape::melt(table(matches$isocat))
          }
          
          mzMode = if(grepl(my_selection$mz, pattern="\\-")) "negative" else "positive"
          
          matches = data.table::as.data.table(MetaboShiny::get_prematches(who = gsub(my_selection$mz, pattern="/.*$|RT.*$", replacement=""),
                                                                          what = "map.query_mz",
                                                                          patdb = lcl$paths$patdb,
                                                                          showadd = result_filters$add[[mzMode]],
                                                                          showdb = result_filters$db,
                                                                          showiso = result_filters$iso))  
          
          if(nrow(matches)>0){
            
            shiny::setProgress(0.2)
            
            matches$compoundname[matches$source != "magicball"] <- tolower(matches$compoundname[matches$source != "magicball"])
            
            # =====
            
            uniques = data.table::as.data.table(unique(data.table::as.data.table(matches)[, -c("source", 
                                                                                               "description",
                                                                                               "identifier"),
                                                                                          with=F]))
            shiny::setProgress(0.6)
            
            # === aggregate ===
            
            info_only = unique(matches[,c("compoundname", 
                                          "source", 
                                          "structure", 
                                          "description",
                                          "identifier"),with=F])
            info_only$description <- paste0("Database ID: ", info_only$identifier, ". ", info_only$description)
            info_only <- unique(info_only[,-"identifier"])
            
            info_no_na <- info_only[!is.na(info_only$structure)]
            info_na <- info_only[is.na(info_only$structure)]
            
            info_aggr <- aggregate(info_no_na, by = list(info_no_na$compoundname), FUN = function(x) paste0(x, collapse = "SEPERATOR"))
            info_aggr <- aggregate(info_aggr, by = list(info_aggr$structure), FUN = function(x) paste0(x, collapse = "SEPERATOR"))
            
            info_aggr <- rbind(info_aggr, info_na, fill=T)
            
            # fix structures
            split_structs <- strsplit(info_aggr$structure, split = "SEPERATOR")
            main_structs <- unlist(lapply(split_structs, function(x) x[[1]]))
            info_aggr$structure <- main_structs
            
            # move extra names to descriptions
            split_names <- strsplit(info_aggr$compoundname, split = "SEPERATOR")
            main_names <- unlist(lapply(split_names, function(x) if(length(x) > 1) x[[1]] else x[1]))
            
            synonyms <- unlist(lapply(split_names, function(x){
              if(length(x)>1){
                paste0("SYNONYMS: ", paste0(unique(x[2:length(x)]), collapse=", "),".SEPERATOR") 
              }else NA
            }))
            
            info_aggr$compoundname <- main_names
            has.syn <- which(!is.na(synonyms))
            
            info_aggr$description[has.syn] <- paste0(synonyms[has.syn], info_aggr$description[has.syn])
            info_aggr <- data.table::as.data.table(info_aggr)
            
            # =================
            
            uniques = uniques[structure %in% info_aggr$structure]
            
            is.no.na.uniq <- which(!is.na(uniques$structure))
            is.no.na.info <- which(!is.na(info_aggr$structure))
            
            uniq.to.aggr <- match(uniques[is.no.na.uniq]$structure, 
                                  info_aggr[is.no.na.info]$structure)
            
            uniques$compoundname[is.no.na.uniq] <- info_aggr$compoundname[is.no.na.info][uniq.to.aggr]
            uniques$structure[is.no.na.uniq] <- info_aggr$structure[is.no.na.info][uniq.to.aggr]
            
            uniques <- unique(uniques)
            
            shown_matches$forward_unique <- uniques[,-grepl(colnames(uniques), pattern = "Group\\.\\d"),with=F]
            shown_matches$forward_full <- info_aggr[,-grepl(colnames(info_aggr), pattern = "Group\\.\\d"),with=F]
            
          }else{
            shown_matches$forward_unique <- data.table::data.table()
            shown_matches$forward_full <- data.table::data.table()
          }
          
          my_selection$name <- ""
          my_selection$form <- ""
          my_selection$struct <- ""  
          
        }else{
          NULL
        }
        
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
      # update star logo
      shinyWidgets::updatePrettyToggle(session, 
                                       "star_mz",
                                       value =  mSet$report$mzStarred[my_selection$mz]$star)
    }
  })
  
  shiny::observe({
    if(my_selection$name != ""){
      if(nrow(shown_matches$forward_full) > 0 ){
        subsec = data.table::as.data.table(shown_matches$forward_full)[compoundname == my_selection$name]
        
        if(grepl("SYNONYMS:", x = subsec$description)){
          has_syn = T
          subsec$source <- paste0("synonymSEPERATOR", subsec$source)
          subsec$description <- gsub(subsec$description, pattern = "SYNONYMS: ", replacement="")
        }else{
          has_syn = F
        }
        
        subsec <- subsec[, .(compoundname, 
                             source = strsplit(source, split = "SEPERATOR")[[1]], 
                             structure = structure, 
                             description = strsplit(description, split = "SEPERATOR")[[1]]
        )
        ]
        
        subsec <- aggregate(subsec, by = list(subsec$source), FUN=function(x) paste0(unique(x), collapse="."))
        
        keep <- which(trimws(subsec$description) %not in% c("","Unknown","unknown", " ",
                                                            "Involved in pathways: . More specifically: . Also associated with compound classes:"))
        subsec <- subsec[keep,]
        
        if(has_syn){
          subsec <- subsec[order(as.numeric(grepl(subsec$source, pattern = "synonym")), decreasing = T),]
        }
        
        if(nrow(subsec) > 0){
          
          # render descriptions seperately
          output$desc_ui <- shiny::renderUI({
            
            lapply(1:nrow(subsec), function(i){
              
              row = subsec[i,]
              
              # icon(s) in one row
              db = row$source
              desc_id = paste0("curr_desc_", db)
              desc = row$description
              #output[[desc_id]] <- shiny::renderText({desc})
              
              if(db == "synonym"){
                ui = shiny::fluidRow(align="center", 
                                     tags$h3("Synonyms:"),
                                     helpText(desc),
                                     shiny::br()
                )
              }else{
                id = gbl$constants$db.build.info[[db]]$image_id
                address = unlist(sapply(gbl$constants$images, function(item) if(item$name == id) item$path else NULL))
                ui = shiny::fluidRow(align="center", 
                                     shiny::tags$button(
                                       id = paste0(db, "_copy_id"),
                                       class = "btn btn-default action-button shiny-bound-input",
                                       shiny::img(src = address,
                                                  height = "30px"),
                                       style = "vertical-align: middle;border-radius: 0px;border-width: 0px;background-color: #ff000000;"
                                     ),
                                     shiny::br(),
                                     helpText(desc),
                                     shiny::br()
                )
                
                shiny::observeEvent(input[[paste0(db, "_copy_id")]], {
                  shiny::showNotification("Saving database identifier to clipboard!")
                  dbID = stringr::str_match(desc, "Database ID: (.*?). ")[,2]
                  clipr::write_clip(dbID, allow_non_interactive = TRUE)
                  shiny::updateTextInput(session,
                                         "wordcloud_searchTerm",
                                         value = dbID)
                })
              }
              return(ui)
            })
          })  
          
        }else{
          output$desc_ui <- shiny::renderUI({
            helpText("No additional info available!")
          })
        }
      }
    }
  })
  
  shiny::observe({
    if(!MetaDBparse::is.empty(my_selection$struct)){
      width = shiny::reactiveValuesToList(session$clientData)$output_empty_width
      if(width > 300) width = 300
      output$curr_struct <- shiny::renderPlot({plot_mol(my_selection$struct,
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
  
  shiny::observeEvent(input$curr_mz, {
    if(input$curr_mz %in% colnames(mSet$dataSet$norm)){
      my_selection$mz <- input$curr_mz   
      plotmanager$make <- "summary"
    }
  })
  
  shiny::observe({
    if(my_selection$struct != "" & input$tab_iden_2 == "molmz"){
      if(!mSet$metshiParams$prematched){
        MetaboShiny::metshiAlert("Please perform pre-matching first to enable this feature!")
        return(NULL)
      }else{
        
        if(lcl$prev_struct != my_selection$struct){
          rev_matches = MetaboShiny::get_prematches(who = my_selection$struct,
                                                    what = "con.structure",
                                                    patdb = lcl$paths$patdb,
                                                    showadd = c(),
                                                    showiso = c(),
                                                    showdb = c())  
          lcl$prev_struct <<- my_selection$struct
          if(nrow(rev_matches) == 0){
            shown_matches$reverse <- data.table::data.table()
            return(NULL)
          }else{
            pieinfo$db <- reshape::melt(table(rev_matches$source))
            pieinfo$add <- reshape::melt(table(rev_matches$adduct))
            pieinfo$iso <- reshape::melt(table(rev_matches$isocat))
          }
        }
        mzMode =if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
        rev_matches = MetaboShiny::get_prematches(who = my_selection$struct,
                                                  what = "con.structure",
                                                  patdb = lcl$paths$patdb,
                                                  showadd = result_filters$add[[mzMode]],
                                                  showiso = result_filters$iso,
                                                  showdb = result_filters$db)
        if(nrow(rev_matches)>0){
          shown_matches$reverse <- unique(rev_matches[,c("query_mz", "adduct", "%iso", "dppm")])
        }else{
          shown_matches$reverse <- data.table::data.table()
        }
      }
    } 
  })
  
  # pie charts
  lapply(c("add", "iso", "db"), function(which_pie){
    output[[paste0("match_pie_", which_pie)]] <- plotly::renderPlotly({
      
      pievec = pieinfo[[which_pie]]
      
      if(nrow(pievec) > 0){
        
        m <- list(
          l = 0,
          r = 0,
          b = 30,
          t = 30)
        
        pulls = rep(0, nrow(pievec))
        lines = rep(1, nrow(pievec))
        
        if(!is.null(pievec)){
          if(which_pie == "add"){
            mzMode = if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
            targets = result_filters$add[[mzMode]]
          }else{
            targets = result_filters[[which_pie]]
          }
          pulls[which(as.character(pievec$Var.1) %in% targets)] <- 0.15  
          lines[which(as.character(pievec$Var.1) %in% targets)] <- 4
        }
        
        myCols <- gbl$functions$color.functions[[lcl$aes$spectrum]](n = nrow(pievec))
        
        if(length(pievec)>0){
          p = plotly::plot_ly(pievec, labels = ~Var.1, 
                              values = ~value, size=~value*10, type = 'pie',
                              textposition = 'inside',
                              textinfo = 'label+percent',
                              insidetextfont = list(colors = ggdark::invert_color(myCols)),
                              hoverinfo = 'text',
                              pull = pulls,
                              text = ~paste0(Var.1, ": ", value, ' matches'),
                              marker = list(colors = myCols,
                                            line = list(color = "gray", 
                                                        width = lines)),
                              #The 'pull' attribute can also be used to create space between the sectors
                              showlegend = FALSE) %>%
            plotly::layout(autosize = T, margin = m,
                           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))  
        }
        return(p)          
      }
    })
  })
  
  # ===== UI SWITCHER ====
  
  shiny::observeEvent(input$heattable,{
    if(!is.null(input$overview)){
      if(input$overview == "heatmap"){
        statsmanager$calculate <- "heatmap"
        plotmanager$make <- "heatmap"
        uimanager$refresh <- "heatmap"
      }  
    }
  })
  
  shiny::observeEvent(input$network_style, {
    if("network" %in% names(mSet$analSet)){
      plotmanager$make <- "network"
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
      list("overview", "volcano"),#8
      list("overview", "venn"),#9
      list("overview", "heatmap"),#10
      list("overview", "power"),#11
      list("permz", "tt"),#12
      list("permz", "fc"),#13
      list("dimred", "tsne"),#14
      list("permz", "corr"),#15
      list("overview", "enrich")#16
      # list("main", "dimred"),#17
      # list("main", "overview"),#18
      # list("main", "ml"),#19
      # list("main", "permz")#20
    )
    
    # check mode of interface (depends on timeseries /yes/no and bivariate/multivariate)
    # then show the relevent tabs
    # TODO: enable multivariate time series analysis
    if(is.null(interface$mode)){
      show.tabs <- hide.tabs[1]
    }else if(interface$mode == '1fb'){
      show.tabs <- hide.tabs[c(1,2,3,7,8,9,10,11,12,13,14,15,16)]
      shiny::updateSelectInput(session, "ml_method",
                               selected = "rf",
                               choices = as.list(gbl$constants$ml.models))
    }else if(interface$mode == '1fm'){
      show.tabs <- hide.tabs[c(1,2,3,6,7,9,10,11,14,15,16)]
      shiny::updateSelectInput(session, "ml_method",
                               selected = "rf",
                               choices = as.list(setdiff(gbl$constants$ml.models,
                                                         gbl$constants$ml.twoonly)))
    }else if(interface$mode == '2f'){
      show.tabs <- hide.tabs[c(1,2,4,6,9,10,11,14,16)]
    }else if(interface$mode == 't1f'){
      show.tabs = hide.tabs[c(1,2,4,5,6,9,10,11,14,16)]
    }else if(interface$mode == 't'){
      show.tabs = hide.tabs[c(1,2,5,6,7,9,10,11,14,15,16)]
    }else{
      show.tabs <- hide.tabs[1]
    }
    
    # hide all the tabs to begin with
    if(length(show.tabs) > 1){
      for(bigtab in c("dimred", "permz", "overview", "ml")){
        shiny::showTab("statistics", bigtab)  
      }
      for(bigtab in c("search", "plot_aes", "switchset", "metadata")){
        shiny::showTab("anal_sidebar", bigtab)  
      }
      shiny::hideTab("anal_sidebar", "start")
      shiny::updateTabsetPanel(session, "anal_sidebar", "switchset")
    }else{
      for(bigtab in c("dimred", "permz", "overview", "ml")){
        shiny::hideTab("statistics", bigtab)  
      }
      for(bigtab in c("search", "plot_aes", "switchset", "metadata")){
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
  
  # print current compound in sidebar
  observe({
    shinyWidgets::updatePickerInput(session, 
                                    "curr_mz", 
                                    selected = my_selection$mz)
  })
  
  
  # make miniplot for sidebar with current compound
  output$curr_plot <- plotly::renderPlotly({
    if(my_selection$mz != ""){
      MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, 
                                 cols = lcl$aes$mycols, cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                 styles = input$ggplot_sum_style,
                                 add_stats = input$ggplot_sum_stats, 
                                 color.fac = input$col_var,
                                 text.fac = input$txt_var,
                                 plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                 font = lcl$aes$font)
    }
  })
  
  shiny::observeEvent(input$undo_match_filt, {
    result_filters$add <- list(positive = c(), negative = c())
    result_filters$db <- c()
    result_filters$iso <- c()
    search$go <- T
  }) 
  
  output$curr_add <- shiny::renderText({
    paste0(unlist(result_filters$add), collapse=", ")
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
      mzMode = if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
      for(mode in mzMode){
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
  # shiny::observe({
  #   if(length(input$curr_mz) > 0){
  #     if(input$curr_mz != my_selection$mz){
  #       print(input$curr_mz %in% colnames(mSet$dataSet$norm))
  #       shinyWidgets::updatePickerInput(session,
  #                                       "curr_mz",
  #                                       selected = my_selection$mz)
  #     }      
  #   }
  # })
  
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
  
  shiny::observeEvent(input$ncores, {
    if(!is.null(session_cl)){
      shiny::showNotification("Stopping threads...")
      parallel::stopCluster(session_cl)
    }
    shiny::showNotification("Starting new threads...")

    session_cl <<- parallel::makeCluster(input$ncores)#,outfile="")#,setup_strategy = "sequential") # leave 1 core for general use and 1 core for shiny session
    # send specific functions/packages to other threads
    parallel::clusterEvalQ(session_cl, {
      library(data.table)
      library(iterators)
      library(MetaboShiny)
      library(MetaDBparse)})
    MetaboShiny::setOption(lcl$paths$opt.loc, "cores", input$ncores)
  })
  
  shiny::observeEvent(input$set_api, {
    MetaboShiny::setOption(lcl$paths$opt.loc, "apikey", input$apikey)
    lcl$apikey <<- input$apikey
    output$api_set <- shiny::renderText("key saved!")
  })
  
  shiny::observeEvent(input$omit_unknown, {
    new.val = if(input$omit_unknown) "yes" else "no"
    MetaboShiny::setOption(lcl$paths$opt.loc, "omit_unknown", new.val)
  })
  
  shiny::observeEvent(input$ml_top_x, {
    if(!is.null(mSet)){
      plotmanager$make <- "ml"
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
      tablemanager$make <- "ml"
      plotmanager$make <- "ml"
    }
  },ignoreNULL = T, ignoreInit = T)
  
  shiny::observeEvent(input$save_mset, {
    # save mset
    if(!is.null(mSet)){
      shiny::withProgress({
        fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
        if(exists("mSet")){
          save(mSet, file = fn)
        }
      })  
    }
  })
  
  shiny::observeEvent(input$stats_type,{
    if(!is.null(mSet)){
      uimanager$refresh <- "statspicker"
    }
  })
  
  # for(store_item in names(mSet$storage)){
  #   print(store_item)
  #   print(object.size(mSet$storage[[store_item]]), unit="Mb")
  # }
  
  shiny::observeEvent(input$load_mset, {
    # load mset
    shiny::withProgress({
      fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
      shiny::showNotification(paste0("Loading existing file: ", fn))
      if(file.exists(fn)){
        mSet <- tryCatch({
          load(fn)
          if(is.list(mSet)){
            metshiAlert("Old save selected! Please re-import data to use in v2.0. Conversion is not possible.")
          }
          mSet <- NULL
          gc()
        },
        error = function(cond){
          mSet <- qs::qload(fn)
          mSet
        })
        if(is.list(mSet)){
          mSet$dataSet$combined.method <- TRUE # FC fix
          mSet <<- mSet
          opts <- MetaboShiny::getOptions(lcl$paths$opt.loc)
          lcl$proj_name <<- opts$proj_name
          lcl$paths$proj_dir <<- file.path(lcl$paths$work_dir, lcl$proj_name)
          lcl$paths$patdb <<- file.path(lcl$paths$proj_dir, paste0(opts$proj_name, ".db"))
          lcl$paths$csv_loc <<- file.path(lcl$paths$proj_dir, paste0(opts$proj_name, ".csv"))
          
          shiny::updateCheckboxInput(session,
                                     "paired",
                                     value = mSet$dataSet$paired)
          
          uimanager$refresh <- c("general","statspicker",if("adducts" %in% names(opts)) "adds" else NULL, "ml")
          plotmanager$make <- "general"  
        }
      }
    })
  })
  
  shiny::observeEvent(input$debug_metshi, {
    assign("lcl", lcl, envir = .GlobalEnv)
    assign("mSet", mSet, envir = .GlobalEnv)
    assign("input", shiny::isolate(shiny::reactiveValuesToList(input)), envir = .GlobalEnv)
    assign("enrich", shiny::isolate(shiny::reactiveValuesToList(enrich)), envir = .GlobalEnv)
    assign("shown_matches", shiny::isolate(shiny::reactiveValuesToList(shown_matches)), envir = .GlobalEnv)
    assign("my_selection", shiny::isolate(shiny::reactiveValuesToList(my_selection)), envir = .GlobalEnv)
    assign("browse_content",  shiny::isolate(shiny::reactiveValuesToList(browse_content)), envir = .GlobalEnv)
    assign("pieinfo",  shiny::isolate(shiny::reactiveValuesToList(pieinfo)), envir = .GlobalEnv)
    assign("result_filters",  shiny::isolate(shiny::reactiveValuesToList(result_filters)), envir = .GlobalEnv)
    assign("report_yes",  shiny::isolate(shiny::reactiveValuesToList(report_yes)), envir = .GlobalEnv)
    assign("venn_yes",  shiny::isolate(shiny::reactiveValuesToList(venn_yes)), envir = .GlobalEnv)
  })
  
  shiny::observeEvent(input$ml_train_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$vectors$ml_train <<- c(input$subset_var,
                               input$subset_group)
    output$ml_train_ss <- shiny::renderText(subset.name)
  })
  
  shiny::observeEvent(input$reset_ml_train, {
    subset.name <- "all"
    lcl$vectors$ml_train <<- NULL
    output$ml_train_ss <- shiny::renderText(subset.name)
  })
  
  shiny::observeEvent(input$ml_test_ss, {
    keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$vectors$ml_test <<- c(input$subset_var, input$subset_group)
    output$ml_test_ss <- shiny::renderText(subset.name)
  })
  
  shiny::observeEvent(input$reset_ml_test, {
    subset.name <- "all"
    lcl$vectors$ml_test <<- NULL
    output$ml_test_ss <- shiny::renderText(subset.name)
  })
  
  shiny::observeEvent(input$subset_var, {
    lvls = levels(as.factor(mSet$dataSet$covars[[input$subset_var]]))
    shiny::updateSelectizeInput(session, "subset_group", choices = lvls)
  },ignoreInit = T)
  
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
  
  observeEvent(c(input$heatsign, input$heatlimits),{
    if(!is.null(mSet)){
      if("heatmap" %in% names(mSet$analSet)){
        plotmanager$make <- "heatmap"  
      }  
    }
  })
  
  observeEvent(input$save_fav_adducts,{
    # save to options
    pasted_adducts <- paste0(input$fav_adducts, collapse = "&")
    MetaboShiny::setOption(lcl$paths$opt.loc, 
                           "adducts", 
                           pasted_adducts)
    uimanager$refresh <- "adducts"
  })
  
  observeEvent(input$nav_general, {
    if(!is.null(mSet)){
      if(input$nav_general == "report"){
        statsmanager$calculate <- "vennrich"
        tablemanager$make <- "vennrich"
        uimanager$refresh <- "vennrich"
      }}
  })
  
  observeEvent(input$statistics, { 
    if(!is.null(mSet)){
      if(!is.null(input$statistics)){
        uimanager$refresh <- input$statistics
        if(input$statistics %in% c("venn", "enrich", "heatmap", "network")){
          statsmanager$calculate <- "vennrich"
          tablemanager$make <- "vennrich"
          uimanager$refresh <- "vennrich"
        }
        if(input$statistics %in% names(mSet$analSet)){
          shinyjs::show(selector = paste0("div.panel[value=collapse_", input$statistics, "_tables]"))
          tablemanager$make <- input$statistics
          shinyjs::show(selector = paste0("div.panel[value=collapse_", input$statistics, "_plots]"))
          plotmanager$make <- input$statistics
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
                "pca", "tsne", "tt", "aov",
                "fc", "volcano", "heatmap", 
                "meba", "asca", "corr", 
                "enrich", "network")
  
  lapply(analyses, function(an){
    shiny::observeEvent(input[[paste0("do_", an)]], {
      try({
        statsmanager$calculate <- an
        shinyjs::show(selector = paste0("div.panel[value=collapse_", an, "_tables]"))
        tablemanager$make <- an
        shinyjs::show(selector = paste0("div.panel[value=collapse_", an, "_plots]"))
        plotmanager$make <- an
        uimanager$refresh <- an
      })
    })    
  })
  
  output$manual_search <- renderUI({
    if(search_button$on){
      tags$button(
        id = "search_mz",
        class = "btn btn-default action-button",
        img(src = "detective.png",
            height = "50px")
      )
    }else{
      fluidRow(align="center", 
               shiny::img(src = "pawprint.png",height = "50px"),
               br(),
               tags$h2("pre-matched")
               )
    }
  })
  
  # ==== LOAD LOGIN UI ====
  
  # init all observer
  for(fp in list.files("reactive", full.names = T)){
    source(fp, local = T)
  }  
  
  observeEvent(input$quit_metshi, {
    if(!is.null(mSet)){
      shinyWidgets::confirmSweetAlert(
        session = session,
        inputId = "save_exit",
        text = tags$div(
          tags$b("Click upper right ", icon("times"), " button to cancel."),br(),
          shiny::img(#class = "rotategem", 
            src = "gemmy_rainbow.png", 
            width = "70px", height = "70px"),
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
      shiny::withProgress({
        fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
        if(exists("mSet")){
          qs::qsave(mSet, fn)
        }
      })
    }
    shiny::stopApp()
    shinyjs::js$closeWindow()
  },ignoreNULL = T)
  
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
                             #shiny::helpText(pkgRes$msg)
                             #shiny::verbatimTextOutput(textID)
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
        MetaboShiny::metshiAlert(shiny::tagList(h3(headerMsg),
                                                hr(),
                                                alertContent), 
                                 session = session,
                                 title = "Notification")
      }
      
    }  
  },silent = T)
  
  onStop(function() {
    print("Closing MetaboShiny ~ヾ(＾∇＾)")
    if(!is.null(session_cl)){
      parallel::stopCluster(session_cl)
    }
    session_cl <<- NULL
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if(all(file.remove(rmv))) NULL
  })
}