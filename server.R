# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

shinyServer(function(input, output, session) {
  
  # ================================= DEFAULTS ===================================
  
  source('./backend/scripts/joanna/shiny_general.R')

  shinyOptions(progress.style="old")
  
  parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
    "mape",
    "flattenlist"
  ))
  
  parallel::clusterEvalQ(session_cl, library(data.table))

  # - - - - - - - - - - - - -
  
  spinnyimg <- reactiveVal("www/electron.png")
  
  lapply(global$constants$default.text, FUN=function(default){
    output[[default$name]] = renderText(default$text)
  })
  
  lapply(global$constants$images, FUN=function(image){
    output[[image$name]] <- renderImage({
      filename <- normalizePath(image$path)
      # Return a list containing the filename and alt text
      list(src = filename, 
           width = image$dimensions[1],
           height = image$dimensions[2])
    }, deleteFile = FALSE)
  })
  
  output$spinny <- renderText({spinnyimg()})
  
  output$taskbar_image <- renderImage({
    list(src = file.path(getwd(), 
                         "www", 
                         getOptions("user_options.txt")$taskbar_image), 
         width = 120,
         height = 120,
         style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;")
  }, deleteFile = FALSE)
  
  output$colorPickers <- renderUI({
    lapply(c(1:global$constants$max.cols), function(i) {
      print(i)
      colourpicker::colourInput(inputId = paste("col", i, sep="_"),
                                label = paste("Choose colour", i),
                                value = global$vectors$mycols[i],
                                allowTransparent = F)
    })
  })
  
  observe({
    values <- unlist(lapply(c(1:global$constants$max.cols), function(i) {
      input[[paste("col", i, sep="_")]]
    }))
    if(!any(is.null(values))){
      set.col.map("user_options.txt", values)
      global$vectors$mycols <<- get.col.map("user_options.txt")
    }
  })
  
  # ===== STARTUP SETTINGS =====
  
  datamanager <- reactiveValues()
  
  observe({
    if(exists("mSet")){
      datamanager$mset_present = TRUE
      # - - -
      if(mSet$dataSet$cls.num <= 1){
        interface$mode <- NULL } 
      else if(mSet$dataSet$cls.num == 2){
        interface$mode <- "bivar"}
      else{
        interface$mode <- "multivar"}
    }else{
      datamanager$mset_present = FALSE
    }
  })
  
  observe({
    if(datamanager$mset_present){
      # update select input bars
      updateSelectInput(session, "first_var", selected = mSet$dataSet$cls.name, choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
      updateSelectInput(session, "second_var", choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
      updateSelectInput(session, "subset_var", choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
      # timecourse button
      if(all(grepl(pattern = "_T\\d", x = rownames(mSet$dataSet$norm)))){
        timebutton$status <- "on"
      }else{
        timebutton$status <- "off"
      }
      if(mSet$dataSet$cls.num == 2 ){
        heatbutton$status <- "ttfc"
      }else{
        heatbutton$status <- NULL
      }
    }else{
      timebutton$status <- "off"
      heatbutton$status <- "asmb"
    }
  })

  # ===== VARIABLE SWITCHER ====
  
  heatbutton <- reactiveValues(status = "ttfc")
  
  timebutton <- reactiveValues(status = "off")
  
  output$heatbutton <- renderUI({
    print(heatbutton$status)
    if(is.null(heatbutton$status)){
      NULL
    }else{
      print(heatbutton$status)
      switch(heatbutton$status,
             asmb = switchButton(inputId = "heatmode",
                                 label = "Use data from:", 
                                 value = TRUE, col = "BW", type = "ASMB"),
             ttfc = switchButton(inputId = "heatmode",
                              label = "Use data from:", 
                              value = TRUE, col = "BW", type = "TTFC")
             )
    }
  })
  
  output$timebutton <- renderUI({
    print(timebutton$status)
    if (is.null(timebutton$status)) {
      NULL
    } else{
      switch(timebutton$status, 
             off = NULL,
             on = switchButton(inputId = "timecourse_trigger",
                               label = "Toggle time course mode?", 
                               value = FALSE, col = "BW", type = "YN"))
    }
  })
  
  observeEvent(input$timecourse_trigger, {
    
    print(input$timecourse_trigger)
    
    if(!("storage" %in% names(mSet))){
      mSet$storage <<- list()
    }
    
    if(input$timecourse_trigger){
      # change to timecourse mode
      # save previous mset
      mSet$storage[[paste0(as.character(levels(mSet$dataSet$cls)), collapse="-")]] <<- mSet$dataSet$analSet
      
      # adjust mset
      SetDesignType(mSet, "time")
      
      # facs
      facA <- as.factor(mSet$dataSet$covars[,mSet$dataSet$cls.name, with=F][[1]])
      facB <- mSet$dataSet$covars[,"time"][[1]]
      
      mSet$dataSet$exp.fac <<- as.factor(facA)
      mSet$dataSet$time.fac <<- as.factor(facB)
      
      interface$mode <- "time"
      
      heatbutton$status <- "asmb"
      #mSet$analSet <<- NULL
      
    }else{
      # change back to normal mode
      # save previous analyses (should be usable in venn diagram later)
      mSet$storage[[paste0("(timecourse)", paste0(as.character(levels(mSet$dataSet$cls)), collapse="-"))]] <<- mSet$dataSet$analSet
      # - - - - - - - - - - - -
      mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,mSet$dataSet$cls.name, with=F][[1]])
      mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
      # remove old analSet
      
      heatbutton$status <- "ttfc"
      #mSet$analSet <<- NULL
      
      # reset interface
      if(mSet$dataSet$cls.num <= 1){
        interface$mode <- NULL } 
      else if(mSet$dataSet$cls.num == 2){
        interface$mode <- "bivar"}
      else{
        interface$mode <- "multivar"}
    }
  }, ignoreInit = TRUE)
  
  interface <- reactiveValues()
  
  observeEvent(input$change_cls, {
    print(input$change_cls)
    print(input$first_var)
    
    if(!("storage" %in% names(mSet))){
      mSet$storage <<- list()
    }
    
    # save previous analyses (should be usable in venn diagram later)
    mSet$storage[[paste0(as.character(levels(mSet$dataSet$cls)), collapse="-")]] <<- mSet$dataSet$analSet
    # - - - - - - - - - - - -
    mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,input$first_var, with=F][[1]])
    mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
    # remove old analSet
    mSet$analSet <<- NULL
    mSet$dataSet$cls.name <<- input$first_var
    # reset interface
    if(mSet$dataSet$cls.num <= 1){
      interface$mode <- NULL } 
    else if(mSet$dataSet$cls.num == 2){
      interface$mode <- "bivar"}
    else{
      interface$mode <- "multivar"}
  })
  
  # ===== UI SWITCHER ====
  
  observe({
    
    hide.tabs <- c("inf", "pca", "plsda", "tt", "fc", "aov", "meba", "asca", "ml", "volc", "heatmap", "enrich")
    
    print(interface$mode)
    
    if (is.null(interface$mode)) {
      show.tabs <- c("inf")
    } else if(interface$mode == 'multivar'){ 
      show.tabs <- c("pca", "aov", "heatmap")
    } else if(interface$mode == 'bivar'){  
      show.tabs <- c("pca", "plsda", "tt", "fc", "volc", "heatmap", "ml")
    } else if(interface$mode == 'time'){
      show.tabs <- c("pca", "asca", "meba", "heatmap")
    }
    else{
      show.tabs <- c("inf")
    }
    # - show/hide - 
    for(tab in hide.tabs){
      print(tab)
      hideTab(inputId = "statistics", tab, session = session)
    }
    i=1
    for(tab in show.tabs){
      print(tab)
      showTab(inputId = "statistics", tab, select = ifelse(i==1, TRUE, FALSE), session = session)
      i = i + 1
    }
  })
  
  # -----------------
  
  
  output$pos_add_tab <-DT::renderDataTable({
    # -------------
    DT::datatable(global$vectors$pos_adducts,
                  selection = list(mode = 'multiple', selected = c(1:3, nrow(global$vectors$pos_adducts)), target = 'row'),
                  options = list(pageLength = 5, dom = 'tp'), 
                  rownames = F)
  })
  
  output$neg_add_tab <-DT::renderDataTable({
    # -------------
    DT::datatable(global$vectors$neg_adducts,
                  selection = list(mode = 'multiple', selected = c(1, 2, 14, 15, nrow(global$vectors$neg_adducts)), target = 'row'),
                  options = list(pageLength = 5, dom = 'tp'), 
                  rownames = F)
  })
  
  observeEvent(input$sel_all_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', selected = c(1:nrow(global$vectors$pos_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$neg_adducts,
                    selection = list(mode = 'multiple', selected = c(1:nrow(global$vectors$neg_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  observeEvent(input$sel_no_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', selected = c(0), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$neg_adducts,
                    selection = list(mode = 'multiple', selected = c(0), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })

  observeEvent(input$sel_comm_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', selected = c(1:3, nrow(global$vectors$pos_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$neg_adducts,
                    selection = list(mode = 'multiple', selected = c(1, 2, 14:15, nrow(global$vectors$neg_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  observeEvent(input$nav_general, {
    # - - - - - -
    pkg_tbl <- get.package.table() #TODO: sort by 'No' first!! (ascending?) - or translate to numeric factor first?
    output$package_tab <- DT::renderDataTable({
    # - - - - - -
      DT::datatable(pkg_tbl,
                    selection = 'none',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(10, 20, 30), pageLength = 10), 
                    rownames = F)
    })
  }
  )
  
  observeEvent(input$update_packages, {
    pacman::p_load(char = global$constants$packages, update = T, character.only = T)
    # - - refresh package list - - 
    pkg_tbl <- get.package.table() 
    output$package_tab <- DT::renderDataTable({
      # - - - - - - - - - - -
      DT::datatable(pkg_tbl,
                    selection = 'none',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(10, 20, 30), pageLength = 10), 
                    rownames = F)
    })
  })

  observe({
    shinyFileChoose(input, 'outlist_pos', roots=global$paths$volumes, filetypes=c('csv'))
    shinyFileChoose(input, 'outlist_neg', roots=global$paths$volumes, filetypes=c('csv'))
    shinyFileChoose(input, 'excel', roots=global$paths$volumes, filetypes=c('xls', 'xlsm', 'xlsx'))
    shinyFileChoose(input, 'database', roots=global$paths$volumes, filetypes=c('sqlite3', 'db', 'sqlite'))
    shinyFileChoose(input, 'taskbar_image_path', roots=global$paths$volumes, filetypes=c('png', 'jpg', 'jpeg', 'bmp'))
  })
  
  observe({
    # - - - - 
    if(!is.list(input$taskbar_image_path)) return()
    img_path <- parseFilePaths(global$paths$volumes, input$taskbar_image_path)$datapath
    new_path <- file.path(getwd(), "www", basename(img_path))
    
    if(img_path != new_path) file.copy(img_path, new_path, overwrite = T)
    # - - -
    output$taskbar_image <- renderImage({
      list(src = new_path, 
           width = 120,
           height = 120,
           style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;")
    }, deleteFile = FALSE)
    # - - -
    setOption('user_options.txt', 'taskbar_image', basename(new_path))
  })
  
  observe({  
    shinyDirChoose(input, "get_db_dir", 
                   roots=global$paths$volumes, 
                   session = session)
     
    test <<- input$get_db_dir
    
    if(typeof(input$get_db_dir) != "list") return()

    given_dir <- parseDirPath(global$paths$volumes, 
                              input$get_db_dir)
    if(is.null(given_dir)) return()
    # - - connect - -
    setOption("user_options.txt", "db_dir", given_dir)
    output$curr_db_dir <- renderText({getOptions('user_options.txt')$db_dir})
  })
  
  observe({
    shinyDirChoose(input, "get_work_dir",
                   roots = global$paths$volumes,
                   session = session)

    if(typeof(input$get_work_dir) != "list") return()
    
    given_dir <- parseDirPath(global$paths$volumes,
                              input$get_work_dir)
    if(is.null(given_dir)) return()
    setOption("user_options.txt", "work_dir", given_dir)
    output$curr_exp_dir <- renderText({getOptions('user_options.txt')$work_dir})
  })
  
  observeEvent(input$set_proj_name, {
    proj_name <<- input$proj_name
    if(proj_name == "") return(NULL)
    global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(proj_name,".db", sep=""))
    # --- connect ---
    setOption("user_options.txt", "proj_name", proj_name)
    output$proj_name <<- renderText(proj_name)
    global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
    
  })
  
  observeEvent(input$set_ppm, {
    ppm <<- input$ppm
    output$ppm <<- renderText(ppm)
    # --- connect ---
    setOption("user_options.txt", "ppm", ppm)
  })
  
  observeEvent(input$color_ramp,{
    output$ramp_plot <- plotly::renderPlotly({

      global$functions$color.functions[[getOptions("user_options.txt")$gspec]] <<- global$functions$color.functions[[input$color_ramp]]
      setOption("user_options.txt", "gspec", input$color_ramp)
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        titlefont = list(size = 20)
      )
      
      # --- re-render --- 
      #TODO: SHOULD BE GENERAL PURPOSE FUNCTION w/ listeners
      plotly::plot_ly(z = volcano, 
                      colors = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](100), 
                      type = "heatmap",
                      showscale=FALSE)  %>%
        layout(xaxis = ax, yaxis = ax)
    })
  })
  
  observeEvent(input$ggplot_theme,{
    
    setOption("user_options.txt", "gtheme", input$ggplot_theme)
    
    output$ggplot_theme_example <- renderPlot({
      p <- ggplot(mtcars) + geom_boxplot(aes(x = wt, y = mpg,
                                           colour = factor(gear)))
      p + global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]]()
    })
  })
  
  observeEvent(input$change_css, {
    # - - - - - - - - - - - - - - -
    setOption("user_options.txt", "col1", input$bar.col.1)
    setOption("user_options.txt", "col2", input$bar.col.2)
    setOption("user_options.txt", "col3", input$bar.col.3)
    setOption("user_options.txt", "col4", input$bar.col.4)
    
    setOption("user_options.txt", "font1", input$font.1)
    setOption("user_options.txt", "font2", input$font.2)
    setOption("user_options.txt", "font3", input$font.3)
    setOption("user_options.txt", "font4", input$font.4)
    
    setOption("user_options.txt", "size1", input$size.1)
    setOption("user_options.txt", "size2", input$size.2)
    setOption("user_options.txt", "size3", input$size.3)
    setOption("user_options.txt", "size4", input$size.4)
  })

  # --- adduct table editing ---
  
  values = reactiveValues()
  
  observeEvent(input$import_adducts, {
    req(input$add_tab)
    # ----------------
    DF = fread(input$add_tab$datapath)
    output$adduct_tab <- rhandsontable::renderRHandsontable({
      if (!is.null(DF))
        rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
    })
    output$adduct_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  
  adduct_tab_data <- reactive({
    if (!is.null(input$adduct_tab)) {
      DF = rhandsontable::hot_to_r(input$adduct_tab)
    } else {
      if (is.null(values[["DF"]]))
        DF = adducts
      else
        DF = values[["DF"]]
    }
    values[["DF"]] = DF
    # ---------------
    DF
  })
  
  output$adduct_tab <- rhandsontable::renderRHandsontable({
    DF = adduct_tab_data()
    if (!is.null(DF))
      rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
  })
  
  observe({
    DF = adduct_tab_data()
    shinyFileSave(input, "save_adducts", roots = c(home = '~'), session=session)
    fileinfo <- parseSavePath(roots = c(home = '~'), input$save_adducts)
    if (nrow(fileinfo) > 0) {
      switch(fileinfo$type,
             csv = fwrite(file = fileinfo$datapath, x = DF)
      )}
  })
  
  # ======================== DB CHECK ============================
  
  # --- check for db files ---
  
  lapply(global$vectors$db_list, FUN=function(db){
    observeEvent(input[[paste0("check_", db)]],{
      db_folder_files <- list.files(getOptions("user_options.txt")$db_dir)
      is.present <- paste0(db, ".full.db") %in% db_folder_files
      check_pic <- if(is.present) "yes.png" else "no.png"
      output[[paste0(db,"_check")]] <- renderImage({
        # When input$n is 3, filename is ./images/image3.jpeg
        filename <- normalizePath(file.path('www', check_pic))
        # Return a list containing the filename and alt text
        list(src = filename, width = 70,
             height = 70)
      }, deleteFile = FALSE)
    })
  })
  
  # --- build db ---
  
  lapply(global$vectors$db_list, FUN=function(db){
    observeEvent(input[[paste0("build_", db)]], {
      # ---------------------------
      library(RCurl)
      library(XML)
      library(SPARQL)
      # ---------------------------
      withProgress({
        
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
        build.base.db(db,
                      outfolder = getOptions("user_options.txt")$db_dir, 
                      cl = session_cl)
        shiny::setProgress(session = session, 0.5)
        build.extended.db(db, 
                          outfolder = getOptions("user_options.txt")$db_dir,
                          adduct.table = adducts, 
                          cl = session_cl, 
                          fetch.limit = 500) #TODO: figure out the optimal fetch limit...
      })
    })
  })
  
  observeEvent(input$score_iso, {

    req(input$iso_score_method)
    
    if(!data.table::is.data.table(global$tables$last_matches)) return(NULL)
    
    if("score" %in% colnames(global$tables$last_matches)){
      global$tables$last_matches <<- global$tables$last_matches[,-"score"]
    }
    
    # - - - 
    
    score_table <- score.isos(global$paths$patdb, method=input$iso_score_method, inshiny=T) 
    
    global$tables$last_matches <<- global$tables$last_matches[score_table, on = c("baseformula", "adduct")]
    
    output$match_tab <-DT::renderDataTable({
      DT::datatable(global$tables$last_matches[,-c("description","structure", "baseformula", "dppm")],
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })  
  })
  
  # ================== DATA IMPORT ===========================
  
  observeEvent(input$create_db,{
    
    global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name, ".db"))
    
    withProgress({
      
      shiny::setProgress(session=session, value= .1)
      
      switch(input$new_proj,
             `From DB` = {
               
               req(input$database, input$excel)
               
               db_path <- parseFilePaths(global$paths$volumes, input$database)$datapath
               excel_path <- parseFilePaths(global$paths$volumes, input$excel)$datapath
               
               file.copy(db_path, global$paths$patdb, overwrite = T)
               
               shiny::setProgress(session=session, value= .30)
               
               exp_vars <- load.excel(excel_path, global$paths$patdb)
               
               shiny::setProgress(session=session, value= .60)
               
             },
             `From CSV` = {
               
               req(input$outlist_neg, input$outlist_pos, input$excel)
               
               build.pat.db(global$paths$patdb,
                            ppm = ppm,
                            pospath = parseFilePaths(global$paths$volumes, input$outlist_pos)$datapath,
                            negpath = parseFilePaths(global$paths$volumes, input$outlist_neg)$datapath,
                            overwrite = T)
               
               shiny::setProgress(session=session, value= .95,message = "Adding excel sheets to database...")
               
               exp_vars <<- load.excel(parseFilePaths(global$paths$volumes, input$excel)$datapath, global$paths$patdb)}
             )
      })
  })
  
  observeEvent(input$import_db, {
    req(input$pat_db)
    global$paths$patdb <<- input$pat_db$datapath
    output$db_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  
  observeEvent(input$import_csv, {
    req(input$pat_csv)
    global$paths$csv_loc <<- input$pat_csv$datapath
    
    output$csv_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  
  # ==================== CREATE CSV =======================
  
  observeEvent(input$create_csv, {

    withProgress({
      
      shiny::setProgress(session = session, value= 1/4)

      tbl <- get.csv(global$paths$patdb,
                     group_adducts = if(length(global$vectors$add_search_list) == 0) F else T,
                     groupfac = input$group_by,
                     which_dbs = global$vectors$add_search_list,
                     which_adducts = selected_adduct_list
      )
      

      if(any(duplicated(tbl$sample))){
        tbl$sample <- paste0(tbl$sample, "_T", tbl$time)
        show.times = T
      }else{
        show.times = F
      }
      
      shiny::setProgress(session=session, value= 2/4)
      global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
      fwrite(tbl, global$paths$csv_loc, sep="\t")

      as.numi <- as.numeric(colnames(tbl)[1:100])
      
      exp.vars <- which(is.na(as.numi))
      
      shiny::setProgress(session=session, value= 3/4)
      
      output$csv_tab <-DT::renderDataTable({
        
      overview_tab <- if(show.times){
        t(data.table(keep.rownames = F,
                     Identifiers = ncol(tbl) - length(exp.vars),
                     Samples = nrow(tbl),
                     Times = length(unique(tbl$time))
        )) 
      }else{
        t(data.table(keep.rownames = F,
                     Identifiers = ncol(tbl) - length(exp.vars),
                     Samples = nrow(tbl)
        )) 
      }
      
        # --- render ---
        DT::datatable(overview_tab, 
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE))
      })
    })
  })
  
  # ===================== METABOSTART ========================
  
  observeEvent(input$check_csv, {
    # ----------------------
    
    csv <- fread(global$paths$csv_loc, 
                 header = T)
    
    as.numi <- as.numeric(colnames(csv)[1:100])
    
    exp.vars <- which(is.na(as.numi))
    
    opts <<- colnames(csv[,..exp.vars])
    
    batch <<- which(sapply(exp.vars, function(x) length(unique(csv[,..x][[1]])) < nrow(csv)))
    
    numi <<- which(sapply(exp.vars, function(x) is.numeric(csv[,..x][[1]])))

    updateSelectInput(session, "samp_var",
                      choices = opts[numi])
    updateSelectizeInput(session, "batch_var",
                         choices = opts[batch],
                         options = list(maxItems = 3L - (length(input$batch_var)))
    )
  })
  
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
  
  observeEvent(input$check_csv, {
    req(global$paths$csv_loc)
    switch(input$norm_type,
           ProbNorm=updateSelectInput(session, "ref_var",
                                      choices = get_ref_vars(fac = "label") # please add options for different times later, not difficult
           ),
           CompNorm=updateSelectInput(session, "ref_var",
                                      choices = get_ref_cpds() # please add options for different times later, not difficult
           ))
    # get excel table stuff.
    
  })
  
  output$ref_select <- renderUI({ref.selector()})
  
  observeEvent(input$initialize,{
    req(input$filt_type)
    req(input$norm_type)
    req(input$trans_type)
    req(input$scale_type)
    withProgress({

      shiny::setProgress(session=session, value= .1)

      # input <- list(batch_var = c("batch", "country"),
      #               exp_type = "stat",
      #               perc_limit = .99,
      #               filt_type = "none",
      #               miss_type = "rf",
      #               norm_type = "SumNorm",
      #               trans_type = "LogNorm",
      #               scale_type = "AutoNorm",
      #               ref_var = "none",
      #               remove_outliers = FALSE
      # )
      
      # input <- list(batch_var = "",
      #               exp_type = "stat",
      #               perc_limit = .99,
      #               filt_type = "none",
      #               miss_type = "rf",
      #               norm_type = "SumNorm",
      #               trans_type = "LogNorm",
      #               scale_type = "AutoNorm",
      #               ref_var = "none",
      #               remove_outliers = FALSE
      # )
      
      # - - check if time series!!! - - 
      
      csv_orig <- fread(global$paths$csv_loc, 
                        data.table = TRUE,
                        header = T)
      
      # Below is your R command history: 
      
      mSet <- InitDataObjects("pktable",
                              "stat",
                              FALSE)
      
      # === Load and re-save CSV ===
      
      csv_orig[,(1:ncol(csv_orig)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)})]
      
      # ============================
      
      csv_orig$sample <- gsub(csv_orig$sample, pattern=" ", replacement="")
      
      as.numi <- as.numeric(colnames(csv_orig)[1:100])
      
      exp.vars <- which(is.na(as.numi))
      
      # === BATCHES ===
    
      batches <- input$batch_var
      
      # === PICK AN INITIAL CONDITION ===|
      
      qc.rows <- which(grepl("QC", csv_orig$sample))
      
      unique.levels <- apply(csv_orig[!qc.rows,..exp.vars, with=F], MARGIN=2, function(col){
        lvls <- levels(as.factor(col))
        # - - - - - -
        length(lvls)
      })
      
      which.default <- unique.levels[which(unique.levels == min(unique.levels[which(unique.levels> 1)]))][1]
      
      condition = names(which.default)
      
      # =================================|
      
      if(is.null(batches)) batches <- ""
      
      batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE
      
      if("batch" %in% batches){ 
        batches = c(batches, "injection")
      }
      
      first_part <- csv_orig[,..exp.vars, with=FALSE]
      first_part[first_part == "" | is.null(first_part)] <- "unknown"
      
      csv <- cbind(first_part[,-c("label")],
                        "label" = first_part[,..condition][[1]],
                        csv_orig[,-..exp.vars,with=FALSE])
      

      # - - - remove outliers? - - -
      
      if(input$remove_outliers){
        sums <- rowSums(csv[,-exp.vars,with=FALSE],na.rm = TRUE)
        names(sums) <- csv$sample
        outliers = c(car::Boxplot(as.data.frame(sums)))
        csv <- csv[!(sample %in% outliers),]
      } 
      
      # - - - remove peaks that are missing in all - - -
      
      csv <- csv[,which(unlist(lapply(csv, function(x)!all(is.na(x))))),with=F]
      
      # - - - low signal samples - - -
      
      complete.perc <- rowMeans(!is.na(csv))
      keep_samps <- csv$sample[which(complete.perc > .2)]
      
      csv <- csv[sample %in% keep_samps,]
      
      covar_table <- first_part[sample %in% keep_samps,]
      
      batchview = if(condition == "batch") TRUE else FALSE
      
      if(any(grepl("QC", csv$sample))){
        samps <- which(!grepl(csv$sample, pattern = "QC"))
        batchnum <- unique(csv[samps, "batch"][[1]])
        keep_samps_post_qc <- covar_table[which(covar_table$batch %in% batchnum),"sample"][[1]]
        covar_table <- covar_table[which(covar_table$batch %in% batchnum),]
        csv <- csv[which(csv$sample %in% keep_samps_post_qc),-"batch"]
      }
      
      colnames(csv)[which( colnames(csv) == "time")] <- "Time"
      
      as.numi <- as.numeric(colnames(csv)[1:100])
      
      exp.vars <- which(is.na(as.numi))
      
      # remove all except sample and time in saved csv
      exp_var_names <- colnames(csv)[exp.vars]
      
      keep_cols <-  c("sample", "label")
      
      remove <- which(!(exp_var_names %in% keep_cols))
      
      print("Removing:")
      print(exp_var_names[remove])
      
      csv_loc_final <- gsub(pattern = "\\.csv", replacement = "_no_out.csv", x = global$paths$csv_loc)
      
      if(file.exists(csv_loc_final)) file.remove(csv_loc_final)
      
      fwrite(csv[,-remove,with=F], file = csv_loc_final)
      
      rownames(covar_table) <- covar_table$sample

      mSet <- Read.TextData(mSet, 
                            filePath = csv_loc_final, 
                            "rowu")  
      
      mSet$dataSet$covars <- covar_table
      
      # - - - sanity check - - -
      
      mSet <- SanityCheckData(mSet)
      
      mSet <- RemoveMissingPercent(mSet, 
                                    percent = input$perc_limit/100)
      
      if(input$miss_type != "none"){
        if(input$miss_type == "pmm"){
          require(mice)
          base <- mSet$dataSet$preproc
          imp <- mice::mice(base, printFlag = TRUE)
        }else if(input$miss_type == "rf"){
          samples <- rownames(mSet$dataSet$preproc)
          
          w.missing <- mSet$dataSet$preproc#[,1:50]
          w.missing <- apply(w.missing, 2, as.numeric)

          #library(doParallel)

          #registerDoParallel(session_cl)
          
          auto.mtry <- floor(sqrt(ncol(mSet$dataSet$preproc)))
          
          mtry <- ifelse(auto.mtry > 100, 100, auto.mtry)
          
          imp <- missForest::missForest(w.missing, 
                                        parallelize = "variables",
                                        verbose = T,
                                        ntree = 10,
                                        mtry = mtry)
          
          mSet$dataSet$procr <- imp$ximp
          rownames(mSet$dataSet$procr) <- rownames(mSet$dataSet$preproc)
          # - - - - - - - - - - - - 
        }else{
          mSet <- ImputeVar(mSet,
                             method = # "knn"
                              input$miss_type
                              )
        }
      }
      # ------------------------
      
      if(input$filt_type != "none"){
        
        mSet <- FilterVariable(mSet,
                                filter = input$filt_type,
                                qcFilter = "F",
                                rsd = 25)
      }
      # ------------------------------------
      
      shiny::setProgress(session=session, value= .2)
      
      if(input$norm_type == "SpecNorm"){
        norm.vec <- mSet$dataSet$covars[match(mSet$dataSet$covars$sample,
                                               rownames(mSet$dataSet$preproc)),][[input$samp_var]]
        norm.vec <- scale(x = norm.vec, center = 1)
        
      }else{
        norm.vec <- rep(1, length(mSet$dataSet$cls))
      }
      
      mSet <- Normalization(mSet,
                             rowNorm = input$norm_type,
                             transNorm = input$trans_type,
                             scaleNorm = input$scale_type,
                             ref = input$ref_var)
      
      shiny::setProgress(session=session, value= .4)
      
      smps <- rownames(mSet$dataSet$norm)
      
      qc_rows <- which(grepl(pattern = "QC", x = smps))
      
      # - - - - - - - - - - - - - - - - - - - - - - -
      
      has.qc <- length(qc_rows) > 0
        
      if(batch_corr){
        
        if("batch" %in% input$batch_var & has.qc){
          
          smpnames = smps[qc_rows]
          
          batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"batch"][[1]]))
          seq.idx = as.numeric(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"injection"][[1]])
          
          # --- QC CORRECTION ---
          
          corr_cols <- pbapply::pblapply(1:ncol(mSet$dataSet$norm), function(i){
            #
            vec = mSet$dataSet$norm[,i]
            #
            corr_vec = BatchCorrMetabolomics::doBC(Xvec = as.numeric(vec), 
                                                   ref.idx = as.numeric(qc_rows), 
                                                   batch.idx = batch.idx,
                                                   seq.idx = seq.idx)
            # ---------
            corr_vec
          })
          
          qc_corr_matrix <- as.data.frame(do.call(cbind, corr_cols))
          
          colnames(qc_corr_matrix) <- colnames(mSet$dataSet$norm)
          rownames(qc_corr_matrix) <- rownames(mSet$dataSet$norm)
          
          mSet$dataSet$norm <- as.data.frame(qc_corr_matrix)
          
        }
        
        if(!batchview){
          mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
          mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
          mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
          mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
        }
        
        left_batch_vars <- grep(input$batch_var, 
                                pattern =  ifelse(has.qc, "batch|injection|sample", "injection|sample"),
                                value = T,
                                invert = T)
        
        if(length(left_batch_vars) > 2){ 
          NULL
        } else if(length(left_batch_vars) == 0){
          NULL
        } else{
          
          smp <- rownames(mSet$dataSet$norm)
          exp_lbl <- mSet$dataSet$cls
          
          csv <- as.data.table(cbind(sample = smp, 
                                     label = mSet$dataSet$cls,
                                     mSet$dataSet$norm))
          
          csv_edata <-t(csv[,!c(1,2)])
          colnames(csv_edata) <- csv$sample
          
          if(length(left_batch_vars) == 1){
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample),left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = c(0),
                                    outcome = as.factor(exp_lbl))     
            batch_normalized= t(sva::ComBat(dat=csv_edata,
                                            batch=csv_pheno$batch1
                                            #mod=mod.pheno,
                                            #par.prior=TRUE
            ))
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }else{
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample), left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample), left_batch_vars[2], with=FALSE][[1]],
                                    outcome = as.factor(exp_lbl))
            batch_normalized = t(limma::removeBatchEffect(x = csv_edata,
                                                          #design = mod.pheno,
                                                          batch = csv_pheno$batch1,
                                                          batch2 = csv_pheno$batch2))
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }
          
          mSet$dataSet$norm <- as.data.frame(batch_normalized)
        }
      } else{
        if(!batchview & has.qc){
          mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
          mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
          mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
          mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
        }
      }
      
      # REORDER COVARS

      mSet$dataSet$covars <- mSet$dataSet$covars[match(rownames(mSet$dataSet$norm), 
                                                        mSet$dataSet$covars$sample),]
      
      mSet$dataSet$cls.name <- condition
      
      # - - set2global - - 
      
      mSet <<- mSet
      
      datamanager$mset_present = TRUE
      
      if(mSet$dataSet$cls.num <= 1){
        interface$mode <- NULL } 
      else if(mSet$dataSet$cls.num == 2){
        interface$mode <- "bivar"}
      else{
        interface$mode <- "multivar"}
      
      # - - - - - - - - - -
      
      shiny::setProgress(session=session, value= .5)
      
      shiny::setProgress(session=session, value= .6)
      varNormPlots <- ggplotNormSummary(mSet)
      output$var1 <- renderPlot(varNormPlots$tl)
      output$var2 <- renderPlot(varNormPlots$bl)
      output$var3 <- renderPlot(varNormPlots$tr)
      output$var4 <- renderPlot(varNormPlots$br)
      
      shiny::setProgress(session=session, value= .7)
      sampNormPlots <-  ggplotSampleNormSummary(mSet)
      output$samp1 <- renderPlot(sampNormPlots$tl)
      output$samp2 <- renderPlot(sampNormPlots$bl)
      output$samp3 <- renderPlot(sampNormPlots$tr)
      output$samp4 <- renderPlot(sampNormPlots$br)
      shiny::setProgress(session=session, value= .8)
      
      mSet$dataSet$adducts <<- selected_adduct_list
      shiny::setProgress(session=session, value= .9)
      
    })
  })
  
  # MAIN EXECUTION OF ANALYSES
  
  observeEvent(input$statistics, {
    
    if(input$nav_general != "analysis") return(NULL)
    
    if(!exists("mSet")) return(NULL)
    
    # get excel table stuff.
    switch(input$statistics,
           pca = {
             if(!"pca" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- PCA.Anal(mSet)
               })
             }
             output$pca_scree <- renderPlot({
               df <- data.table(
                 pc = 1:length(names(mSet$analSet$pca$variance)),
                 var = mSet$analSet$pca$variance)
               p <- ggplot2::ggplot(data=df[1:20,]) + ggplot2::geom_line(mapping = aes(x=pc, y=var, colour=var), cex=3) + 
                 global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]](base_size = 10) +
                 ggplot2::scale_colour_gradientn(colours = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](256)) +
                 theme(axis.text=element_text(size=10),
                       axis.title=element_text(size=19,face="bold"),
                       #legend.title=element_text(size=15),
                       #legend.text=element_text(size=12),
                       legend.position="none")
               # - - - - - 
               p
               #ggplotly(p)
             })
             output$plot_pca <- plotly::renderPlotly({
               plotPCA.3d(mSet, global$vectors$mycols,
                          pcx = input$pca_x,
                          pcy = input$pca_y,
                          pcz = input$pca_z,
                          shape.fac = input$second_var)
             })
             # === LEGEND ===
             output$pca_legend <- plotly::renderPlotly({
               frame <- data.table(x = c(1), 
                                   y = fac.lvls)
               p <- ggplot(data=frame,
                           aes(x, 
                               y, 
                               color=factor(y),
                               fill=factor(y))) + 
                 geom_point(shape = 21, size = 5, stroke = 5) +
                 scale_colour_manual(values=chosen.colors) +
                 theme_void() + 
                 theme(legend.position="none")
               # --- return ---
               ggplotly(p, tooltip = NULL) %>% config(displayModeBar = F)
             })
             # ==============
             output$pca_tab <-DT::renderDataTable({
               pca.table <- as.data.table(round(mSet$analSet$pca$variance * 100.00,
                                                digits = 2),
                                          keep.rownames = T)
               colnames(pca.table) <- c("Principal Component", "% variance")
               
               DT::datatable(pca.table, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
             
             output$pca_load_tab <-DT::renderDataTable({
               pca.loadings <- mSet$analSet$pca$rotation[,c(input$pca_x,
                                                            input$pca_y,
                                                            input$pca_z)]
               #colnames(pca.loadings)[1] <- "m/z"
               DT::datatable(pca.loadings, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 10))
             })
             },
           meba = {
             if("MB" %not in% names(mSet$analSet)){
               mSet <<- performMB(mSet, 10)
             }
             output$meba_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$MB$stats, 
                             selection = 'single',
                             colnames = c("Compound", "Hotelling/T2 score"),
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
           },
           asca = {
             if("asca" %not in% names(mSet$analSet)){
               mSet <<- Perform.ASCA(mSet, 1, 1, 2, 2)
               mSet <<- CalculateImpVarCutoff(mSet, 0.05, 0.9)
             }
             output$asca_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$asca$sig.list$Model.ab, 
                             selection = 'single',
                             colnames = c("Compound", "Leverage", "SPE"),
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
             
           },
           heatmap = {
             
             withProgress({
               
               output$heatmap <- plotly::renderPlotly({
                 
                 if(interface$mode == "bivar"){
                   if(input$heatmode){
                     tbl = as.data.frame(mSet$analSet$tt$sig.mat)
                     used.values <- "p.value"
                     decreasing = F
                   }else{
                     tbl = as.data.frame(mSet$analSet$fc$sig.mat)
                     tbl$abs_log2 <- abs(tbl$`log2(FC)`)
                     used.values <- "abs_log2"
                     decreasing = T
                   }
                 }else if(interface$mode == "multivar"){
                   tbl = as.data.frame(mSet$analSet$aov$sig.mat)
                   used.values = "p.value"
                   decreasing = F
                 }else{
                   if(input$heatmode){
                     tbl = as.data.frame(mSet$analSet$asca$sig.list$Model.ab)
                     used.values = "Leverage"
                   }else{
                     tbl = as.data.frame(mSet$analSet$MB$stats)
                     used.values = "Hotelling-T2"
                   }
                   decreasing = T
                 }
                 
                 topn = if(length(tbl[[used.values]]) < input$heatmap_topn) length(tbl[[used.values]]) else input$heatmap_topn
                 
                 mzorder <- order(tbl[[used.values]], decreasing = decreasing)
                 mzsel <- rownames(tbl)[mzorder][1:topn]
                 
                 x <- mSet$dataSet$norm[,mzsel]
                 final_matrix <<- t(x)
                 
                 sample_order <- match(colnames(final_matrix), rownames(mSet$dataSet$norm))
                 
                 if(timebutton$status == "on"){
                   if(input$timecourse_trigger){
                     translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$exp.fac[sample_order], Time=mSet$dataSet$time.fac[sample_order])
                     hmap.lvls <- c(levels(mSet$dataSet$exp.fac), levels(mSet$dataSet$time.fac))

                     split.translator <- split(translator, by = c("Time"))
                     split.translator.ordered <- lapply(split.translator, function(tbl) tbl[order(tbl$Group)])
                     translator <- rbindlist(split.translator.ordered)
                     
                     final_matrix <<- final_matrix[,match(translator$Sample, colnames(final_matrix))]
                     
                     my_order=F
                     
                   }else{
                     translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$cls[sample_order])
                     hmap.lvls <- levels(mSet$dataSet$cls)
                     my_order = T
                   }
                 }else{
                   translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$cls[sample_order])
                   hmap.lvls <- levels(mSet$dataSet$cls)
                   my_order = T
                 } 

                 color.mapper <- {
                  classes <- hmap.lvls
                  cols <- sapply(1:length(classes), function(i) global$vectors$mycols[i])
                  names(cols) <- classes
                  # - - -
                  cols
                 }
                
                 hmap <- heatmaply::heatmaply(final_matrix,
                                              Colv=my_order, 
                                              Rowv=T,
                                              branches_lwd = 0.3,
                                              margins = c(60, 0, NA, 50),
                                              colors = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](256),
                                              col_side_colors = translator[,!1],
                                              col_side_palette = color.mapper,
                                              subplot_widths = c(.9,.1),
                                              subplot_heights = if(my_order) c(.1, .05, .85) else c(.05,.95),
                                              column_text_angle = 90,
                                              xlab = "Sample",
                                              ylab = "m/z",
                                              showticklabels = c(T,F)
                                              #label_names = c("m/z", "sample", "intensity") #breaks side colours
                 )
                 
                 hmap_mzs <<- hmap$x$layout$yaxis3$ticktext
                 
                 # - - 
                 
                 hmap
                 
               })               
             })
           },
           tt = {
             if(!"tt" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- Ttests.Anal(mSet,
                                      nonpar = FALSE, 
                                      threshp = 0.05, 
                                      paired = FALSE,
                                      equal.var = TRUE
                                      #  multicorr = "BH"
                 )
               })
             }
             res <<- mSet$analSet$tt$sig.mat
             if(is.null(res)) res <<- data.table("No significant hits found")
             output$tt_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
             output$tt_overview_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotTT(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           },
           fc = {
             if(!"fc" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- FC.Anal.unpaired(mSet,
                                           2.0, 
                                           1)                 
               })
             }
             res <<- mSet$analSet$fc$sig.mat
             if(is.null(res)) res <<- data.table("No significant hits found")
             output$fc_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
             output$fc_overview_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotFC(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           },
           aov = {
             if(!"aov" %in% names(mSet$analSet)){
               withProgress({
                 ANOVA.Anal(mSet, thresh=0.05,nonpar = F)
               # TODO: Fix ANOVA2
               # mSet <<- ifelse(timebutton_trigger,
               #                 {mSet$dataSet$facA <- mSet$dataSet$exp.fac;
               #                  mSet$dataSet$facB <- mSet$dataSet$time.fac;
               #                  ANOVA2.Anal(mSet, 0.05, "fdr", "time", 3, 1)
               #                   },
               #                 ANOVA.Anal(mSet, thresh=0.05,nonpar = F))
               })
             }
    
             output$aov_tab <-DT::renderDataTable({
              
               mSet$analSet$aov$sig.mat
               
               # -------------
               DT::datatable(if(is.null(mSet$analSet$aov$sig.mat)) data.table("No significant hits found") else mSet$analSet$aov$sig.mat, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
           },
           volc = {
             if(!"volc" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- Volcano.Anal(mSet,FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
               })
             }
             output$volc_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$volc$sig.mat, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
             output$volc_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotVolc(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           }
    )
  })
  
  observeEvent(input$do_plsda, {
    
    library(e1071)
    
    switch(input$plsda_type,
           normal={
               require(caret)

               mSet <<- PLSR.Anal(mSet)
               mSet <<- PLSDA.CV(mSet, methodName=if(nrow(mSet$dataSet$norm) < 50) "L" else "T",compNum = 5)
               mSet <<- PLSDA.Permut(mSet,num = 100, type = "accu")
               })
  })
  
  observe({
    
    if(!exists("mSet")) return()
    
    # - - - - - - - - - - -
    
    if("plsda" %in% names(mSet$analSet)){
      output$plsda_cv_plot <- renderPlot({
        ggPlotClass(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
      })
      output$plsda_perm_plot <- renderPlot({
        ggPlotPerm(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
      })
      output$plot_plsda_3d <- plotly::renderPlotly({
        plotPCA.3d(mSet, cols = global$vectors$mycols, input$second_var, input$plsda_x, input$plsda_y, input$plsda_z, mode = "plsda")
      })
      output$plsda_tab <- DT::renderDataTable({
        # - - - -
        plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                           / mSet$analSet$plsr$Xtotvar 
                                           * 100.0,
                                           digits = 2),
                                     keep.rownames = T)
        colnames(plsda.table) <- c("Principal Component", "% variance")
        plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
        # -------------
        DT::datatable(plsda.table, 
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })
      output$plsda_load_tab <-DT::renderDataTable({
        plsda.loadings <- mSet$analSet$plsda$vip.mat
        colnames(plsda.loadings) <- paste0("PC", c(1:ncol(plsda.loadings)))
        # -------------
        DT::datatable(plsda.loadings[, c(input$plsda_x, input$plsda_y, input$plsda_z)],
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })
      output$plot_plsda <- plotly::renderPlotly({
        plotPCA.3d(mSet, global$vectors$mycols,
                   pcx = input$plsda_x,
                   pcy = input$plsda_y,
                   pcz = input$plsda_z, mode = "plsda",
                   shape.fac = input$second_var)
      })
    }
  })
  
  observe({
    
    if(!exists("mSet")) return()
    
    # - - - - - - - - - - -
    
    if("plsda" %in% names(mSet$analSet)){
      output$pla_cv_plot <- renderPlot({
        ggPlotClass(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
      })
      output$plsda_perm_plot <- renderPlot({
        ggPlotPerm(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
      })
      output$plot_plsda_3d <- plotly::renderPlotly({
        plotPCA.3d(mSet, cols = global$vectors$mycols, input$second_var, input$plsda_x, input$plsda_y, input$plsda_z, mode = "plsda")
      })
      output$plsda_tab <- DT::renderDataTable({
        # - - - -
        plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                           / mSet$analSet$plsr$Xtotvar 
                                           * 100.0,
                                           digits = 2),
                                     keep.rownames = T)
        colnames(plsda.table) <- c("Principal Component", "% variance")
        plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
        # -------------
        DT::datatable(plsda.table, 
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })
      output$plsda_load_tab <-DT::renderDataTable({
        plsda.loadings <- mSet$analSet$plsda$vip.mat
        colnames(plsda.loadings) <- paste0("PC", c(1:ncol(plsda.loadings)))
        # -------------
        DT::datatable(plsda.loadings[, c(input$plsda_x, input$plsda_y, input$plsda_z)],
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })
      output$plot_plsda <- plotly::renderPlotly({
        plotPCA.3d(mSet, global$vectors$mycols,
                   pcx = input$plsda_x,
                   pcy = input$plsda_y,
                   pcz = input$plsda_z, mode = "plsda",
                   shape.fac = input$second_var)
      })
    }
  })
  
  observeEvent(input$pca_2d3d,{

    if(!exists("mSet")) return()
    
    if(!("pca" %in% names(mSet$analSet))) return()
    
    mode <- if("timecourse_trigger" %in% names(input)){
      if(input$timecourse_trigger){
        "ipca"
      }else{
        "pca"
      } 
    }else{
      "pca"
      }
    
    # - - - - -
    if(input$pca_2d3d){
      # 2d
      output$plot_pca <- plotly::renderPlotly({
        plotPCA.2d(mSet, global$vectors$mycols,
                   pcx = input$pca_x,
                   pcy = input$pca_y, mode = mode,
                   shape.fac = input$second_var)
      })
    }else{
      # 3d
      output$plot_pca <- plotly::renderPlotly({
        plotPCA.3d(mSet, global$vectors$mycols,
                   pcx = input$pca_x,
                   pcy = input$pca_y,
                   pcz = input$pca_z, mode = mode,
                   shape.fac = input$second_var)
      })
    }
  })
  
  observeEvent(input$plsda_2d3d,{
    
    if(!exists("mSet")) return()
    
    if(!("plsda" %in% names(mSet$analSet) | "plsr" %in% names(mSet$analSet))) return()
    
    if(input$plsda_2d3d){
      # 2d
      output$plot_plsda <- plotly::renderPlotly({
        plotPCA.2d(mSet, global$vectors$mycols,
                   pcx = input$plsda_x,
                   pcy = input$plsda_y, mode = "plsda",
                   shape.fac = input$second_var)
      })
    }else{
      # 3d
      output$plot_plsda <- plotly::renderPlotly({
        plotPCA.3d(mSet, global$vectors$mycols,
                   pcx = input$plsda_x,
                   pcy = input$plsda_y,
                   pcz = input$plsda_z, mode = "plsda",
                   shape.fac = input$second_var)
      })
    }
  })
  
  observeEvent(input$do_ml, {
    
    # do_ml, ml_attempts, ml_train_perc
    
    
    # -- NOTE : CHECK IF AUC INCREASES W/ INCREASING TRAINING SET SIZE --- !!!!!!!!!!!!!!!!
    # 'saturation point'
    #switchButton(inputId = "ml_saturation", label = "SatMode", value = FALSE, col = "BW", type = "OO"),
    
    withProgress({
      
      setProgress(value = 0)
      # prepare matric
      
      curr <- as.data.table(mSet$dataSet$preproc)
      
      curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
      
      config <- mSet$dataSet$covars[match(mSet$dataSet$covars$sample,rownames(mSet$dataSet$preproc)),]
      config <- config[!is.na(config$sample),]
      config$label <- mSet$dataSet$cls
      config <- cbind(config, label=mSet$dataSet$cls)
      
      # bye bye NAs
      
      config <- config[,apply(!is.na(config), 2, any), with=FALSE]
      
      
      # - - - - - -
      keep_curr <- match(mSet$dataSet$covars$sample,rownames(mSet$dataSet$preproc))
      
      curr <- cbind(config, curr[keep_curr])
      
      curr <- curr[which(!grepl(config$sample,
                                  pattern = "QC"))]        
      

      # remove cols with all NA
      
      curr <- curr[,colSums(is.na(curr))<nrow(curr),with=FALSE]
      
      configCols <- which(!(colnames(curr) %in% colnames(mSet$dataSet$norm)))
      mzCols <- which(colnames(curr) %in% colnames(mSet$dataSet$norm))
      
      curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
      curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]
      
      goes = as.numeric(input$ml_attempts)
      
      repeats <- pbapply::pblapply(1:goes, function(i){
        
        # train / test
        #ml_train_regex=ml_test_regex=""
        #ml_train_perc=.8
        ml_train_regex <<- input$ml_train_regex
        ml_test_regex <<- input$ml_test_regex
        
        ml_train_perc <- input$ml_train_perc/100
        
        if(ml_train_regex == "" & ml_test_regex == ""){ # BOTH ARE NOT DEFINED
          test_idx = caret::createDataPartition(y = curr$label, p = ml_train_perc, list = FALSE)
          train_idx = setdiff(1:nrow(curr), test_idx)
          inTrain = train_idx
          inTest = test_idx
        }else if(ml_train_regex != ""){ #ONLY TRAIN IS DEFINED
          train_idx = grep(config$sample, pattern = ml_train_regex)
          test_idx = setdiff(1:nrow(curr), train_idx)
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc)
          inTrain <- train_idx[reTrain$Resample1]
          inTest = test_idx
        }else{ # ONLY TEST IS DEFINED
          test_idx = grep(config$sample, pattern = ml_test_regex)
          train_idx = setdiff(1:nrow(curr), test_idx)
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc)
          inTrain <- train_idx[reTrain$Resample1]
          
          inTest <- test_idx
          #reTest <- caret::createDataPartition(y = config[test_idx, label], p = ml_train_perc)
          #inTest <- test_idx[reTest$Resample1]
        }
        
        # - - - re-split - - -
        #reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc)
        #inTrain <- train_idx[reTrain$Resample1]
        #reTest <- caret::createDataPartition(y = config[test_idx, label], p = ml_train_perc)
        #inTest <- test_idx[reTest$Resample1]
        
        # - - divide - -
        
        predictor = "label"
        
        trainY <- curr[inTrain, 
                       ..predictor][[1]]
        testY <- curr[inTest,
                      ..predictor][[1]]
        
        group.cols <- grep(colnames(curr), pattern = "^group",value = T)
        
        remove.cols <- c("sample", "label", group.cols, "Stool_condition") #TODO: make group column removed at start
        remove.idx <- which(colnames(curr) %in% remove.cols)
        training <- curr[inTrain,-remove.idx, with=FALSE]
        testing <- curr[inTest,-remove.idx, with=FALSE]
        
        predIdx <- which(colnames(curr) %in% colnames(config))
        
        #trainX <- apply(training, 2, as.numeric)
        #testX <- apply(testing, 2, as.numeric)
        
        training <- data.matrix(gdata::drop.levels(training))
        testing <- data.matrix(gdata::drop.levels(testing))
        
        setProgress(value = i/goes)
        
        switch(input$ml_method,
               rf = {
                 model = randomForest::randomForest(x = training, 
                                                    y = trainY, 
                                                    ntree = 500,
                                                    importance=TRUE)
                 
                 prediction <- stats::predict(model, 
                                              testing, 
                                              "prob")[,2]
                 
                 importance = as.data.table(model$importance, keep.rownames = T)
                 rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
                 rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
                 rf_tab <- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
                 list(type="rf",
                      feats = as.data.table(rf_tab, keep.rownames = T), 
                      model = model,
                      prediction = prediction,
                      labels = testY)
               }, 
               ls = {

                 nfold = switch(input$ml_folds, 
                                "5" = 5,
                                "10" = 10,
                                "20" = 20,
                                "50" = 50,
                                "LOOCV" = length(trainY))
                 family = "binomial"
                 
                 cv1 <- glmnet::cv.glmnet(training, trainY, family = family, type.measure = "auc", alpha = 1, keep = TRUE, nfolds=nfold)
                 cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[["lambda.min"]]], lambda = cv1[["lambda.min"]], alpha = 1)
                 
                 model <- glmnet::glmnet(as.matrix(training), trainY, family = family, lambda = cv2$lambda, alpha = cv2$alpha)
                 
                 prediction <- stats::predict(model,
                                              type = "response", 
                                              newx = testing, 
                                              s = "lambda.min")#[,2] # add if necessary
                 list(type = "ls",
                      model = model,
                      prediction = prediction, 
                      labels = testY)
               }, 
               gls = {
                 NULL
               })
      })
      
      if(!"ml" %in% names(mSet$analSet)){
        mSet$analSet$ml <<- list(ls=list(), rf=list())
      }
      
      #input <- list(ml_attempts = 50, ml_method = "ls", ml_name = "test")
      
      xvals <<- list(type = {unique(lapply(repeats, function(x) x$type))},
                     models = {lapply(repeats, function(x) x$model)},
                     predictions = {lapply(repeats, function(x) x$prediction)},
                     labels = {lapply(repeats, function(x) x$labels)})
      
      mSet$analSet$ml[[input$ml_method]][[input$ml_name]] <<- list("roc" = xvals,
                                                                  "bar" = repeats)
      
      output$ml_roc <- plotly::renderPlotly({plotly::ggplotly(ggPlotROC(xvals, input$ml_attempts, global$functions$color.functions[[getOptions("user_options.txt")$gspec]]))})
      output$ml_bar <- plotly::renderPlotly({plotly::ggplotly(ggPlotBar(repeats, input$ml_attempts, global$functions$color.functions[[getOptions("user_options.txt")$gspec]], input$ml_top_x, input$ml_name))})
      
    })
  })
  
  observe({
    # ---------------------------------
    db_search_list <- lapply(global$vectors$db_list, 
                             FUN = function(db){
                               # -----------------------
                               if(!input[[paste0("search_", db)]]){
                                 c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db")))
                               }
                               else{NA}
                             }
    )
    global$vectors$db_search_list <<- db_search_list[!is.na(db_search_list)]
  })  
  
  observe({
    # ---------------------------------
    add_search_list <- lapply(global$vectors$db_list, 
                              FUN = function(db){
                                if(is.null(input[[paste0("add_", db)]])){
                                  return(NA)
                                }else{
                                  # -----------------------
                                  if(!input[[paste0("add_", db)]]){
                                    c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db")))
                                  }
                                  else{NA}
                                }
                              }        
    )
    global$vectors$add_search_list <<- add_search_list[!is.na(add_search_list)]
  })  
  
  observe({
    # ---------------------------------
    enrich_db_list <- lapply(global$vectors$db_list, 
                             FUN = function(db){
                               if(is.null(input[[paste0("enrich_", db)]])){
                                 return(NA)
                               }else{
                                 # -----------------------
                                 if(!input[[paste0("enrich_", db)]]){
                                   c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db")))
                                 }
                                 else{NA}
                               }
                             }
    )
    global$vectors$enrich_db_list <<- enrich_db_list[!is.na(enrich_db_list)]
  })  
  
  observe({
    # --------------
    wanted.adducts.pos <- global$vectors$pos_adducts[input$pos_add_tab_rows_selected, "Name"]
    wanted.adducts.neg <- global$vectors$neg_adducts[input$neg_add_tab_rows_selected, "Name"]
    # ---------
    selected_adduct_list <<- rbind(wanted.adducts.neg, 
                                   wanted.adducts.pos)$Name
  })
  
  res.update.tables <<- c("tt", 
                          "fc", 
                          "aov",
                          "rf",
                          "asca", 
                          "meba",
                          "pca_load",
                          "plsda_load",
                          "enrich_pw",
                          "ml",
                          paste0("alpha_", seq(0,1,0.2), "_lasnet"),
                          paste0("alpha_stab_", seq(0,1,0.2), "_lasnet"))
  
  observeEvent(input$enriched_rows_selected, {
    curr_row = input$enriched_rows_selected
    # LOAD MEMBERS IN GROUP
    curr_pw <- enrich_tab[curr_row,]
    pw_memb_here <- gsaRes$gsc[[curr_pw$Name]]
    pw_memb_tab <<- as.data.frame(pw_memb_here)
    rownames(pw_memb_tab) <- pw_memb_here
    pw_memb_tab[,1] <- c("Yes")
    colnames(pw_memb_tab) <- "Present"
    # get missing members
    # pw_memb_all <- gset_proc$gsc[[curr_pw$Name]]
    # pw_memb_all_tab <<- as.data.frame(pw_memb_all)
    # rownames(pw_memb_all_tab) <- pw_memb_all
    # pw_memb_all_tab[,1] <- c("No")
    # colnames(pw_memb_all_tab) <- "Present"
    # join
    enrich_overview_tab <<- pw_memb_tab #rbind(pw_members_tab, pw_memb_all_tab)
    #
    output$enrich_pw_tab <-DT::renderDataTable({
      DT::datatable(enrich_overview_tab,
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), 
                                   pageLength = 5))
    })  
  })
  
  lapply(unique(res.update.tables), FUN=function(table){
    observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
      curr_row = input[[paste0(table, "_tab_rows_selected")]]
      # do nothing if not clicked yet, or the clicked cell is not in the 1st column
      if (is.null(curr_row)) return()
      
      if(grepl(pattern = "lasnet",x = table)){
        alpha <- gsub("alpha|stab|_|lasnet", "", table)
        if(!grepl(pattern = "stab",x = table)){
          table <- "lasnet_a"
        }else{ 
          table <- "lasnet_b"
        }
      }
      curr_cpd <<- data.table::as.data.table(switch(table,
                                                    tt = mSet$analSet$tt$sig.mat,
                                                    fc = mSet$analSet$fc$sig.mat,
                                                    pca_load = mSet$analSet$pca$rotation,
                                                    plsda_load = mSet$analSet$plsda$vip.mat,
                                                    ml = ml_tab,
                                                    asca = mSet$analSet$asca$sig.list$Model.ab,
                                                    aov = mSet$analSet$aov$sig.mat,
                                                    rf = vip.score,
                                                    enrich_pw = enrich_overview_tab,
                                                    meba = mSet$analSet$MB$stats,
                                                    plsda_vip = plsda_tab)
                                             , keep.rownames = T)[curr_row, rn]
      
      if(grepl(pattern = "lasnet",x = table)){
        outplot_name = "lasnet_specific_plot"
      }else{
        outplot_name = paste0(table, "_specific_plot")
      }
      
      output$curr_cpd <- renderText(curr_cpd)

      output$curr_plot <- plotly::renderPlotly({
        # --- ggplot ---
        ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
      })
      
      output[[outplot_name]] <- plotly::renderPlotly({
        # --- ggplot ---
        if(table == 'meba'){
          ggplotMeba(curr_cpd, draw.average = T, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
        }else if(table == 'asca'){
          ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]], mode = "ts")
        }else{
          ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
        }
      })
    })
  })
  

  observeEvent(plotly::event_data("plotly_click"),{ # WHAT HAPPENS ON CLICK
    d <- plotly::event_data("plotly_click")
    req(d)
    if(input$statistics %in% c("tt", "fc", "rf", "aov", "volc", "lasnet")){
        if('key' %not in% colnames(d)) return(NULL)
        mzs <- switch(input$statistics, 
                      tt= names(mSet$analSet$tt$p.value),
                      fc = names(mSet$analSet$fc$fc.log),
                      aov = names(mSet$analSet$aov$p.value),
                      volc = rownames(mSet$analSet$volcano$sig.mat)
        )
        if(d$key %not in% mzs) return(NULL)
        curr_cpd <<- d$key
        # - return -
        output[[paste0(input$statistics, "_specific_plot")]] <- plotly::renderPlotly({
          # --- ggplot ---
          ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
        })
      }else if(input$statistics == "pca"){
        if(!"z" %in% names(d)){
          which_group = 1
          which_group <- d$curveNumber + 1
          traceLoc <- length(unique(mSet$dataSet$cls)) + 1
          scatter <- pca_plot$x$attrs[[traceLoc]]
          idx <- scatter$color == which_group
          if(pca_plot$x$data[[which_group]]$visible == "legendonly"){
            pca_plot$x$data[[which_group]]$visible = TRUE
            scatter$visible[idx] <- T
          }else{ # hide
            pca_plot$x$data[[which_group]]$visible = "legendonly"
            scatter$visible[idx] <- F
          }
          pca_plot$x$attrs[[traceLoc]] <- scatter
          pca_plot <<- pca_plot
          output$plot_pca <- plotly::renderPlotly({pca_plot})
        }}else if(input$statistics == "ml"){
          switch(input$ml_results, roc = {
            attempt = d$curveNumber - 1
            if(attempt > 1){
              ml_type <- xvals$type[[1]]
              model <- xvals$models[[attempt]]
              output$ml_tab <- switch(ml_type,
                                      rf = {
                                        importance = as.data.table(model$importance, keep.rownames = T)
                                        rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
                                        rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
                                        # - - - return - - -
                                        ml_tab <<- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
                                        DT::renderDataTable({
                                          DT::datatable(rf_tab,
                                                        selection = 'single',
                                                        autoHideNavigation = T,
                                                        options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                                        })
                                      },
                                      ls = {
                                        tab = model$beta
                                        keep = which(tab[,1] > 0)
                                        tab_new = data.frame("beta" = tab[keep,1],
                                                             "absbeta" = abs(tab[keep,1]),
                                                             row.names = rownames(tab)[keep])
                                        colnames(tab_new) <- c("beta", "abs_beta")
                                        ml_tab <<- tab_new[order(tab_new[,1],decreasing = TRUE),]
                                        DT::renderDataTable({
                                          DT::datatable(ml_tab,
                                                        selection = 'single',
                                                        autoHideNavigation = T,
                                                        options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                                        })
                                      })
              
            }
          }, bar = {
            #
            #d = list(x=1)
            curr_cpd <<- as.character(global$tables$ml_bar_tab[d$x,"mz"][[1]])
            
            output$ml_specific_plot <- plotly::renderPlotly({
              # --- ggplot ---
              ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
            })
          })}else if(grepl(pattern = "heatmap", x = input$statistics)){
            if(!exists("hmap_mzs")) return(NULL)
            if(d$y > length(hmap_mzs)) return(NULL)
            curr_cpd <<- hmap_mzs[d$y]
          }
  
    output$curr_plot <- plotly::renderPlotly({
      # --- ggplot ---
      ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
    })
    
    # ----------------------------
    output$curr_cpd <- renderText(curr_cpd)
  })
  
  # --- find matches ---
  
  output$find_mol_icon <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('www/search.png'))
    # Return a list containing the filename and alt text
    list(src = filename,
         width=70,
         height=70)
  }, deleteFile = FALSE)
  
  observeEvent(input$search_cpd, {
    req(global$vectors$db_search_list)
    # ----------------
    if(length(global$vectors$db_search_list) > 0){
      global$tables$last_matches <<- unique(multimatch(curr_cpd, global$vectors$db_search_list))
      output$match_tab <- DT::renderDataTable({
        DT::datatable(global$tables$last_matches[,-c("description","structure", "baseformula", "dppm")],
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })  
      }
  })
  
  observeEvent(input$match_tab_rows_selected,{
    curr_row = input$match_tab_rows_selected
    curr_row <<- input$match_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_def <<- global$tables$last_matches[curr_row,'description']
    output$curr_definition <- renderText(curr_def$description)
    curr_struct <<- global$tables$last_matches[curr_row,'structure'][[1]]
    output$curr_struct <- renderPlot({plot.mol(curr_struct,style = "cow")})
    curr_formula <<- global$tables$last_matches[curr_row,'baseformula'][[1]]
    output$curr_formula <- renderText({curr_formula})
  })
  
  observeEvent(input$browse_db,{
    #req(db_search_list)
    # -------------------
    cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
      browse_db(match.table)
    })
    # ------------------
    browse_table <<- unique(as.data.table(rbindlist(cpd_list)))
    
    output$browse_tab <-DT::renderDataTable({
      DT::datatable(browse_table[,-c("Description", "Charge")],
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
    })
  })
  
  observeEvent(input$revsearch_cpd, {
    #req(db_search_list)
    req(input$browse_tab_rows_selected)
    # -------------------
    search_cmd <- browse_table[curr_row,c('Formula', 'Charge')]
    # -------------------
    cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
      get_mzs(search_cmd$Formula, search_cmd$Charge, match.table)})
    # ------------------
    hits_table <<- unique(as.data.table(rbindlist(cpd_list)))
    output$hits_tab <-DT::renderDataTable({
      DT::datatable(hits_table,
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
    })
  })
  
  observeEvent(input$browse_tab_rows_selected,{
    curr_row = input$browse_tab_rows_selected
    curr_cpd = browse_table[curr_row, Formula]
    curr_row <<- input$browse_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_def <<- browse_table[curr_row,'Description']
    output$browse_definition <- renderText(curr_def$Description)
    # --- search ---
    #TODO: this should be a function and not re-written
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
  })
  
  observeEvent(input$hits_tab_rows_selected,{
    curr_row = input$hits_tab_rows_selected
    curr_row <<- input$hits_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_cpd <<- hits_table[curr_row, mzmed.pgrp]
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
  })
  
  observeEvent(input$go_enrich,{
    #TODO: fix, currently broken (build in m/z > cpd search?)
    enrich_db_list <- paste0("./backend/db/", c("kegg", 
                                                "wikipathways", 
                                                "smpdb",
                                                "metacyc"), ".full.db")
    withProgress({
      all_pathways <- lapply(enrich_db_list, FUN=function(db){
        conn <- RSQLite::dbConnect(RSQLite::SQLite(), db) # change this to proper var later
        dbname <- gsub(basename(db), pattern = "\\.full\\.db", replacement = "")
        #RSQLite::dbGetQuery(conn, "SELECT distinct * FROM pathways limit 10;")
        gset <- RSQLite::dbGetQuery(conn,"SELECT DISTINCT c.baseformula AS cpd,
                                    p.name as name
                                    FROM pathways p
                                    JOIN base c
                                    ON c.pathway = p.identifier 
                                    ")
        RSQLite::dbDisconnect(conn)
        # --- returny ---
        gset$name <- paste(gset$name, " (", dbname, ")", sep="")
        gset
      })
      shiny::setProgress(session=session, value= 0.33)
      # --- only anova top 100 for now ---
      used.analysis <- input$enrich_stats
      if(used.analysis == "rf"){
        ranking <- vip.score[accuracyDrop > 0,]
        ranking <- ranking[order(-accuracyDrop),]
        vec <- ranking$accuracyDrop
        names(vec) <- ranking$rn
        sigvals <- switch(input$enrich_vals,
                          sig=vec[1:500],
                          t50=vec[1:50],
                          t100=vec[1:100],
                          t200=vec[1:200],
                          t500=vec[1:500]
        )
      }else{
        stat.tab <- switch(input$enrich_stats,
                           tt="p.value",
                           aov="p.value",
                           fc="fc.all")
        sigvals <- switch(input$enrich_vals,
                          sig=mSet$analSet[[used.analysis]][[stat.tab]][mSet$analSet[[used.analysis]]$inx.imp],
                          t50=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:50],
                          t100=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:100],
                          t200=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:200],
                          t500=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:500]
        )        
      }
      
      # -----------------------
      gset <- rbindlist(all_pathways)
      gset_proc <<- piano::loadGSC(gset)
      shiny::setProgress(session=session, value= 0.66)
      
      # - - - - - - 
      
      sigvals <- rownames(lasnet_tables[[1]])
      matches <- rbindlist(lapply(sigvals, function(mz){
        subtables <- lapply(enrich_db_list, function(db){
          matches = get_matches(mz, db, F, "mz")
        }) 
        tab <- rbindlist(subtables)
        tab$mz = c(mz)
        tab
      }))
      
      lasnet_vals <- lasnet_tables[[1]]
      lasnet_vals$mz <- rownames(lasnet_vals)
      matches <- merge(matches, lasnet_vals, by = "mz")
      sigvals <- unique(matches[, c("Baseformula", "beta")])
      
      # - - - - - -
      
      gsaRes <<- piano::runGSA(sigvals, 
                               gsc = gset_proc)
      enrich_tab <<- piano::GSAsummaryTable(gsaRes)[,1:3]
      # --- render ---
      output$enriched <- DT::renderDataTable({
        # -------------
        DT::datatable(enrich_tab, 
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
        
      })
      shiny::setProgress(session=session, value= 1)
    })
  })
  
  
  # === META ANALYSIS ===
  
  # --- combine tab ---
  
  observeEvent(input$venn_members, {
    
    if("ls" %in% input$venn_members | "rf" %in% input$venn_members | "plsda" %in% input$venn_members){
      uiElements <- lapply(intersect(c("rf", "ls", "plsda"), input$venn_members), function(name){
        options <- switch(name,
                          tt = NULL,
                          fc = NULL,
                          ls = names(mSet$analSet$ml$ls),
                          rf = names(mSet$analSet$ml$rf),
                          volc = NULL,
                          plsda = c("PC1", "PC2", "PC3"))
        
        shorthand = switch(name,
                           ls = "lasso",
                           rf = "random forest",
                           plsda = "PLSDA")
        selectInput(inputId = paste0(name,"_choice"), 
                    label = gsubfn::fn$paste("Which $shorthand group(s):"), 
                    multiple=TRUE,
                    choices = options)
      })
      
      output$venn_ml_ui <- renderUI({uiElements})
    }else{
      output$venn_ml_ui <- renderUI({br()})
    }
  })
  
  observeEvent(input$build_venn, {
    
    # - - - cats - - -
    
    top = input$venn_tophits
    
    categories <- input$venn_members
    
    #top = 200
    #categories = c("tt", "fc")
    
    tables <- lapply(categories, function(name){
      
      tbls <- switch(name,
                     tt = list(rownames(mSet$analSet$tt$sig.mat[order(mSet$analSet$tt$sig.mat[,2], decreasing = F),])),
                     fc = list(rownames(mSet$analSet$fc$sig.mat[order(abs(mSet$analSet$fc$sig.mat[,2]), decreasing = F),])),
                     ls = {
                       tbls_ls <- lapply(input$ls_choice, function(name){
                         mSet$analSet$ml$ls[[name]][order(mSet$analSet$ml$ls[[name]]$count, decreasing = T),]$mz})
                       names(tbls_ls) <- input$ls_choice
                       # - - - 
                       tbls_ls
                     },
                     rf = {
                       tbls_rf <- lapply(input$rf_choice, function(name){
                         mSet$analSet$ml$rf[[name]][order(mSet$analSet$ml$rf[[name]]$mda, decreasing = T),]$mz})
                       names(tbls_rf) <- input$rf_choice
                       # - - -
                       tbls_rf
                     },
                     plsda = {
                       tbls_plsda <- lapply(input$plsda_choice, function(name){
                         compounds_pc <- as.data.table(mSet$analSet$plsda$vip.mat,keep.rownames = T)
                         colnames(compounds_pc) <- c("rn", paste0("PC", 1:3))
                         ordered_pc <- setorderv(compounds_pc, name, -1)
                         ordered_pc[, c("rn")][[1]]               
                       })
                       names(tbls_plsda) <- paste0("pls_", input$plsda_choice)
                       # - - -
                       tbls_plsda
                     },
                     volc = list(rownames(mSet$analSet$volcano$sig.mat)))
      
      # - - - - - - - 
      # FIX NAMES!!!
      
      tbls_top <- lapply(tbls, function(tbl){
        if(length(tbl) < top){
          tbl
        }else{
          tbl[1:top]
        }
      })
      
      # - - return - -
      tbls_top
    })
    
    names(tables) <- categories
    # - - unlist - -
    
    flattened <<- flattenlist(tables)
    names(flattened) <<- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
    
    circles = length(flattened)
    
    # - - - - - - --
    
    # TODO
    # Hypergeometric testing?
    # SO : 'Calculate venn diagram hypergeometric p value using R'
    
    venn.plot <- VennDiagram::venn.diagram(x = flattened,
                                           filename = NULL)
    
    # - - - - - - - - -
    
    items <- strsplit(as.character(venn.plot), split = ",")[[1]]
    
    circ_values <<- data.frame(
      id = 1:length(grep(items, pattern="polygon"))
      #,value = c(3, 3.1, 3.1, 3.2, 3.15, 3.5)
    )
    
    txt_values <- data.frame(
      id = grep(items, pattern="text"),
      value = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$label))
    )
    
    txt_values$value <- gsub(x = txt_values$value, pattern = "(.*\\.)(.*$)", replacement = "\\2")
    categories <- c(categories, input$rf_choice, input$ls_choice, input$plsda_choice)
    
    x_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$x))
    y_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$y))
    
    x_t = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$x))
    y_t = unlist(lapply(grep(items, pattern="text"), function(i)venn.plot[[i]]$y))
    
    positions_c <- data.frame(
      id = rep(circ_values$id, each = length(x_c)/length(circ_values$id)),
      x = x_c,
      y = y_c
    )
    
    positions_t <- data.frame(
      id = rep(txt_values$id, each = length(x_t)/length(txt_values$id)),
      x = x_t,
      y = y_t
    )
    
    datapoly <- merge(circ_values, positions_c, by=c("id"))
    datatxt <- merge(txt_values, positions_t, by=c("id"))
    
    numbers <- datatxt[!(datatxt$value %in% names(flattened)),]
    headers <- datatxt[(datatxt$value %in% names(flattened)),]
    
    if(circles == 2){
      occur <- table(numbers$y)
      newy <- names(occur[occur == max(occur)])
      # - - -
      numbers$y <- as.numeric(c(newy))
    }
    
    p <- ggplot(datapoly, 
                aes(x = x, 
                    y = y)) + geom_polygon(colour="black", alpha=0.5, aes(fill=id, group=id)) +
      geom_text(mapping = aes(x=x, y=y, label=value), data = numbers, size = 5) +
      geom_text(mapping = aes(x=x, y=y, label=value), data = headers, fontface="bold", size = 7) +
      theme_void() +
      theme(legend.position="none") + 
      scale_fill_gradientn(colours = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](circles)) +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
    
    
    output$venn_plot <- plotly::renderPlotly({
      
      ggplotly(p, tooltip = "label")
      #%>% layout(plot_bgcolor='transparent')
      #%>% layout(paper_bgcolor='transparent')
    
      })
    updateSelectizeInput(session, "intersect_venn", choices = names(flattened))
  })
  
  observeEvent(input$intersect_venn, {

    if(length(input$intersect_venn) > 1){
      venn_overlap <<- Reduce("intersect", lapply(input$intersect_venn, function(x){
        flattened[[x]]
      })
      )
    }
    output$venn_tab <- DT::renderDataTable({
      # -------------
      DT::datatable(data.table(mz = venn_overlap), 
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })
  })
  
  observeEvent(input$venn_tab_rows_selected, {
    curr_cpd <<- venn_overlap[input$venn_tab_rows_selected]
    output$curr_cpd <- renderText(curr_cpd)
  })
  
  
  
  # - - ON CLOSE - -
  
  observe({
    if (input$nav_general == "stop"){
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
                               if(any(!is.na(session_cl))) parallel::stopCluster(session_cl)
                               R.utils::gcDLLs() # flush dlls
                               stopApp() 
                             }else{
                               NULL
                             }})
    }
  })
})
