# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

shinyServer(function(input, output, session) {
  
  # ================================= DEFAULTS ===================================
  
  source('./backend/scripts/joanna/shiny_general.R')

  # set progress bar style to 'old' (otherwise it's not movable with CSS)
  shinyOptions(progress.style="old")
  
  # send specific functions/packages to other threads
  parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
    "mape",
    "flattenlist"
  ))
  parallel::clusterEvalQ(session_cl, library(data.table))

  # create default text objects in UI
  lapply(global$constants$default.text, FUN=function(default){
    output[[default$name]] = renderText(default$text)
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
      print(i)
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
  
  # ===== STARTUP SETTINGS =====
  
  # create listener for what mode we're currently working in (bivariate, multivariate, time series...)
  datamanager <- reactiveValues()
  
  # check if a dataset is already loaded in
  # change mode according to how many levels the experimental variable has
  # change interface based on that
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
  
  # reload multiple ui interfaces if mset is already present
  # reload figures for pca/plsda
  # re-render drop down bars for picking which experimental variable is chosen
  observe({
    if(datamanager$mset_present){
	  # reload pca, plsda (make datamanager do that)
      datamanager$reload <- "pca"
      datamanager$reload <- "plsda"
      # update select input bars with current variable and covariables defined in excel
      updateSelectInput(session, "first_var", selected = mSet$dataSet$cls.name, choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
      updateSelectInput(session, "second_var", choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
      updateSelectInput(session, "subset_var", choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
      # if _T in sample names, data is time series. This makes the time series swap button visible. 
      if(all(grepl(pattern = "_T\\d", x = rownames(mSet$dataSet$norm)))){
        timebutton$status <- "on"
      }else{
        timebutton$status <- "off"
      }
	  # show a button with t-test or fold-change analysis if data is bivariate. hide otherwise.
	  # TODO: add button for anova/other type of sorting...
      if(mSet$dataSet$cls.num == 2 ){
        heatbutton$status <- "ttfc"
      }else{
        heatbutton$status <- NULL
      }
    }else{
	  # hide time series button
      timebutton$status <- "off"
      heatbutton$status <- "asmb"
    }
  })

  # ===== VARIABLE SWITCHER ====
  
  # set default mode for heatmap top hits pick button (tt/fc or asca/meba)
  heatbutton <- reactiveValues(status = "ttfc")
  
  # set default timemode switcher button mode to hidden
  timebutton <- reactiveValues(status = "off")
  
  # render heatmap button
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
  
  # render time series swap button
  output$timebutton <- renderUI({
    if (is.null(timebutton$status)) {
      NULL
    } else{
      switch(timebutton$status, 
             off = NULL,
             on = switchButton(inputId = "timecourse_trigger", # this is used to trigger time series dataswap
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
	  # TODO: use this in venn diagram creation
      mSet$storage[[paste0(as.character(levels(mSet$dataSet$cls)), collapse="-")]] <<- mSet$analSet
      
      # adjust mset design type (necessary for metaboanalystr)
      SetDesignType(mSet, "time")
      
      # rename some factors of interest (your experimental variable, and 'time') as A and B
      facA <- as.factor(mSet$dataSet$covars[,mSet$dataSet$cls.name, with=F][[1]])
      facB <- mSet$dataSet$covars[,"time"][[1]]
      
	  # change mSet experimental factors (these are used in ASCA/MEBA etc.)
      mSet$dataSet$exp.fac <<- as.factor(facA)
      mSet$dataSet$time.fac <<- as.factor(facB)
      mSet$dataSet$facA <<- mSet$dataSet$exp.fac;
      mSet$dataSet$facB <<- mSet$dataSet$time.fac;
      
	  # change interface to timeseries mode (make 'interface' manager do it)
      interface$mode <- "time"
      
	  # change heatmap chooser to asca/meba because those are timeseries-specific
      heatbutton$status <- "asmb"
	  
	  # REMOVE PREVIOUS ANALYSIS TO TRIGGER RELOAD (or the PCA won't reload)
      mSet$analSet <<- NULL
      
    }else{
      # change back to normal mode
      # save previous analyses (should be usable in venn diagram later)
      mSet$storage[[paste0("(timecourse)", paste0(as.character(levels(mSet$dataSet$cls)), collapse="-"))]] <<- mSet$analSet
      # - - - - - - - - - - - -
	  # rename experimental factors
      mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,mSet$dataSet$cls.name, with=F][[1]])
      mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
      # remove old analSet
      mSet$analSet <<- NULL
      
	  # change heatmap button back to t-test/fold-change bivariate mode
	  heatbutton$status <- "ttfc"

      # reset interface
      if(mSet$dataSet$cls.num <= 1){
        interface$mode <- NULL } 
	# check bivariate/multivariate and change interface back
      else if(mSet$dataSet$cls.num == 2){
        interface$mode <- "bivar"}
      else{
        interface$mode <- "multivar"}
    }
  }, ignoreInit = TRUE)
  
  # create interface mode storage object.
  interface <- reactiveValues()
  
  # triggers when the 'change variable' dropdown menu is filled and button is clicked
  observeEvent(input$change_cls, {
    
	# check if previous analysis storage already exists, if not, make it
    if(!("storage" %in% names(mSet))){
      mSet$storage <<- list()
    }
    
    # save previous analyses (should be usable in venn diagram later)
    mSet$storage[[paste0(as.character(levels(mSet$dataSet$cls)), collapse="-")]] <<- mSet$analSet
    
	# change current variable of interest to user pick from covars table
    mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,input$first_var, with=F][[1]])
	
	# adjust bivariate/multivariate (2, >2)...
    mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
    
    # remove old analSet
    mSet$analSet <<- NULL
    
	# adjust name of experimental variable
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
  
  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  observe({
    
	# hide all tabs by default, easier to hide them and then make visible selectively
    hide.tabs <- c("inf", "pca", "plsda", "tt", "fc", "aov", "meba", "asca", "ml", "volc", "heatmap", "enrich")
  
    # check mode of interface (depends on timeseries /yes/no and bivariate/multivariate)
	# then show the relevent tabs
	# TODO: enable multivariate time series analysis
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
      show.tabs <- c("inf") # 'info' tab that loads when no data is loaded currently
    }
	
	# hide all the tabs to begin with
    for(tab in hide.tabs){
      print(tab)
      hideTab(inputId = "statistics", tab, session = session)
    }
    i=1
	# show the relevant tabs
    for(tab in show.tabs){
      print(tab)
      showTab(inputId = "statistics", tab, select = ifelse(i==1, TRUE, FALSE), session = session)
      i = i + 1
    }
  })
  
  # -----------------
  
  # generate positive and negative adduct picker tabs (for csv creation)
  # defaults are in the huge global object :-)
  observe({
    modes = c("pos", "neg")
	lapply(modes, function(mode){
		output[[paste0(mode, "add_tab")]] <- DT::renderDataTable({
			DT::datatable(global$vectors$pos_adducts,
					selection = list(mode = 'multiple', 
					selected = global$vectors[[paste0(mode, "selected_adducts")]], target="row"),
					options = list(pageLength = 5, dom = 'tp'), 
					rownames = F)
		})
	})
  })
  
  # toggles when 'select all adducts' is pressed (filled circle)
  observeEvent(input$sel_all_adducts, {
    global$vectors$neg_selected_adducts <<- c(1:nrow(global$vectors$pos_adducts))
	global$vectors$pos_selected_adducts <<- c(1:nrow(global$vectors$neg_adducts))
  })
  
  # triggers when 'select no adducts' is selected
  observeEvent(input$sel_no_adducts, {
    global$vectors$neg_selected_adducts <<- c(0)
	global$vectors$pos_selected_adducts <<- c(0)
  })

  # triggers when common adducts are to be selected
  observeEvent(input$sel_comm_adducts, {
    global$vectors$neg_selected_adducts <<- c(1:3, nrow(global$vectors$pos_adducts))
	global$vectors$pos_selected_adducts <<- c(1, 2, 14:15, nrow(global$vectors$neg_adducts))
  })
  
  # if general tab selection changes, load package table
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
  
  # triggers when 'update' button is pressed on the packages table tab
  observeEvent(input$update_packages, {
	# use pacman package to update packages (should swap between normal installed and bioconductor)
    pacman::p_load(char = global$constants$packages, update = T, character.only = T)
    # refresh package list 
    pkg_tbl <- get.package.table() 
	# reload table
    output$package_tab <- DT::renderDataTable({
      # - - - - - - - - - - -
      DT::datatable(pkg_tbl,
                    selection = 'none',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(10, 20, 30), pageLength = 10), 
                    rownames = F)
    })
  })

  # listener for all the file pickers, trigger a window if they are clicked
  observe({
    shinyFileChoose(input, 'outlist_pos', roots=global$paths$volumes, filetypes=c('csv'))
    shinyFileChoose(input, 'outlist_neg', roots=global$paths$volumes, filetypes=c('csv'))
    shinyFileChoose(input, 'excel', roots=global$paths$volumes, filetypes=c('xls', 'xlsm', 'xlsx'))
    shinyFileChoose(input, 'database', roots=global$paths$volumes, filetypes=c('sqlite3', 'db', 'sqlite'))
    shinyFileChoose(input, 'taskbar_image_path', roots=global$paths$volumes, filetypes=c('png', 'jpg', 'jpeg', 'bmp'))
  })
  
  # observes if a new taskbar image is chosen by user
  observe({
    # - - - - 
    if(!is.list(input$taskbar_image_path)) return() # if nothing is chosen, do nothing
    img_path <- parseFilePaths(global$paths$volumes, input$taskbar_image_path)$datapath
    new_path <- file.path(getwd(), "www", basename(img_path)) # set path to copy to
    
	# copy image to the www folder
    if(img_path != new_path) file.copy(img_path, new_path, overwrite = T)
    # - - -
	# render taskbar image preview
    output$taskbar_image <- renderImage({
      list(src = new_path, 
           width = 120,
           height = 120,
           style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;")
    }, deleteFile = FALSE)
    # change chosen taskbar image in user option file
    setOption('user_options.txt', 'taskbar_image', basename(new_path))
  })
  
  # observes if user is choosing a different database storage folder
  observe({  
	# trigger window
    shinyDirChoose(input, "get_db_dir", 
                   roots=global$paths$volumes, 
                   session = session)
         
    if(typeof(input$get_db_dir) != "list") return() # if nothing selected or done, ignore

	# parse the file path given based on the possible base folders (defined in global)
    given_dir <- parseDirPath(global$paths$volumes, 
                              input$get_db_dir)
							  
    if(is.null(given_dir)) return()
    # change db storage directory in user options file
    setOption("user_options.txt", "db_dir", given_dir)
	
	# render current db location in text
    output$curr_db_dir <- renderText({getOptions('user_options.txt')$db_dir})
  })
  
  # see above, but for working directory. CSV/DB files with user data are stored here.
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
  
  # triggers if user changes their current project name
  observeEvent(input$set_proj_name, {
    proj_name <<- input$proj_name
    if(proj_name == "") return(NULL) # if empty, ignore
	# change path of current db in global
    global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(proj_name,".db", sep=""))
    # change project name in user options file
    setOption("user_options.txt", "proj_name", proj_name)
	# print the changed name in the UI
    output$proj_name <<- renderText(proj_name)
	# change path CSV should be / is saved to in session
    global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
    
  })
  
  # change ppm accuracy, ONLY USEFUL if loading in from CSV
  # TODO: make it possible to change this and re-make user database (mzranges table specifically)
  observeEvent(input$set_ppm, {
    ppm <<- input$ppm
	# show ppm amount in UI
    output$ppm <- renderText(ppm)
    # change in options file
    setOption("user_options.txt", "ppm", ppm)
  })
  
  
  # triggers when a new color spectrum is chosen
  observeEvent(input$color_ramp,{
	# render preview plot
    output$ramp_plot <- plotly::renderPlotly({

	  # change the current spectrum in global
      global$functions$color.functions[[getOptions("user_options.txt")$gspec]] <<- global$functions$color.functions[[input$color_ramp]]
      # change the current spectrum in user options file
	  setOption("user_options.txt", "gspec", input$color_ramp)
      
	  # create plot background (no grid, no lines, just color ;) )
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        titlefont = list(size = 20)
      )
      
      # re-render preview plot with the new options (general heatmap using R standard volcano dataset)
      plotly::plot_ly(z = volcano, 
                      colors = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](100), 
                      type = "heatmap",
                      showscale=FALSE)  %>%
        layout(xaxis = ax, yaxis = ax)
    })
  })
  
  # triggers when new plot theme is picked
  observeEvent(input$ggplot_theme,{
    
	# change default plot theme in user settings
    setOption("user_options.txt", "gtheme", input$ggplot_theme)
    
	# change preview plot (uses mtcars default R dataset)
    output$ggplot_theme_example <- renderPlot({
      p <- ggplot(mtcars) + geom_boxplot(aes(x = wt, y = mpg,
                                           colour = factor(gear)))
      p + global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]]()
    })
  })
  
  # triggers when changes in interface aesthetics are applied
  observeEvent(input$change_css, {

	# set default user color options
    setOption("user_options.txt", "col1", input$bar.col.1)
    setOption("user_options.txt", "col2", input$bar.col.2)
    setOption("user_options.txt", "col3", input$bar.col.3)
    setOption("user_options.txt", "col4", input$bar.col.4)
    
	# set default user font options
    setOption("user_options.txt", "font1", input$font.1)
    setOption("user_options.txt", "font2", input$font.2)
    setOption("user_options.txt", "font3", input$font.3)
    setOption("user_options.txt", "font4", input$font.4)
    
	# set default user font size options
    setOption("user_options.txt", "size1", input$size.1)
    setOption("user_options.txt", "size2", input$size.2)
    setOption("user_options.txt", "size3", input$size.3)
    setOption("user_options.txt", "size4", input$size.4)
  })

  # adduct table editing from settings tab
  
  values = reactiveValues()
  
  # TODO: fix and re-docment this
  observeEvent(input$import_adducts, {
    req(input$add_tab)
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
  
  # uses rhandsontable for live table editing...
  # TODO: fix
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
  
  # everything below uses the dblist defined in global
  # as well as the logos defined here
  # if you add a db, both the name and associated logo need to be added
  
  # create checkcmarks if database is present
  lapply(global$vectors$db_list, FUN=function(db){
	# creates listener for if the 'check db' button is pressed
    observeEvent(input[[paste0("check_", db)]],{
	  # see which db files are present in folder
      db_folder_files <- list.files(getOptions("user_options.txt")$db_dir)
      is.present <- paste0(db, ".full.db") %in% db_folder_files
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
  lapply(global$vectors$db_list, FUN=function(db){
    observeEvent(input[[paste0("build_", db)]], {
      # ---------------------------
      library(RCurl)
      library(XML)
      library(SPARQL)
      # ---------------------------
      withProgress({
        
		# send necessary functions and libraries to parallel threads
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
        shiny::setProgress(session = session, 0.1)
		# build base db (differs per db, parsers for downloaded data)
        build.base.db(db,
                      outfolder = getOptions("user_options.txt")$db_dir, 
                      cl = session_cl)
        shiny::setProgress(session = session, 0.5)
		# extend base db (identical per db, makes adduct and isotope variants of downloaded compounds)
        build.extended.db(db, 
                          outfolder = getOptions("user_options.txt")$db_dir,
                          adduct.table = adducts, 
                          cl = session_cl, 
                          fetch.limit = 500) #TODO: figure out the optimal fetch limit...
      })
    })
  })
  
  # triggers if isotope scoring is clicked after finding db matches
  observeEvent(input$score_iso, {

	# check if the matches table even exists
    if(!data.table::is.data.table(global$tables$last_matches)) return(NULL)
    
	# check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
    if("score" %in% colnames(global$tables$last_matches)){
      global$tables$last_matches <<- global$tables$last_matches[,-"score"]
    }
    
	# get table including isotope scores
	# as input, takes user method for doing this scoring
    score_table <- score.isos(global$paths$patdb, method=input$iso_score_method, inshiny=T) 
    
	# update the match table available to the rest of metaboshiny
    global$tables$last_matches <<- global$tables$last_matches[score_table, on = c("baseformula", "adduct")]
    
	# re-render match table
    output$match_tab <-DT::renderDataTable({
	  # don't show some columns but keep them in the original table, so they can be used
	  # for showing molecule descriptions, structure
      DT::datatable(global$tables$last_matches[,-c("description","structure", "baseformula", "dppm")],
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })  
  })
  
  # triggers when user wants to create database from .db and excel or 2 csv files and excel
  observeEvent(input$create_db,{
    
	# update the path to patient db
    global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name, ".db"))
    
    withProgress({
      
      shiny::setProgress(session=session, value= .1)
      
      switch(input$new_proj,
			 # if loading in a .db file... (FAST, MOSTLY FOR ADMINS USING HPC)
             `From DB` = {
               
               req(input$database, input$excel)
               
			   # get the db and excel path from the UI elemnts
               db_path <- parseFilePaths(global$paths$volumes, input$database)$datapath
               excel_path <- parseFilePaths(global$paths$volumes, input$excel)$datapath
               
			   # copy the user selected db to the processing folder under proj_name renaming
               file.copy(db_path, global$paths$patdb, overwrite = T)
               
               shiny::setProgress(session=session, value= .30)
               
			   # add excel file to this database
               exp_vars <- load.excel(excel_path, global$paths$patdb)
               
               shiny::setProgress(session=session, value= .60)
               
             },
			 # if loading in .csv files...
             `From CSV` = {
               
               req(input$outlist_neg, input$outlist_pos, input$excel)
               
			   # build patient db from csv files with a given ppm error margin
               build.pat.db(global$paths$patdb,
                            ppm = ppm,
                            pospath = parseFilePaths(global$paths$volumes, input$outlist_pos)$datapath,
                            negpath = parseFilePaths(global$paths$volumes, input$outlist_neg)$datapath,
                            overwrite = T)
               
               shiny::setProgress(session=session, value= .95,message = "Adding excel sheets to database...")
               
			   # add excel file to .db file generated in previous step
               exp_vars <<- load.excel(parseFilePaths(global$paths$volumes, input$excel)$datapath, global$paths$patdb)}
             )
      })
  })
  
  # imports existing db file
  # TODO: is deprecated, fix!!
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
  
  # imports existing csv file
  observeEvent(input$import_csv, {
    req(input$pat_csv)
	# change path to current csv file to user given path
    global$paths$csv_loc <<- input$pat_csv$datapath
    
	# show checkmark underneath select csv button
    output$csv_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  

  # is triggered when the create csv button is clicked
  observeEvent(input$create_csv, {

    withProgress({
      
      shiny::setProgress(session = session, value= 1/4)

	  # create csv table from patient database and user chosen settings in that pane
      tbl <- get.csv(global$paths$patdb,
                     group_adducts = if(length(global$vectors$add_search_list) == 0) F else T, # group by addicts?
                     groupfac = input$group_by, # group by mz or formula
                     which_dbs = global$vectors$add_search_list, # used databases
                     which_adducts = selected_adduct_list # used adducts
      )
      
	  # check if any samples are duplicated in the resulting table 
	  # IF SO: rename to time series data format because that is the only one that should have duplicate sample names
      if(any(duplicated(tbl$sample))){
        tbl$sample <- paste0(tbl$sample, "_T", tbl$time)
        show.times = T # make sure the 'times' column shows in the csv making result table
      }else{
        show.times = F
      }
      
      shiny::setProgress(session=session, value= 2/4)
	  
	  # change location of csv in global
      global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
	  
	  # write csv file to new location
      fwrite(tbl, global$paths$csv_loc, sep="\t")

	  # find the experimental variables by checking which column names wont convert to numeric
      as.numi <- as.numeric(colnames(tbl)[1:100])
      exp.vars <- which(is.na(as.numi))
      
      shiny::setProgress(session=session, value= 3/4)
      
	  # render overview table
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
  
  # triggers when 'get options' is clicked in the normalization pane
  observeEvent(input$check_csv, {
    # ----------------------
    
	# read in csv 
	# TODO: only read in the first x -rows- to save time??
    csv <- fread(global$paths$csv_loc, 
                 header = T)
    
	# find experimental variables by checking which wont convert to numeric
    as.numi <- as.numeric(colnames(csv)[1:100])
    exp.vars <- which(is.na(as.numi))
    
	# get the names of those experimental variables
    opts <<- colnames(csv[,..exp.vars])
    
	# get columns that can be used for batch correction (need to be non-unique)
    batch <<- which(sapply(exp.vars, function(x) length(unique(csv[,..x][[1]])) < nrow(csv)))
    
	# get sample columns
    numi <<- which(sapply(exp.vars, function(x) is.numeric(csv[,..x][[1]])))

	# update the possible options in the UI
    updateSelectInput(session, "samp_var",
                      choices = opts[numi])
    updateSelectizeInput(session, "batch_var",
                         choices = opts[batch],
                         options = list(maxItems = 3L - (length(input$batch_var)))
    )
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
  
  # this triggers when the user wants to normalize their data to proceed to statistics
  observeEvent(input$initialize,{
  
    withProgress({

      shiny::setProgress(session=session, value= .1)

      # read in original CSV file
      csv_orig <- fread(global$paths$csv_loc, 
                        data.table = TRUE,
                        header = T)
      
      # create empty mSet with 'stat' as default mode
	  mSet <- InitDataObjects("pktable",
                              "stat",
                              FALSE)
      
      # convert all 0's to NA so metaboanalystR will recognize them
	  csv_orig[,(1:ncol(csv_orig)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)})]
      
      # remove whitespace
	  csv_orig$sample <- gsub(csv_orig$sample, pattern=" ", replacement="")
      
	  # find experimental variables by converting to numeric
      as.numi <- as.numeric(colnames(csv_orig)[1:100])
      exp.vars <- which(is.na(as.numi))
      
      # load batch variable chosen by user
      batches <- input$batch_var
      
      # locate qc containing rows in csv
	  qc.rows <- which(grepl("QC", csv_orig$sample))
      
	  # for the non-qc samples, check experimental variables. Which have at least 2 different factors, but as little as possible?
      unique.levels <- apply(csv_orig[!qc.rows,..exp.vars, with=F], MARGIN=2, function(col){
        lvls <- levels(as.factor(col))
        # - - - - - -
        length(lvls)
      })
      
	  # use this as the default selected experimental variable (user can change later)
      which.default <- unique.levels[which(unique.levels == min(unique.levels[which(unique.levels> 1)]))][1]
      condition = names(which.default)
      
      # =================================|
      
	  # if nothing is selected for batch, give empty
      if(is.null(batches)) batches <- ""
      
	  # only turn on batch correction if user says so
      batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE
      
	  # if 'batch' is selected, 'injection' is often also present
	  # TODO: i can imagine this doesn't work for all  users, please disable this...
      if("batch" %in% batches){ 
        batches = c(batches, "injection")
      }
      
	  # get the part of csv with only the experimental variables
      first_part <- csv_orig[,..exp.vars, with=FALSE]
	  
	  # set NULL or missing levels to "unknown"
      first_part[first_part == "" | is.null(first_part)] <- "unknown"
      
	  # re-make csv with the corrected data
      csv <- cbind(first_part[,-c("label")], # if 'label' is in excel file remove it, it will clash with the metaboanalystR 'label'
                        "label" = first_part[,..condition][[1]], # set label as the initial variable of interest
                        csv_orig[,-..exp.vars,with=FALSE])
      

      # remove outliers by making a boxplot and going from there
      if(input$remove_outliers){
        sums <- rowSums(csv[,-exp.vars,with=FALSE],na.rm = TRUE)
        names(sums) <- csv$sample
        outliers = c(car::Boxplot(as.data.frame(sums)))
        csv <- csv[!(sample %in% outliers),]
      } 
      
      # remove peaks that are missing in all 
      csv <- csv[,which(unlist(lapply(csv, function(x)!all(is.na(x))))),with=F]
      
      # remove samples with really low numbers of peaks
	  complete.perc <- rowMeans(!is.na(csv))
      keep_samps <- csv$sample[which(complete.perc > .2)]
      csv <- csv[sample %in% keep_samps,]
      
	  # also remove them in the table with covariates
      covar_table <- first_part[sample %in% keep_samps,]
      
	  # if the experimental condition is batch, make sure QC samples are not removed at the end for analysis
	  # TODO: this is broken with the new system, move this to the variable switching segment of code
      batchview = if(condition == "batch") TRUE else FALSE
      
	  # if QC present, only keep QCs that share batches with the other samples (may occur when subsetting data/only loading part of the samples)
      if(any(grepl("QC", csv$sample))){
        samps <- which(!grepl(csv$sample, pattern = "QC"))
        batchnum <- unique(csv[samps, "batch"][[1]])
        keep_samps_post_qc <- covar_table[which(covar_table$batch %in% batchnum),"sample"][[1]]
        covar_table <- covar_table[which(covar_table$batch %in% batchnum),]
        csv <- csv[which(csv$sample %in% keep_samps_post_qc),-"batch"]
      }
      
	  # rename time column or metaboanalyst won't recognize it
      colnames(csv)[which( colnames(csv) == "time")] <- "Time"
      
      # deduplicate columns
	  # TODO: remove the source of the duplicated columns ealier, might already be fixed
	  remove = which( duplicated( t(csv[,..exp.vars] )))
      colnms <- colnames(csv)[remove]
      remove.filt <- setdiff(colnms,c("sample","label"))
      csv <- csv[ , -remove.filt, with = FALSE ]

      # find experimental variables
	  as.numi <- as.numeric(colnames(csv)[1:100])
      exp.vars <- which(is.na(as.numi))
      
      # remove all except sample and time in saved csv
      exp_var_names <- colnames(csv)[exp.vars]
      keep_cols <-  c("sample", "label")
      remove <- which(!(exp_var_names %in% keep_cols))
      
      # define location to write processed csv to
      csv_loc_final <- gsub(pattern = "\\.csv", replacement = "_no_out.csv", x = global$paths$csv_loc)
      
	  # remove file if it already exists
      if(file.exists(csv_loc_final)) file.remove(csv_loc_final)
      
      # write new csv to new location
      fwrite(csv[,-remove,with=F], file = csv_loc_final)
      
	  # rename row names of covariant table to the sample names
      rownames(covar_table) <- covar_table$sample

	  # load new csv into empty mSet!
      mSet <- Read.TextData(mSet, 
                            filePath = csv_loc_final, 
                            "rowu")  # rows contain samples
      
	  # add covars to the mSet for later switching and machine learning
      mSet$dataSet$covars <- covar_table
      
      # sanity check data 
	  mSet <- SanityCheckData(mSet)
      
	  # remove metabolites with more than user defined perc missing
      mSet <- RemoveMissingPercent(mSet, 
                                    percent = input$perc_limit/100)
      
	  # missing value imputation
      if(input$miss_type != "none"){
        if(input$miss_type == "rowmin"){ # use sample minimum
          new.mat <- apply(mSet$dataSet$preproc, 1, function(x) {
            if (sum(is.na(x)) > 0) {
              x[is.na(x)] <- min(x, na.rm = T)/2
            }
            x
          })
          mSet$dataSet$procr <- t(new.mat)
        }
        else if(input$miss_type == "pmm"){ # use predictive mean matching
		# TODO: re-enable, it's very slow
          require(mice)
          base <- mSet$dataSet$preproc
          imp <- mice::mice(base, printFlag = TRUE)
        }else if(input$miss_type == "rf"){ # random forest
          samples <- rownames(mSet$dataSet$preproc)
          
		  # convert all to as numeric
		  # TODO: remove, should be automatic
          w.missing <- mSet$dataSet$preproc
          w.missing <- apply(w.missing, 2, as.numeric)

		  # register other threads as parallel threads
          doParallel::registerDoParallel(session_cl)
          
		  # set amount of tries (defined by missforest package)
          auto.mtry <- floor(sqrt(ncol(mSet$dataSet$preproc)))
          
          mtry <- ifelse(auto.mtry > 100, 100, auto.mtry)
          
		  # impute missing values with random forest
          imp <- missForest::missForest(w.missing, 
                                        parallelize = "variables", # parallelize over variables, 'forests' is other option
                                        #verbose = T,
                                        #ntree = 10,
                                        mtry = mtry)
          
          mSet$dataSet$procr <- imp$ximp
          rownames(mSet$dataSet$procr) <- rownames(mSet$dataSet$preproc)
          # - - - - - - - - - - - - 
        }else{
		  # use built in imputation methods, knn means etc.
          mSet <- ImputeVar(mSet,
                             method = # "knn"
                              input$miss_type
                              )
        }
      }
      # ------------------------
      
	  # if user picked a data filter
      if(input$filt_type != "none"){
        # filter dataset
		# TODO; add option to only keep columns that are also in QC ('qcfilter'?)
        mSet <- FilterVariable(mSet,
                                filter = input$filt_type,
                                qcFilter = "F",
                                rsd = 25)
      }
      # ------------------------------------
      
      shiny::setProgress(session=session, value= .2)
      
	  # if normalizing by a factor, do the below
      if(input$norm_type == "SpecNorm"){
        norm.vec <- mSet$dataSet$covars[match(mSet$dataSet$covars$sample,
                                               rownames(mSet$dataSet$preproc)),][[input$samp_var]]
        norm.vec <- scale(x = norm.vec, center = 1) # normalize scaling factor
        
      }else{
        norm.vec <- rep(1, length(mSet$dataSet$cls)) # empty
      }
      
	  # normalize dataset with user settings(result: mSet$dataSet$norm)
      mSet <- Normalization(mSet,
                             rowNorm = input$norm_type,
                             transNorm = input$trans_type,
                             scaleNorm = input$scale_type,
                             ref = input$ref_var)
      
      shiny::setProgress(session=session, value= .4)
      
	  # get sample names
      smps <- rownames(mSet$dataSet$norm)
      
	  # get which rows are QC samples
      qc_rows <- which(grepl(pattern = "QC", x = smps))
      
      # if at least one row has a QC in it, batch correct
      has.qc <- length(qc_rows) > 0
      
	  # lowercase all the covars table column names
      colnames(mSet$dataSet$covars) <- tolower(colnames(mSet$dataSet$covars))
      
	  # remove special characters in sample names 
      # TODO: find a better fix for samples w/ special characters in name...
      mSet$dataSet$covars$sample <- gsub(mSet$dataSet$covars$sample, pattern = "\\(|\\)", replacement="")
      
      if(batch_corr){
        
        if("batch" %in% input$batch_var & has.qc){
          
		  # get batch for each sample
          batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"batch"][[1]]))
		  
		  # get injection order for samples
          seq.idx = as.numeric(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$sample),"injection"][[1]])
          
          # go through all the metabolite columns
          corr_cols <- pbapply::pblapply(1:ncol(mSet$dataSet$norm), function(i){
            # fetch non-corrected values
            vec = mSet$dataSet$norm[,i]
            # correct values using QCs and injectiono rder
            corr_vec = BatchCorrMetabolomics::doBC(Xvec = as.numeric(vec), 
                                                   ref.idx = as.numeric(qc_rows), 
                                                   batch.idx = batch.idx,
                                                   seq.idx = seq.idx,
                                                   result = "correctedX",
                                                   minBsamp = 1) # at least one QC necessary
            corr_vec
          })
          
		  # cbind the corrected columns to re-make table
          qc_corr_matrix <- as.data.frame(do.call(cbind, corr_cols))
          # fix rownames to old rownames
          colnames(qc_corr_matrix) <- colnames(mSet$dataSet$norm)
          rownames(qc_corr_matrix) <- rownames(mSet$dataSet$norm)
          # save to mSet
          mSet$dataSet$norm <- as.data.frame(qc_corr_matrix)
          
        }
        
		# remove QC samples if user doesn't use batch as condition
        if(!batchview & has.qc){
          mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
          mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
          mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
          mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
        }
        
		# check which batch values are left after initial correction
        left_batch_vars <- grep(input$batch_var, 
                                pattern =  ifelse(has.qc, "batch|injection|sample", "injection|sample"),
                                value = T,
                                invert = T)
        
        if(length(left_batch_vars) > 2){ 
          NULL  # no option for more than 2 other batch variables yet
        } else if(length(left_batch_vars) == 0){
          NULL # if none left, continue after this
        } else{
          # get sample names and classes
          smp <- rownames(mSet$dataSet$norm)
          exp_lbl <- mSet$dataSet$cls
          
		  # create csv for comBat
          csv <- as.data.table(cbind(sample = smp, 
                                     label = mSet$dataSet$cls,
                                     mSet$dataSet$norm))
          
		  # transpose for combat
          csv_edata <-t(csv[,!c(1,2)])
          colnames(csv_edata) <- csv$sample
          
          if(length(left_batch_vars) == 1){
			# create a model table
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample),left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = c(0),
                                    outcome = as.factor(exp_lbl)) 
			# batch correct with comBat
            batch_normalized= t(sva::ComBat(dat=csv_edata,
                                            batch=csv_pheno$batch1
                                            #mod=mod.pheno,
                                            #par.prior=TRUE
            ))
			# fix row names
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }else{
			# create a model table
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample), left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$sample), left_batch_vars[2], with=FALSE][[1]],
                                    outcome = as.factor(exp_lbl))
			# batch correct with limma and two batches
            batch_normalized = t(limma::removeBatchEffect(x = csv_edata,
                                                          #design = mod.pheno,
                                                          batch = csv_pheno$batch1,
                                                          batch2 = csv_pheno$batch2))
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }
          # save normalized table to mSet
          mSet$dataSet$norm <- as.data.frame(batch_normalized)
        }
      } else{
		# if qcs presnt and user doesn't want to analyse qc samples
        if(!batchview & has.qc){
		# remove QC rows and associated data from mSet
          mSet$dataSet$norm <- mSet$dataSet$norm[-qc_rows,]
          mSet$dataSet$cls <- mSet$dataSet$cls[-qc_rows, drop = TRUE]
          mSet$dataSet$covars <- mSet$dataSet$covars[-grep("QC", mSet$dataSet$covars$sample),]
          mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
        }
      }
      
      # make sure covars order is consistent with mset$..$norm order
      mSet$dataSet$covars <- mSet$dataSet$covars[match(rownames(mSet$dataSet$norm), 
                                                        mSet$dataSet$covars$sample),]
      
	  # set name of variable of interest
      mSet$dataSet$cls.name <- condition
      
      # save mset to global environment
      mSet <<- mSet
      
	  # tell datamanager that an mset is present 
      datamanager$mset_present = TRUE
      
	  # change interface mode based on the current condition
      if(mSet$dataSet$cls.num <= 1){
        interface$mode <- NULL } 
      else if(mSet$dataSet$cls.num == 2){
        interface$mode <- "bivar"}
      else{
        interface$mode <- "multivar"}
      
      # - - - - - - - - - -
      
      shiny::setProgress(session=session, value= .5)
      
      shiny::setProgress(session=session, value= .6)
	  
	  # generate summary plots and render them in UI
	  
      varNormPlots <- ggplotNormSummary(mSet)
      output$var1 <- renderPlot(varNormPlots$tl)
      output$var2 <- renderPlot(varNormPlots$bl)
      output$var3 <- renderPlot(varNormPlots$tr)
      output$var4 <- renderPlot(varNormPlots$br)
      
      sampNormPlots <-  ggplotSampleNormSummary(mSet)
      output$samp1 <- renderPlot(sampNormPlots$tl)
      output$samp2 <- renderPlot(sampNormPlots$bl)
      output$samp3 <- renderPlot(sampNormPlots$tr)
      output$samp4 <- renderPlot(sampNormPlots$br)
      shiny::setProgress(session=session, value= .8)
      
	  # save the used adducts to mSet
      mSet$dataSet$adducts <<- selected_adduct_list
      shiny::setProgress(session=session, value= .9)
      
    })
  })
  
  # triggered when user enters the statistics tab
  observeEvent(input$statistics, {
    
	# check if we're in the right tab, otherwise abort
    if(input$nav_general != "analysis") return(NULL)
    
	# check if an mset is present, otherwise abort
    if(!exists("mSet")) return(NULL)
    
    # depending on the present tab, perform analyses accordingly
    switch(input$statistics,
           pca = {
             if(!"pca" %in% names(mSet$analSet)){ # if PCA already has been done, don't redo it
               withProgress({
                 mSet <<- PCA.Anal(mSet) # perform PCA analysis
               })
             }
             datamanager$reload <- "pca" # reload pca plots
             },
           meba = {
             if("MB" %not in% names(mSet$analSet)){ # if already done, don't redo
               mSet <<- performMB(mSet, 10) # perform MEBA analysis
             }
			 # render results table for UI
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
             if("asca" %not in% names(mSet$analSet)){ # if already done, don't redo
			   # perform asca analysis
               mSet <<- Perform.ASCA(mSet, 1, 1, 2, 2)
               mSet <<- CalculateImpVarCutoff(mSet, 0.05, 0.9)
             }
             output$asca_tab <-DT::renderDataTable({ # render results table for UI
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
               
			   # rendering within output object, so it reloads when changing colours etc.
               output$heatmap <- plotly::renderPlotly({
                 
				 # change top hits used in heatmap depending on time series / bivariate / multivariate mode
				 # reordering of hits according to most significant at the top
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
                 
				 # check top x used (slider bar in UI), if more than total matches use total matches
                 topn = if(length(tbl[[used.values]]) < input$heatmap_topn) length(tbl[[used.values]]) else input$heatmap_topn
                 mzorder <- order(tbl[[used.values]], decreasing = decreasing)
                 mzsel <- rownames(tbl)[mzorder][1:topn]
                 
				 # reorder matrix used
                 x <- mSet$dataSet$norm[,mzsel]
                 final_matrix <<- t(x) # transpose so samples are in columns
                 
				 # check if the sample order is correct - mSet$..$ norm needs to match the matrix
                 sample_order <- match(colnames(final_matrix), rownames(mSet$dataSet$norm))
                 
                 if(timebutton$status == "on"){ # check if time series
                   if(input$timecourse_trigger){
					 # create convenient table with the ncessary info
                     translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$exp.fac[sample_order], Time=mSet$dataSet$time.fac[sample_order])
                     hmap.lvls <- c(levels(mSet$dataSet$exp.fac), levels(mSet$dataSet$time.fac))

					 # reorder first by time, then by sample
                     split.translator <- split(translator, by = c("Time"))
                     split.translator.ordered <- lapply(split.translator, function(tbl) tbl[order(tbl$Group)])
                     translator <- rbindlist(split.translator.ordered)
                     
					 # ensure correct sample order
                     final_matrix <<- final_matrix[,match(translator$Sample, colnames(final_matrix))]
                     
					 # disable automatic ordering of samples through clustering
                     my_order=F
                     
                   }else{
					 # no complicated reordering necessary
                     translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$cls[sample_order])
                     hmap.lvls <- levels(mSet$dataSet$cls)
                     my_order = T # enable sorting through dendrogram
                   }
                 }else{
				   # no complicated reordering necessary
                   translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$cls[sample_order])
                   hmap.lvls <- levels(mSet$dataSet$cls)
                   my_order = T # enable sorting through dendrogram
                 } 

				 # create name - to - color mapping vector for the plotting functions
                 color.mapper <- {
                  classes <- hmap.lvls
                  cols <- sapply(1:length(classes), function(i) global$vectors$mycols[i]) # use user-defined colours
                  names(cols) <- classes
                  # - - -
                  cols
                 }
                
				 # create heatmap object :- )
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
                 
				 # save the order of mzs for later clicking functionality
                 hmap_mzs <<- hmap$x$layout$yaxis3$ticktext
                 
                 # - - 
                 
                 hmap
                 
               })               
             })
           },
           tt = {
             if(!"tt" %in% names(mSet$analSet)){ # if already done, don't redo
               withProgress({
                 mSet <<- Ttests.Anal(mSet,
                                      nonpar = FALSE, 
                                      threshp = 0.05, # TODO: make the threshold user defined...
                                      paired = FALSE,
                                      equal.var = TRUE
                                      #  multicorr = "BH"
                 )
               })
             }
			 # save results to table
             res <<- mSet$analSet$tt$sig.mat 
             if(is.null(res)) res <<- data.table("No significant hits found")
			 # render results table for UI
             output$tt_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
			 # render manhattan-like plot for UI
             output$tt_overview_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotTT(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           },
           fc = {
             if(!"fc" %in% names(mSet$analSet)){ # if already done, don't redo
               withProgress({
                 mSet <<- FC.Anal.unpaired(mSet,
                                           2.0, # TODO: make this threshold user defined
                                           1)                 
               })
             }
			 # save results table
             res <<- mSet$analSet$fc$sig.mat 
			 # if none found, give the below table...
             if(is.null(res)) res <<- data.table("No significant hits found")
			 # render result table for UI
             output$fc_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
			 # render manhattan-like plot for UI
             output$fc_overview_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotFC(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           },
           aov = {
             if(!"aov" %in% names(mSet$analSet)){ # if done, don't redo
               withProgress({
                 mSet <<- ANOVA.Anal(mSet, thresh=0.05,nonpar = F) # TODO: make threshold user-defined
               # TODO: Fix ANOVA2, this works for time series data!! error likely due to something missing in conversion from stat -> time mode
               # mSet <<- ifelse(timebutton_trigger,
               #                 {mSet$dataSet$facA <- mSet$dataSet$exp.fac;
               #                  mSet$dataSet$facB <- mSet$dataSet$time.fac;
               #                  ANOVA2.Anal(mSet, 0.05, "fdr", "time", 3, 1)
               #                   },
               #                 ANOVA.Anal(mSet, thresh=0.05,nonpar = F))
               })
             }
			 # render results table for UI
             output$aov_tab <-DT::renderDataTable({
               DT::datatable(if(is.null(mSet$analSet$aov$sig.mat)) data.table("No significant hits found") else mSet$analSet$aov$sig.mat, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
           },
           volc = {
             if(!"volc" %in% names(mSet$analSet)){ # if done, don't redo
               withProgress({
                 mSet <<- Volcano.Anal(mSet,FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw") # TODO: make thresholds user-defined
               })
             }
			 # render results table (currently not visible...)
             output$volc_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$volc$sig.mat, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
			 # render volcano plot with user defined colours
             output$volc_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotVolc(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           }
    )
  })
  
  # triggers when the 'go' button is pressed on the PLS-DA tab
  observeEvent(input$do_plsda, {
    
    require(e1071)
    
	# depending on type, do something else
	# TODO: enable sparse and orthogonal PLS-DA
    switch(input$plsda_type,
           normal={
               require(caret)
               mSet <<- PLSR.Anal(mSet) # perform pls regression
               mSet <<- PLSDA.CV(mSet, methodName=if(nrow(mSet$dataSet$norm) < 50) "L" else "T",compNum = 3) # cross validate
               mSet <<- PLSDA.Permut(mSet,num = 300, type = "accu") # permute
               })
	# reload pls-da plots
    datamanager$reload <- "plsda"
  })
  
  # preload pca/plsda
  observe({
    if(is.null(datamanager$reload)){
      NULL # if not reloading anything, nevermind
    }else{
     switch(datamanager$reload,
            "pca" = {
              if("pca" %in% names(mSet$analSet)){
                # create PCA legend plot
				# TODO: re-enable this plot, it was clickable so you could filter out certain groups
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
			  # render PCA variance per PC table for UI
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
              # render PCA loadings tab for UI
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
			  # chekc which mode we're in
              mode <- if("timecourse_trigger" %in% names(input)){
                if(input$timecourse_trigger){ # if time series mode
                  "ipca" # interactive PCA (old name, i like tpca more :P )
                }else{
                  "pca" # normal pca
                } 
              }else{
                "pca"
              }
              
              # - - - - -
              if(input$pca_2d3d){ # check if switch button is in 2d or 3d mode
                # render 2d plot
                output$plot_pca <- plotly::renderPlotly({
                  plotPCA.2d(mSet, global$vectors$mycols,
                             pcx = input$pca_x,
                             pcy = input$pca_y, mode = mode,
                             shape.fac = input$second_var)
                })
              }else{
                # render 3d plot
                output$plot_pca <- plotly::renderPlotly({
                  plotPCA.3d(mSet, global$vectors$mycols,
                             pcx = input$pca_x,
                             pcy = input$pca_y,
                             pcz = input$pca_z, mode = mode,
                             shape.fac = input$second_var)
                })
              }
              }else{NULL} # do nothing
            },
            "plsda" = {
              
              if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
                
			  # render cross validation plot
              output$plsda_cv_plot <- renderPlot({
                ggPlotClass(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
              })
			  # render permutation plot
              output$plsda_perm_plot <- renderPlot({
                ggPlotPerm(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
              })
			  # render table with variance per PC
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
			  # render table with PLS-DA loadings
              output$plsda_load_tab <-DT::renderDataTable({
                plsda.loadings <- mSet$analSet$plsda$vip.mat
                colnames(plsda.loadings) <- paste0("PC", c(1:ncol(plsda.loadings)))
                # -------------
                DT::datatable(plsda.loadings[, c(input$plsda_x, input$plsda_y, input$plsda_z)],
                              selection = 'single',
                              autoHideNavigation = T,
                              options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
              })
			  # see PCA - render 2d or 3d plots, just with plsda as mode instead
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
              }else{NULL}
            })
      # - - - - 
      datamanager$reload <- NULL # set reloading to 'off'
    }
  })


  # triggers if 'go' is pressed in the machine learning tab
  observeEvent(input$do_ml, {
    withProgress({
      
      setProgress(value = 0)
      # prepare matrix
      
	  # get base table to use for process
      curr <- as.data.table(mSet$dataSet$preproc) # the filtered BUT NOT IMPUTED table, ML should be able to deal w/ missing values
      
	  # replace NA's with zero
      curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
      
	  # conv to data frame
      curr <- as.data.frame(curr)
      rownames(curr) <- rownames(mSet$dataSet$preproc)
      
	  # find the qc rows
      is.qc <- grepl("QC|qc", rownames(curr))
      curr <- curr[!is.qc,]
      
	  # reorder according to covars table (will be used soon)
      order <- match(mSet$dataSet$covars$sample,rownames(curr))

	  # get covariates to add to the metabolite table
      config <- mSet$dataSet$covars[order, -"label"] # reorder so both halves match up later
      config <- cbind(config, label=mSet$dataSet$cls[order]) # add current experimental condition
      config <- config[,apply(!is.na(config), 2, any), with=FALSE]
      
      # remove ones w/ every row being different(may be used to identify...)
      covariates <- lapply(1:ncol(config), function(i) as.factor(config[,..i][[1]]))
      names(covariates) <- colnames(config)
      
	  # remove ones with na present
      has.na <- sapply(covariates, function(x) any(is.na(x)))
      has.all.unique <- sapply(covariates, function(x) length(unique(x)) == length(x))
      
	  # rename the variable of interest to 0-1-2 etc.
      char.lbl <- as.character(covariates$label)
      uniques <- unique(char.lbl)
      uniques_new_name <- c(1:length(uniques))
      names(uniques_new_name) = uniques
      
      remapped.lbl <- uniques_new_name[char.lbl]
      
	  # find which variables are covariant with label, they will be removed
      covariant.with.label <- sapply(covariates, function(x){
        char.x <- as.character(x)
        uniques <- unique(char.x)
        uniques_new_name <- c(1:length(uniques))
        names(uniques_new_name) = uniques
        remapped.x = uniques_new_name[char.x]
        res = if(length(remapped.x) == length(remapped.lbl)){
          all(remapped.x == remapped.lbl)
        }else{ FALSE }
        res
      })
      
	  # now filter out unique, covariate and with missing columns from $covars
      keep_configs <- which(!(names(config) %in% names(config)[unique(c(which(has.na), which(has.all.unique), which(covariant.with.label)))]))
      keep_configs <- c(keep_configs, which(names(config) == "label"))
      
      #print("Removing covariates and unique columns. Keeping non-mz variables:")
      #print(names(config)[keep_configs])
      
      config <- config[,..keep_configs,with=F]
      # - - - - - - - - - - - - - - - - - - - - - - -
      
	  # join halves together, user variables and metabolite data
      curr <- cbind(config, curr)
      curr <- as.data.table(curr)

      # remove cols with all NA
      curr <- curr[,colSums(is.na(curr))<nrow(curr), with=FALSE]
      
	  # identify which columns are metabolites and which are config/covars
      configCols <- which(!(colnames(curr) %in% colnames(mSet$dataSet$norm)))
      mzCols <- which(colnames(curr) %in% colnames(mSet$dataSet$norm))
      
	  # make the covars factors and the metabolites numeric.
      curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
      curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]
      
	  # how many models will be built? user input
      goes = as.numeric(input$ml_attempts)
      
	  # get results for the amount of attempts chosen
      repeats <- pbapply::pblapply(1:goes, function(i){
        
        # get regex user input for filtering testing and training set
        ml_train_regex <<- input$ml_train_regex
        ml_test_regex <<- input$ml_test_regex
        
		# get user training percentage
        ml_train_perc <- input$ml_train_perc/100
        
        if(ml_train_regex == "" & ml_test_regex == ""){ # BOTH ARE NOT DEFINED
          test_idx = caret::createDataPartition(y = curr$label, p = ml_train_perc, list = FALSE) # partition data in a balanced way (uses labels)
          train_idx = setdiff(1:nrow(curr), test_idx) #use the other rows for testing
          inTrain = train_idx
          inTest = test_idx
        }else if(ml_train_regex != ""){ #ONLY TRAIN IS DEFINED
          train_idx = grep(config$sample, pattern = ml_train_regex) # get training sample ids with regex
          test_idx = setdiff(1:nrow(curr), train_idx) # use the other rows for testing
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
          inTrain <- train_idx[reTrain$Resample1]
          inTest = test_idx
        }else{ # ONLY TEST IS DEFINED
          test_idx = grep(config$sample, pattern = ml_test_regex) # get training sample ids with regex
          train_idx = setdiff(1:nrow(curr), test_idx) # use the other rows for testing
          reTrain <- caret::createDataPartition(y = config[train_idx, label], p = ml_train_perc) # take a user-defined percentage of the regexed training set
          inTrain <- train_idx[reTrain$Resample1] 
          inTest <- test_idx
        }
        
		# choose predictor "label" (some others are also included but cross validation will be done on this)
        predictor = "label"
        
		# split training and testing data
        trainY <- curr[inTrain, 
                       ..predictor][[1]]
        testY <- curr[inTest,
                      ..predictor][[1]]
        
		# remove predictive column from training set 
        remove.cols <- c("label") #TODO: make group column removed at start
        remove.idx <- which(colnames(curr) %in% remove.cols)
        training <- curr[inTrain,-remove.idx, with=FALSE]
        testing <- curr[inTest,-remove.idx, with=FALSE]
        
		# get covar indices
        predIdx <- which(colnames(curr) %in% colnames(config))
       
        # remove unused levels in the factor part of the table after filtering
        training <- data.matrix(gdata::drop.levels(training))
        testing <- data.matrix(gdata::drop.levels(testing))
        
        setProgress(value = i/goes)
        
		# train and cross validate model
        switch(input$ml_method,
               rf = { # random forest
                 model = randomForest::randomForest(x = training, 
                                                    y = trainY, 
                                                    ntree = 500, # amount of trees made TODO: make user choice
                                                    importance=TRUE) # include variable importance in model
                 
                 prediction <- stats::predict(model, 
                                              testing, 
                                              "prob")[,2]
                 # get importance table
                 importance = as.data.table(model$importance, keep.rownames = T)
                 rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
                 rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)] # reorder for convenience
                 rf_tab <- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
				 # return list with model, prediction on test data etc.
                 list(type="rf",
                      feats = as.data.table(rf_tab, keep.rownames = T), 
                      model = model,
                      prediction = prediction,
                      labels = testY)
               }, 
               ls = { # lasso
				 # user x cross validation input
                 nfold = switch(input$ml_folds, 
                                "5" = 5,
                                "10" = 10,
                                "20" = 20,
                                "50" = 50,
                                "LOOCV" = length(trainY)) # leave one out CV
								
                 family = "binomial" #TODO: enable multinomial for multivariate data!!
                 
				 # 
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
    enrich_overview_tab <<- pw_memb_tab
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
  
  observe({
      datamanager$reload <- input$statistics
  })
  
  observeEvent(input$pca_2d3d, {
    datamanager$reload <- "pca"
  },ignoreInit = TRUE, ignoreNULL = T)
  
  observeEvent(input$plsda_2d3d, {
    datamanager$reload <- "plsda"
  },ignoreInit = TRUE, ignoreNULL = T)

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
                         mSet$analSet$ml$ls[[name]]$tophits[order(mSet$analSet$ml$ls[[name]]$tophits$count, decreasing = T),]$mz})
                       names(tbls_ls) <- input$ls_choice
                       # - - - 
                       tbls_ls
                     },
                     rf = {
                       tbls_rf <- lapply(input$rf_choice, function(name){
                         mSet$analSet$ml$rf[[name]]$tophits[order(mSet$analSet$ml$rf[[name]]$tophits$mda, decreasing = T),]$mz})
                       names(tbls_rf) <- input$rf_choice
                       # - - -
                       tbls_rf
                     },
                     plsda = {
                       tbls_plsda <- lapply(input$plsda_choice, function(name){
                         compounds_pc <- as.data.table(mSet$analSet$plsda$vip.mat,keep.rownames = T)
                         colnames(compounds_pc) <- c("rn", paste0("PC", 1:(ncol(compounds_pc)-1)))
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
