# triggers when user wants to create database from .db and excel or 2 csv files and excel
observeEvent(input$create_db,{
  
  # update the path to patient db
  global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name, ".db"))
  
  withProgress({
    
    shiny::setProgress(session=session, value= .1)
    
    switch(input$new_proj,
           # if loading in a .db file... (FAST, MOSTLY FOR ADMINS USING HPC)
           `From DB` = {
             
             # get the db and excel path from the UI elemnts
             db_path <- parseFilePaths(global$paths$volumes, input$database)$datapath
             excel_path <- parseFilePaths(global$paths$volumes, input$excel)$datapath
             
             # copy the user selected db to the processing folder under proj_name renaming
             file.copy(db_path, global$paths$patdb, overwrite = T)
             
             shiny::setProgress(session=session, value= .30)
             
             # add metadata file to .db file generated in previous step
             metadata_path <- parseFilePaths(global$paths$volumes, input$metadata)$datapath

             if(grepl(metadata_path, pattern = "csv")){
               exp_vars <<- load.metadata.csv(metadata_path, global$paths$patdb)
             }else{
               exp_vars <<- load.metadata.excel(metadata_path, global$paths$patdb)
             }
             
             shiny::setProgress(session=session, value= .60)
             
           },
           # if loading in .csv files...
           `From CSV` = {
             
             proj_name = input$proj_name_new
             updateSelectizeInput(session = session, 
                                  inputId = "proj_name", 
                                  choices = c(global$vectors$project_names, proj_name))
             
             updateSelectizeInput(session = session,
                                  inputId = "proj_name", 
                                  selected = proj_name)
             
             global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(proj_name,".db", sep=""))
             # change project name in user options file
             setOption("user_options.txt", "proj_name", proj_name)
             # print the changed name in the UI
             output$proj_name <<- renderText(proj_name)
             # change path CSV should be / is saved to in session
             global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
             
             # build patient db from csv files with a given ppm error margin
             build.pat.db(global$paths$patdb,
                          ppm = input$ppm,
                          pospath = parseFilePaths(global$paths$volumes, input$outlist_pos)$datapath,
                          negpath = parseFilePaths(global$paths$volumes, input$outlist_neg)$datapath,
                          overwrite = T)
             
             shiny::setProgress(session=session, value= .95,message = "Adding metadata to database...")
             
             # add excel file to .db file generated in previous step
             metadata_path <- parseFilePaths(global$paths$volumes, input$metadata)$datapath

             if(grepl(metadata_path, pattern = "csv")){
               exp_vars <<- load.metadata.csv(metadata_path, global$paths$patdb)
             }else{
               exp_vars <<- load.metadata.excel(metadata_path, global$paths$patdb)
             }
          }
    )
  })
})

# imports existing db file
# TODO: is deprecated, fix!!
observeEvent(input$import_db, {
  
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
                   group_adducts = if(length(global$vectors$db_add_list) == 0) F else T, # group by addicts?
                   groupfac = input$group_by, # group by mz or formula
                   which_dbs = global$vectors$db_add_list, # used databases
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