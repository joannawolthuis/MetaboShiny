
# observeEvent(input$metadata, {
#   if("files" %in% names(input$metadata)){
#     meta_path <- parseFilePaths(gbl$paths$volumes, input$metadata)$datapath
#     type = tools::file_ext(meta_path)
#     choices = switch(type,
#                      csv = colnames(fread(meta_path)),
#                      xlsx = c(colnames(openxlsx::read.xlsx(meta_path, sheet="Setup"), openxlsx::read.xlsx(meta_path,sheet="Individual Data")))
#                      )
#     updateSelectInput(session, inputId = "meta_dupli_col", choices = choices)
#   }
# })

# create checkcmarks if database is present
lapply(c("merge", "db", "csv"), FUN=function(col){
  # creates listener for if the 'check db' button is pressed
  observe({
    # see which db files are present in folder
    folder_files <- list.files(lcl$paths$work_dir)
    is.present <- switch(col,
                         merge = {
                           if(!is.null(input$importmode)){
                             switch(input$importmode,
                                    db = {
                                      is.list(input$database) & is.list(input$metadata) 
                                    },
                                    csv = {
                                      is.list(input$metadata) & is.list(input$outlist_pos) & is.list(input$outlist_neg)
                                    })  
                           }else{
                             F
                           }
                         },
                         db = paste0(input$proj_name_new, ".db") %in% folder_files,
                         csv = paste0(input$proj_name_new, ".csv") %in% folder_files)
    check_pic <- if(is.present) "yes.png" else "no.png"
    # generate checkmark image objects
    output[[paste0("proj_", col, "_check")]] <- renderImage({
      filename <- normalizePath(file.path('www', check_pic))
      list(src = filename, width = 70,
           height = 70)
    }, deleteFile = FALSE)# <- this is important or the checkmark file is deleted, haha
  })
})

# triggers when user wants to create database from .db and excel or 2 csv files and excel
observeEvent(input$create_db,{

  files.present = switch(input$importmode,
         db = {
           is.list(input$database) & is.list(input$metadata) 
         },
         csv = {
           is.list(input$metadata) & is.list(input$outlist_pos) & is.list(input$outlist_neg)
         })
  
  if(!files.present) return(NULL)
  
  # update the path to patient db
  lcl$paths$patdb <<- file.path(lcl$paths$work_dir, paste0(lcl$proj_name, ".db"))

  withProgress({

    shiny::setProgress(session=session, value= .1)

    proj_name = input$proj_name_new

    updateSelectizeInput(session = session,
                         inputId = "proj_name",
                         choices = c(lcl$vectors$project_names, proj_name))

    updateSelectizeInput(session = session,
                         inputId = "proj_name",
                         selected = proj_name)

    lcl$proj_name <<- proj_name
    lcl$paths$patdb <<- file.path(lcl$paths$work_dir, paste0(lcl$proj_name, ".db"))
    # change project name in user options file
    setOption(lcl$paths$opt.loc, key="proj_name", value=lcl$proj_name)
    # print the changed name in the UI
    output$proj_name <<- renderText(proj_name)
    # change path CSV should be / is saved to in session
    #lcl$paths$csv_loc <<- file.path(lcl$paths$work_dir, paste0(lcl$proj_name,".csv"))

    switch(input$importmode,
           # if loading in a .db file... (FAST, MOSTLY FOR ADMINS USING HPC)
           db = {

             # get the db and excel path from the UI elemnts
             db_path <- parseFilePaths(gbl$paths$volumes, input$database)$datapath
             excel_path <- parseFilePaths(gbl$paths$volumes, input$metadata)$datapath

             # copy the user selected db to the processing folder under proj_name renaming
             file.copy(db_path, lcl$paths$patdb, overwrite = T)
             
             shiny::setProgress(session=session, value= .30)

             # add metadata file to .db file generated in previous step
             metadata_path <- parseFilePaths(gbl$paths$volumes, input$metadata)$datapath

             if(grepl(metadata_path, pattern = "csv")){
               exp_vars <<- load.metadata.csv(metadata_path, lcl$paths$patdb, ppm=input$ppm)
             }else{
               exp_vars <<- load.metadata.excel(metadata_path, lcl$paths$patdb, ppm=input$ppm)
             }

             shiny::setProgress(session=session, value= .60)

           },
           # if loading in .csv files...
           csv = {

             # build patient db from csv files with a given ppm error margin
             build.pat.db(lcl$paths$patdb,
                          ppm = input$ppm,
                          pospath = parseFilePaths(gbl$paths$volumes, input$outlist_pos)$datapath,
                          negpath = parseFilePaths(gbl$paths$volumes, input$outlist_neg)$datapath,
                          overwrite = T)

             shiny::setProgress(session=session, value= .95,message = "Adding metadata to database...")

             # add excel file to .db file generated in previous step
             metadata_path <- parseFilePaths(gbl$paths$volumes, input$metadata)$datapath

             if(grepl(metadata_path, pattern = "csv")){
               exp_vars <<- load.metadata.csv(metadata_path, lcl$paths$patdb)
             }else{
               exp_vars <<- load.metadata.excel(metadata_path, lcl$paths$patdb)
             }
          }
    )
    
    output$proj_db_check <- renderImage({
      filename <- normalizePath(file.path('www', "yes.png"))
      list(src = filename, width = 70,
           height = 70)
      },deleteFile = FALSE)
    })
})

# imports existing db file
# TODO: is deprecated, fix!!
observeEvent(input$import_db, {

  lcl$paths$patdb <<- input$pat_db$datapath
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
  lcl$paths$csv_loc <<- input$pat_csv$datapath

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
    tbl <- get.csv(lcl$paths$patdb,
                   group_adducts = F, # if(length(lcl$vectors$db_add_list) == 0) F else T, # group by addicts?
                   groupfac = "mz" #input$group_by, # group by mz or formula
                   #which_dbs = lcl$vectors$db_add_list, # used databases
                   #which_adducts = selected_adduct_list # used adducts
    )

    # check if any samples are duplicated in the resulting table
    # IF SO: rename to time series data format because that is the only one that should have duplicate sample names
    if(any(duplicated(tbl$sample))){
      tbl$sample <- paste0(tbl$sample, 
                           "_T", 
                           tbl$time)
      show.times = T # make sure the 'times' column shows in the csv making result table
    }else{
      show.times = F
    }

    shiny::setProgress(session=session, value= 2/4)

    # change location of csv in global
    lcl$paths$csv_loc <<- file.path(lcl$paths$work_dir, paste0(lcl$proj_name,".csv"))

    # write csv file to new location
    fwrite(tbl, lcl$paths$csv_loc, sep="\t")

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
    output$proj_csv_check <- renderImage({
      filename <- normalizePath(file.path('www', "yes.png"))
      list(src = filename, width = 70,
           height = 70)
    },deleteFile = FALSE)
  })
})

# triggers when 'get options' is clicked in the normalization pane
observeEvent(input$check_csv, {
  # ----------------------

  # read in csv
  # TODO: only read in the first x -rows- to save time??
  csv <- fread(lcl$paths$csv_loc,
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
