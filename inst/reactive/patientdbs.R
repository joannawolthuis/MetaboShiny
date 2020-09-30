# create checkcmarks if database is present
lapply(c("merge", 
         "db", 
         "csv"), FUN=function(col){
  # creates listener for if the 'check db' button is pressed
  shiny::observe({
    # see which db files are present in folder
    folder_files <- list.files(lcl$paths$proj_dir)
    is.present <- switch(col,
                         merge = {
                           if(!is.null(input$ms_modes)){
                             all(sapply(input$ms_modes, function(ionMode){
                               is.list(input[[paste0("outlist_", ionMode)]])
                             }) & is.list(input$metadata))  
                           }else FALSE
                         },
                         csv = paste0(input$proj_name_new, ".csv") %in% folder_files,
                         db = paste0(input$proj_name_new, ".db") %in% folder_files)
    check_pic <- if(is.present) "yes.png" else "no.png"
    # generate checkmark image objects
    output[[paste0("proj_", col, "_check")]] <- shiny::renderImage({
      filename <- normalizePath(file.path('www', check_pic))
      list(src = filename, height = "70px")
    }, deleteFile = FALSE)# <- this is important or the checkmark file is deleted, haha
  })
})

shiny::observeEvent(input$ms_modes, {
  pickerUI = lapply(input$ms_modes, function(mode){
    char=if(mode == "pos") "+" else "-"
    shinyFiles::shinyFilesButton(paste0('outlist_',mode), 
                                 gsubfn::fn$paste('upload $char peaks'), 
                                 gsubfn::fn$paste('Upload $char mode peaks'), 
                                 FALSE)
  })
  output$outlist_pickers <-shiny::renderUI(pickerUI)
})

# create reactive object
# display should be: "Mz left: 100/120000
missValues <- shiny::reactiveValues(pos = c(),
                                    neg = c())

output$missMzRatio <- shiny::renderUI({
  if(length(missValues$pos) > 0){
    posTotal = length(missValues$pos)
    posKeep = length(which(missValues$pos < input$perc_limit_mz))
    posText = paste0("Remaining (+): ", posKeep, " / ", posTotal)
  }else{
    posText = NULL
  }
  if(length(missValues$neg) > 0){
    negTotal = length(missValues$neg)
    negKeep = length(which(missValues$neg < input$perc_limit_mz))
    negText = paste0("Remaining (-):", negKeep, " / ", negTotal)
  }else{
    negText = NULL
  }
  list(if(!is.null(posText)) shiny::helpText(posText),
       if(!is.null(posText)) shiny::helpText(negText))
})

output$missMzPlot <- shiny::renderPlot({
  if(length(missValues$pos) > 0 | length(missValues$neg > 0)){
    # plot(density(missValues$pos),col="blue",main="Missing value distribution",xlab = "% missing", lwd=2)
    # lines(density(missValues$neg),col="red",lwd=2)
    # abline(v = input$perc_limit_mz, col="black", lwd=1.5, lty=2)
    # text(x = input$perc_limit_mz, y = 0.5, "Current threshold",pos = 4,cex = 1.5)  
    data = data.table::data.table(variable = c(rep("pos", length(missValues$pos)),
                                               rep("neg", length(missValues$neg))),
                                  "Missing percentage" = c(missValues$pos, 
                                                           missValues$neg))
    ggplot2::ggplot(data = data) + 
      ggplot2::geom_density(aes(x=`Missing percentage`, fill=variable, y=..scaled..), alpha=0.5) +
      ggplot2::geom_vline(xintercept=input$perc_limit_mz, size=0.7, linetype="dotted") +
      # ggplot2::geom_label(y=0.5,
      #                     x=input$perc_limit_mz,
      #                     label="Current threshold")+
      ggplot2::annotate(geom = "label", y=0.5, x = input$perc_limit_mz, label = "Current threshold") +
      gbl$functions$plot.themes[[lcl$aes$theme]](base_size = 15) + 
      ggplot2::theme(legend.position="none",
                     axis.line = ggplot2::element_line(colour = 'black', size = .5),
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        vjust = 0.1,
                                                        size=lcl$aes$font$title.size*1.2),
                     text = ggplot2::element_text(family = lcl$aes$font$family))
      
  }
})

# triggers when user wants to check missing values profile
shiny::observeEvent(input$checkMiss, {
  files.present = all(sapply(input$ms_modes, function(ionMode){
    is.list(input[[paste0("outlist_", ionMode)]])
  }))
  
  file.pos = shinyFiles::parseFilePaths(gbl$paths$volumes, input$outlist_pos)$datapath
  file.neg = shinyFiles::parseFilePaths(gbl$paths$volumes, input$outlist_neg)$datapath
  
  if(files.present){
    if(is.list(input$outlist_pos)){
      nrows = length(vroom::vroom_lines(file.pos, altrep = TRUE, progress = TRUE)) - 1L
      missPos = getMissing(file.pos, nrow=nrows)
      missValues$pos = missPos/nrows * 100
    }
    if(is.list(input$outlist_neg)){
      nrows = length(vroom::vroom_lines(file.neg, altrep = TRUE, progress = TRUE)) - 1L
      missNeg = getMissing(file.neg, nrow=nrows)
      missValues$neg = missNeg/length(missPos) * 100
    }
  }else{
    MetaboShiny::metshiAlert("Please select files first!")
  }
})


# triggers when user wants to create database from .db and excel or 2 csv files and excel
shiny::observeEvent(input$create_csv,{

  files.present = all(sapply(input$ms_modes, function(ionMode){
    is.list(input[[paste0("outlist_", ionMode)]])
  }) & is.list(input$metadata))
  
  if(!files.present) return(NULL)
  
  # update the path to patient csv
  shiny::withProgress({

    success=F
    try({
      shiny::setProgress(session=session, value= .1)
      
      proj_name = input$proj_name_new
      
      shiny::updateSelectizeInput(session = session,
                                  inputId = "proj_name",
                                  choices = c(lcl$vectors$project_names, proj_name))
      
      shiny::updateSelectizeInput(session = session,
                                  inputId = "proj_name",
                                  selected = proj_name)
      
      lcl$proj_name <<- proj_name
      lcl$paths$proj_dir <<- file.path(lcl$paths$work_dir, proj_name)
      
      if(dir.exists(lcl$paths$proj_dir)) unlink(lcl$paths$proj_dir)
      dir.create(lcl$paths$proj_dir,showWarnings = F)
      
      lcl$paths$csv_loc <<- file.path(lcl$paths$proj_dir, paste0(proj_name, ".csv"))
      lcl$paths$patdb <<- file.path(lcl$paths$proj_dir, paste0(proj_name, ".db"))
      
      # change project name in user options file
      MetaboShiny::setOption(lcl$paths$opt.loc, key="proj_name", value=proj_name)
      # print the changed name in the UI
      output$proj_name <<- shiny::renderText(proj_name)
      
      hasPos = "pos" %in% input$ms_modes
      hasNeg = "neg" %in% input$ms_modes
      
      # if loading in .csv files...
      import.pat.csvs(ppm = input$ppm,
                      pospath = if(hasPos) shinyFiles::parseFilePaths(gbl$paths$volumes, input$outlist_pos)$datapath else c(),
                      negpath = if(hasNeg) shinyFiles::parseFilePaths(gbl$paths$volumes, input$outlist_neg)$datapath else c(),
                      metapath = shinyFiles::parseFilePaths(gbl$paths$volumes, input$metadata)$datapath,
                      wipe.regex = input$wipe_regex,
                      missperc.mz = input$perc_limit_mz,
                      missperc.samp = input$perc_limit_samp,
                      missList = missValues,
                      csvpath = lcl$paths$csv_loc,
                      overwrite = T,
                      inshiny=F,
                      roundMz=input$roundMZ)
      
      success=T
      output$proj_csv_check <- shiny::renderImage({
        filename <- normalizePath(file.path('www', "yes.png"))
        list(src = filename, width = 70,
             height = 70)
      },deleteFile = FALSE)
    })
    if(!success){
      MetaboShiny::metshiAlert("Something is wrong with your input data!")
      gc()
    }
    })
})

# imports existing db file
# TODO: is deprecated, fix!!
shiny::observeEvent(input$import_db, {
  lcl$paths$patdb <<- input$pat_db$datapath
  output$db_upload_check <- shiny::renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath('www/yes.png')
    # Return a list containing the filename and alt text
    list(src = filename, width = 20,
         height = 20)
  }, deleteFile = FALSE)
})


shiny::observeEvent(input$metadata_new_add, {
  
  meta_path <- shinyFiles::parseFilePaths(gbl$paths$volumes, input$metadata_new)$datapath
  success = F
  try({
    new_meta <- data.table::fread(meta_path)
    new_meta <- MetaboShiny::reformat.metadata(new_meta)
    colnames(new_meta) <- tolower(colnames(new_meta))
    mSet <- MetaboShiny::store.mSet(mSet)
    mSet <- MetaboShiny::reset.mSet(mSet,
                                    fn = file.path(lcl$paths$proj_dir, 
                                                   paste0(lcl$proj_name,
                                                          "_ORIG.metshi")))
    mSet$dataSet$covars <- plyr::join(mSet$dataSet$covars[,"sample"], new_meta, type = "left")
    for(project in names(mSet$storage)){
      if("dataSet" %in% names(mSet$storage[[project]])){
        mSet$storage[[project]]$dataSet$covars <- plyr::join(mSet$storage[[project]]$dataSet$covars[,"sample"], 
                                                             new_meta,
                                                             type = "left")
      }
    }
    success = T
  })
  if(success){
    mSet <<- mSet
    # overwrite orig?
    save(mSet, file = file.path(lcl$paths$proj_dir, 
                                paste0(lcl$proj_name,"_ORIG.metshi")))
    mSet$dataSet$missing <- mSet$dataSet$orig <- mSet$dataSet$start <- NULL 
    fn <- paste0(tools::file_path_sans_ext(lcl$paths$csv_loc), ".metshi")
    save(mSet, file = fn)
    shiny::showNotification("Updated metadata!")
    uimanager$refresh <- "general"
  }else{
    shiny::showNotification("Something went wrong! :(")
  }
})

# triggers when 'get options' is clicked in the normalization pane
shiny::observeEvent(input$check_csv, {
  # ----------------------

  # get the names of those experimental variables
  header = data.table::fread(lcl$paths$csv_loc, nrows = 2, header=T)
  met.cols = getColDistribution(header)$meta
  opts = colnames(header)[met.cols]
  metadata = data.table::fread(lcl$paths$csv_loc, header=T, select=colnames(header)[met.cols])
  
  # get columns that can be used for batch correction (need to be non-unique)
  batch <<- which(sapply(opts, function(x) length(unique(metadata[,..x][[1]])) < nrow(metadata)))

  # update the possible options in the UI
  shiny::updateSelectInput(session, "samp_var",
                    choices = opts)
  shiny::updateSelectizeInput(session, "batch_var",
                       choices = opts[batch],
                       options = list(maxItems = 3L - (length(input$batch_var)))
  )
})
