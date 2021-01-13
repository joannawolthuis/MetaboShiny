shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    
    if(!is.null(mSet)){
      
      try({
        mSet <- store.mSet(mSet) # save analyses
        success = F
        
        if(mSetter$do == "load" & !is.null(mSet$storage[[input$storage_choice]]$data)){
           mSet <- load.mSet(mSet, input$storage_choice)
        }else{
          
          oldSettings <- mSet$settings
          
          mSet <- reset.mSet(mSet,
                             fn = file.path(lcl$paths$proj_dir, 
                                            paste0(lcl$proj_name,
                                                   "_ORIG.metshi")))
          
          orig.count <- mSet$metshiParams$orig.count
          
          if(!(mSetter$do %in% c("unsubset"))){
            mSet.settings <- if(mSetter$do == "load") mSet$storage[[input$storage_choice]]$settings else oldSettings
            if(length(mSet.settings$subset) > 0){
              subs = mSet.settings$subset
              subs = subs[!(names(subs) %in% c("sample", "mz"))]
              if(length(subs) > 0){
                for(i in 1:length(subs)){
                  mSet <- subset_mSet(mSet, 
                                      subset_var = names(subs)[i], 
                                      subset_group = subs[[i]])  
                }  
              }
            }
          }else{
            mSet.settings <- oldSettings
          }
          
          mSet$settings <- mSet.settings
          
          mSet <- switch(mSetter$do,
                         refresh = {
                           mSet$dataSet$ispaired <- mSet.settings$ispaired
                           mSet
                         },
                         load = {
                           mSet$dataSet$ispaired <- mSet.settings$ispaired
                           mSet
                         },
                         change = {
                           mSet$dataSet$ispaired <- if(input$stats_type %in% c("t", "t1f") | input$paired) TRUE else FALSE
                           mSet
                         },
                         subset = {
                           mSet <- subset_mSet(mSet,
                                               subset_var = input$subset_var, 
                                               subset_group = input$subset_group)
                           mSet$dataSet$ispaired <- mSet.settings$ispaired
                           mSet
                         },
                         unsubset = {
                           mSet$dataSet$ispaired <- mSet.settings$ispaired
                           mSet$settings$subset <- list()
                           mSet
                         }) 
          
          mSet$analSet <- list(type = "stat")
          mSet$analSet$type <- "stat"
          
          if(input$redo_upon_change){
            mSet$dataSet$orig <- mSet$dataSet$start
            mSet$dataSet$start <- mSet$dataSet$preproc <- mSet$dataSet$proc <- mSet$dataSet$prenorm <- NULL
            mSet <- metshiProcess(mSet)
          }
          
          if(mSetter$do == "change"){
            if(input$omit_unknown & grepl("^1f", input$stats_type)){
              shiny::showNotification("omitting 'unknown' labeled samples...")
              knowns = mSet$dataSet$covars$sample[which(mSet$dataSet$covars[ , input$stats_var, with=F][[1]] != "unknown")]
              if(length(knowns) > 0){
                mSet <- subset_mSet(mSet,
                                    subset_var = "sample", 
                                    subset_group = knowns) 
              }
            }
            mSet <- change.mSet(mSet, 
                                stats_var = input$stats_var, 
                                stats_type = input$stats_type, 
                                time_var = input$time_var)
          }else{
            if(input$omit_unknown & grepl("^1f", mSet$settings$exp.type)){
              shiny::showNotification("omitting 'unknown' labeled samples...")
              knowns = mSet$dataSet$covars$sample[which(mSet$dataSet$covars[ , mSet$settings$exp.var, with=F][[1]] != "unknown")]
              if(length(knowns) > 0){
                mSet <- subset_mSet(mSet,
                                    subset_var = "sample", 
                                    subset_group = knowns) 
              }
            }
            mSet <- change.mSet(mSet, 
                                stats_var = mSet.settings$exp.var, 
                                time_var =  mSet.settings$time.var,
                                stats_type = mSet.settings$exp.type)
          }
          
          new.name = if(mSetter$do == "load") input$storage_choice else name.mSet(mSet)
          
          if(new.name %in% names(mSet$storage)){
            mSet <- load.mSet(mSet, new.name)
          }
          
          mSet$settings$cls.name <- new.name
          
          if(mSet$dataSet$ispaired){
            mSet$settings$ispaired <- TRUE
            mSet <- pair.mSet(mSet)
          }else{
            mSet.settings$ispaired <- FALSE
          }
          if(grepl(mSet$settings$exp.type, pattern = "^1f")){
            if(mSet$dataSet$cls.num == 2){
              mSet$settings$exp.type <- "1fb"
            }else{
              mSet$settings$exp.type <- "1fm"
            }  
          }
        }   
        success=T
      })
      
      if(success){
        if(is.ordered.mSet(mSet)){
          msg = "mSet class label order still correct! :)"
          try({
            shiny::showNotification(msg) 
          })
          print(msg)
        }else{
          msg = "ordering went wrong with the mset sample order! :("
          try({
            shiny::showNotification(msg)
          })
          print(msg)
        }
        mSet <<- mSet
        lcl$has_changed <<- TRUE
        uimanager$refresh <- c("general", "ml")
      }else{
        metshiAlert("Failed! Restoring old mSet...")
      }
      mSetter$do <- NULL
    }
  }
})