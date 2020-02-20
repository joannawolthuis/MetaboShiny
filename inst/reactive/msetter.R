shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    
    if(!is.null(mSet)){
      
      mSet <- MetaboShiny::store.mSet(mSet) # save analyses
      mSet <- MetaboShiny::reset.mSet(mSet) # reset dataset
      
      success = F
      
      try({
        if(!(mSetter$do %in% c("unsubset"))){
          mSet.settings <- if(mSetter$do == "load") mSet$storage[[input$storage_choice]]$settings else mSet$settings
          if(length(mSet.settings$subset) > 0){
            subs = mSet.settings$subset
            subs = subs[names(subs) != "sample"]
            if(length(subs) > 0){
              for(i in 1:length(subs)){
                mSet <- MetaboShiny::subset.mSet(mSet, 
                                                 subset_var = names(subs)[i], 
                                                 subset_group = subs[[i]])  
              }  
            }
          }
        }else{
          mSet.settings <- mSet$settings
        }
        
        # TODO: fix mismatch between cls order, covars order and table sample order
        # use match() somehow?
        # should be fixed in #CHANGE.MSET
        mSet <- switch(mSetter$do,
                       load = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.check$settings$exp.var, 
                                                          time_var =  mSet.check$settings$time.var,
                                                          stats_type = mSet.check$settings$exp.type)
                         mSet$dataSet$paired <- mSet.settings$paired
                         mSet
                       },
                       change = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = input$stats_var, 
                                                          stats_type = input$stats_type, 
                                                          time_var = input$time_var)
                         mSet$dataSet$paired <- if(input$stats_type %in% c("t", "t1f") | input$paired) TRUE else FALSE
                         if(input$omit_unknown & grepl("^1f", input$stats_type)){
                           shiny::showNotification("omitting 'unknown' labeled samples...")
                           knowns = mSet$dataSet$covars$sample[which(mSet$dataSet$covars[ , input$stats_var, with=F][[1]] != "unknown")]
                           if(length(knowns) > 0){
                             mSet <- MetaboShiny::subset.mSet(mSet,
                                                              subset_var = "sample", 
                                                              subset_group = knowns) 
                           }
                         }
                         mSet
                       },
                       subset = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.check$settings$exp.var, 
                                                          time_var =  mSet.check$settings$time.var,
                                                          stats_type = mSet.check$settings$exp.type)
                         mSet <- MetaboShiny::subset.mSet(mSet,
                                                          subset_var = input$subset_var, 
                                                          subset_group = input$subset_group)
                         mSet$dataSet$paired <- mSet.settings$paired
                         mSet
                       },
                       unsubset = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.check$settings$exp.var, 
                                                          time_var =  mSet.check$settings$time.var,
                                                          stats_type = mSet.check$settings$exp.type)
                         mSet$dataSet$paired <- mSet.settings$paired
                         mSet
                       }
        ) 
        
        if(mSet$dataSet$paired){
          mSet$settings$paired <- TRUE
          mSet <- MetaboShiny::pair.mSet(mSet)
        }
        if(grepl(mSet$settings$exp.type, pattern = "^1f")){
          if(mSet$dataSet$cls.num == 2){
            mSet$dataSet$exp.type <- "1fb"
          }else{
            mSet$dataSet$exp.type <- "1fm"
          }  
        }
        
        mSet$settings$exp.type = mSet$dataSet$exp.type 
        new.name = if(mSetter$do == "load") input$storage_choice else MetaboShiny::name.mSet(mSet)
        
        if(new.name %in% names(mSet$storage)){
          mSet <- MetaboShiny::load.mSet(mSet, new.name)
        }else{
          mSet$analSet <- list(type = "stat")
        }
        
        mSet$dataSet$cls.name <- new.name

        shiny::updateCheckboxInput(session, 
                                   "paired", 
                                   value = mSet$dataSet$paired) 
        lcl$last_mset <- mSet$dataSet$cls.name
        success = T
        shiny::showNotification("Changed mSet...")
      })
      
      if(success){
        if(MetaboShiny::is.ordered.mSet(mSet)){
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
      }else{
        MetaboShiny::metshiAlert("Failed! Restoring old mSet...")
      }
      
      mSetter$do <- NULL
    
      }
  }
})