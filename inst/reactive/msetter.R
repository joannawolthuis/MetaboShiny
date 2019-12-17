shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    
    if(!is.null(mSet)){
      
      print(mSetter$do)
      
      mSet <- MetaboShiny::store.mSet(mSet) # save analyses
      
      mSet.old <- mSet
      
      mSet <- MetaboShiny::reset.mSet(mSet) # reset dataset
      
      success = F
      
      try({
        if(mSetter$do != "unsubset"){
          if(length(mSet.old$dataSet$subset) > 0){
            subs = mSet.old$dataSet$subset
            subs = subs[names(subs) != "sample"]
            if(length(subs) > 0){
              for(i in 1:length(subs)){
                mSet <- MetaboShiny::subset.mSet(mSet, 
                                                 subset_var = names(subs)[i], 
                                                 subset_group = subs[[i]])  
              }  
            }
          }
        }
        
        # TODO: fix mismatch between cls order, covars order and table sample order
        # use match() somehow?
        # should be fixed in #CHANGE.MSET
        
        mSet <- switch(mSetter$do,
                       change = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = input$stats_var, 
                                                          stats_type = input$stats_type, 
                                                          time_var = input$time_var)
                         mSet$dataSet$paired <- if(input$stats_type %in% c("t", "t1f") | input$paired) TRUE else FALSE
                         mSet
                       },
                       subset = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.old$dataSet$exp.var, 
                                                          time_var =  mSet.old$dataSet$time.var,
                                                          stats_type = mSet.old$dataSet$exp.type)
                         mSeta <- MetaboShiny::subset.mSet(mSet,
                                                          subset_var = input$subset_var, 
                                                          subset_group = input$subset_group)
                         mSet$dataSet$paired <- mSet.old$dataSet$paired
                         mSet
                       },
                       unsubset = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.old$dataSet$exp.var, 
                                                          time_var =  mSet.old$dataSet$time.var,
                                                          stats_type = mSet.old$dataSet$exp.type)
                         mSet$dataSet$paired <- mSet.old$dataSet$paired
                         mSet
                       }
        ) 
        
        if(mSet$dataSet$paired){
          mSet <- MetaboShiny::pair.mSet(mSet)
        }
        
        
        # ===== FILTERING ======
        if(mSet$metshiParams$filt_type != "none"){
          shiny::showNotification("Filtering dataset...")
          mSet$analSet <- mSet$storage$orig$analysis
          # TODO; add option to only keep columns that are also in QC ('qcfilter'?)
          mSet <- MetaboAnalystR::FilterVariable(mSet,
                                                 filter = mSet$metshiParams$filt_type,
                                                 qcFilter = "F",
                                                 rsd = 25)
          keep.mz <- colnames(mSet$dataSet$filt)
          mSet <- MetaboShiny::filt.mSet(mSet, keep.mz)
        }
        
        # =======================
        
        shiny::updateCheckboxInput(session, 
                                   "paired", 
                                   value = mSet$dataSet$paired) 
        
        if(grepl(mSet$dataSet$exp.type, pattern = "^1f")){
          if(mSet$dataSet$cls.num == 2){
            mSet$dataSet$exp.type <- "1fb"
          }else{
            mSet$dataSet$exp.type <- "1fm"
          }  
        }
        
        new.name = MetaboShiny::name.mSet(mSet)
        
        if(new.name %in% names(mSet$storage)){
          mSet <- MetaboShiny::load.mSet(mSet, new.name)
        }else{
          mSet$analSet <- list()
        }
        
        mSet$dataSet$cls.name <- new.name
        lcl$last_mset <- mSet$dataSet$cls.name
        shiny::showNotification("Changed mSet...")
        
        if(typeof(mSet) != "double"){
          success = T
        }
        
      })
      
      if(success){
        mSet <<- mSet
      }else{
        mSet.debug <<- mSet
        MetaboShiny::metshiAlert("Failed! Restoring old mSet...")
        mSet <<- mSet.old
      }
      mSetter$do <- NULL
    }
  }
})