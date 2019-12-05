shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    
    if(!is.null(mSet)){
      
      print(mSetter$do)
      
      mSet.old <- mSet
      
      mSet <- MetaboShiny::store.mSet(mSet)
      
      mSet <- MetaboShiny::prev.mSet(mSet, 
                                     what = mSetter$do,
                                     input = list(stats_var = input$stats_var, 
                                                  stats_type = input$stats_type, 
                                                  time_var = input$time_var,
                                                  subset_var = input$subset_var, 
                                                  subset_group = input$subset_group,
                                                  paired = input$paired))
      
      fut.name = MetaboShiny::name.mSet(mSet)
      
      success = F
      
      try({
        if(fut.name %in% names(mSet$storage)){
          mSet <- MetaboShiny::load.mSet(mSet, fut.name)
        }else{
          # - - - - - - - - - - - - -
          mSet <- switch(mSetter$do,
                         change = {
                           was.paired <- mSet.old$dataSet$paired
                           was.timeseries <- grepl(mSet.old$dataSet$exp.type, "t") # also was paired
                           if(was.paired){
                             if(was.timeseries){
                               mSet <- MetaboShiny::load.mSet(mSet, "orig") 
                             }else{
                               nonpaired.name <- gsub(mSet.old$dataSet$cls.name, pattern = if(was.timeseries) "\\(timeseries\\)" else  "\\(paired\\)", replacement = "")
                               mSet <- MetaboShiny::load.mSet(mSet, nonpaired.name) 
                             }
                           }
                           
                           mSet <- MetaboShiny::change.mSet(mSet, 
                                                            stats_var = input$stats_var, 
                                                            stats_type = input$stats_type, 
                                                            time_var = input$time_var)
                           if(mSet$dataSet$paired){
                             mSet <- MetaboShiny::pair.mSet(mSet, fut.name)
                           }
                           mSet$dataSet$cls.name <- fut.name
                           mSet
                         },
                         subset = {
                           mSet <- MetaboShiny::subset.mSet(mSet, 
                                                            subset_var = input$subset_var, 
                                                            subset_group = input$subset_group)
                           mSet$dataSet$paired <- FALSE
                           mSet$dataSet$cls.name <- fut.name
                           mSet
                         },
                         unsubset = {
                           mSet <- MetaboShiny::load.mSet(mSet,"orig")
                           mSet <- MetaboShiny::change.mSet(mSet, 
                                                            stats_var = mSet.old$dataSet$exp.var, 
                                                            time_var =  mSet.old$dataSet$time.var,
                                                            stats_type = mSet.old$dataSet$exp.type)
                           mSet$dataSet$paired <- FALSE
                           mSet$dataSet$cls.name <- fut.name
                           mSet
                         }
          ) 
          
          # if(mSetter$do %in% c("subset", "unsubset")){
          #   mSet <- MetaboAnalystR::Normalization(mSet,
          #                                         rowNorm = mSet$metshiParams$norm_type,
          #                                         transNorm = mSet$metshiParams$trans_type,
          #                                         scaleNorm = mSet$metshiParams$scale_type,
          #                                         ref = mSet$metshiParams$ref_var)
          #   
          # }
          
          shiny::updateCheckboxInput(session, "paired", value = mSet$dataSet$paired) 
          
          mSet$dataSet$cls.name <- fut.name
          if(grepl(mSet$dataSet$exp.type, pattern = "^1f")){
            if(mSet$dataSet$cls.num == 2){
              mSet$dataSet$exp.type <- "1fb"
            }else{
              mSet$dataSet$exp.type <- "1fm"
            }  
          }
          mSet$analSet <- list()
        }
        
        lcl$last_mset <- mSet$dataSet$cls.name
        shiny::showNotification("Changed mSet...")
        
        if(typeof(mSet) != "double"){
          success = T
        }
        
      })
      
      if(success){
        mSet <<- mSet
      }else{
        MetaboShiny::metshiAlert("Failed! Restoring old mSet...")
        mSet <<- mSet.old
      }
      mSetter$do <- NULL
    }
  }
})