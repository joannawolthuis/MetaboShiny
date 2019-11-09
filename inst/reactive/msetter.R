shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    isolate({
      
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
            mSet <- switch(mSetter$do,
                           change = {
                             mSet <- MetaboShiny::change.mSet(mSet, 
                                                              stats_var = input$stats_var, 
                                                              stats_type = input$stats_type, 
                                                              time_var = input$time_var)
                             mSet$dataSet$cls.name <- fut.name
                             
                             if(grepl("t", input$stats_type)){
                               mSet <- MetaboShiny::pair.mSet(mSet, fut.name) 
                               shiny::updateCheckboxInput(session, "paired", value = T)
                             }
                             mSet
                           },
                           subset = {
                             mSet <- MetaboShiny::subset.mSet(mSet, 
                                                              subset_var = input$subset_var, 
                                                              subset_group = input$subset_group)
                             mSet$dataSet$cls.name <- fut.name
                             
                             mSet
                           },
                           unsubset = {
                             mSet.old <- mSet
                             mSet <- MetaboShiny::load.mSet(mSet,"orig")
                             mSet <- MetaboShiny::change.mSet(mSet, 
                                                              stats_var = mSet.old$dataSet$exp.var, 
                                                              time_var =  mSet.old$dataSet$time.var,
                                                              stats_type = mSet.old$dataSet$exp.type)
                             mSet$dataSet$cls.name <- fut.name
                             
                             mSet
                           },
                           pair = {
                             # paired
                             mSet <- MetaboShiny::pair.mSet(mSet, fut.name)
                             mSet
                           }
            ) 
            
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
          
          success = T
        
          })

        if(success){
          mSet <<- mSet
        }else{
          MetaboShiny::metshiAlert("Failed! Restoring old mSet...")
          mSet <<- mSet.old
        }
        mSetter$do <- NULL
      }  
    })
    }
})