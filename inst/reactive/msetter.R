shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    isolate({
      
      if(!is.null(mSet)){
        
        print(mSetter$do)
        
        mSet <- MetaboShiny::store.mSet(mSet)
        
        mSet <- prev.mSet(mSet, 
                          what = mSetter$do,
                          input = list(stats_var = input$stats_var, 
                                       stats_type = input$stats_type, 
                                       time_var = input$time_var,
                                       subset_var = input$subset_var, 
                                       subset_group = input$subset_group,
                                       paired = input$paired))
        
        fut.name = name.mSet(mSet)
        
        print(fut.name)
        
        if(fut.name %in% names(mSet$storage)){
          mSet <- load.mSet(mSet, fut.name)
        }else{
          mSet <- switch(mSetter$do,
                         change = {
                           mSet <- MetaboShiny::change.mSet(mSet, 
                                                            stats_var = input$stats_var, 
                                                            stats_type = input$stats_type, 
                                                            time_var = input$time_var)
                           mSet$dataSet$cls.name <- fut.name
                           
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
                           overview.tbl <- data.table::as.data.table(cbind(sample = mSet$dataSet$covars$sample,
                                                                           individual = mSet$dataSet$covars$individual, 
                                                                           variable = as.character(mSet$dataSet$cls)))
                           uniq.overv <- data.table::as.data.table(unique(overview.tbl[,-1]))
                           keep.indiv <- uniq.overv$individual[c(which(duplicated(uniq.overv$individual, fromLast=T)), 
                                                                 which(duplicated(uniq.overv$individual)))]
                           samp.overv <- mSet$dataSet$covars[individual %in% keep.indiv, c("individual","sample")]
                           keep.samples = caret::downSample(samp.overv$sample,as.factor(samp.overv$individual))$x
                           fin.overv <- overview.tbl[sample %in% keep.samples,]
                           keep.samples = caret::downSample(fin.overv$sample,as.factor(fin.overv$variable))$x
                           
                           if(length(keep.samples)> 2){
                             mSet <- MetaboShiny::subset.mSet(mSet, 
                                                              subset_var = "sample", 
                                                              subset_group = keep.samples)
                             mSet$dataSet$cls.name <- fut.name
                           }else{
                             print("Not enough samples for paired analysis!! Doing regular switching...")
                             mSet$dataSet$paired <- FALSE  
                           }  
                           
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
          mSet$analSet <- NULL
        }
        
        lcl$last_mset <- mSet$dataSet$cls.name
        print("changed mSet...")
        mSet <<- mSet
        mSetter$do <- NULL
      }  
    })
    }
})