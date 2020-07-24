shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    
    if(!is.null(mSet)){
      
      mSet <- MetaboShiny::store.mSet(mSet) # save analyses
      
      oldSettings <- mSet$settings
      
      # reset dataset
      mSet <- MetaboShiny::reset.mSet(mSet,
                                      fn = file.path(lcl$paths$proj_dir, 
                                                     paste0(lcl$proj_name,
                                                            "_ORIG.metshi")))
      
      orig.count <- mSet$settings$orig.count
      
      success = F
      
      try({
        if(!(mSetter$do %in% c("unsubset"))){
          mSet.settings <- if(mSetter$do == "load") mSet$storage[[input$storage_choice]]$settings else oldSettings
          if(length(mSet.settings$subset) > 0){
            subs = mSet.settings$subset
            subs = subs[!(names(subs) %in% c("sample", "mz"))]
            if(length(subs) > 0){
              for(i in 1:length(subs)){
                mSet <- MetaboShiny::subset_mSet(mSet, 
                                                 subset_var = names(subs)[i], 
                                                 subset_group = subs[[i]])  
              }  
            }
          }
        }else{
          mSet.settings <- oldSettings
        }
        
        mSet <- switch(mSetter$do,
                       refresh = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.settings$exp.var, 
                                                          time_var =  mSet.settings$time.var,
                                                          stats_type = mSet.settings$exp.type)
                         mSet$dataSet$paired <- mSet.settings$paired
                         mSet
                       },
                       load = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.settings$exp.var, 
                                                          time_var =  mSet.settings$time.var,
                                                          stats_type = mSet.settings$exp.type)
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
                             mSet <- MetaboShiny::subset_mSet(mSet,
                                                              subset_var = "sample", 
                                                              subset_group = knowns) 
                           }
                         }
                         mSet
                       },
                       subset = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.settings$exp.var, 
                                                          time_var =  mSet.settings$time.var,
                                                          stats_type = mSet.settings$exp.type)
                         mSet <- MetaboShiny::subset_mSet(mSet,
                                                          subset_var = input$subset_var, 
                                                          subset_group = input$subset_group)
                         mSet$dataSet$paired <- mSet.settings$paired
                         mSet
                       },
                       unsubset = {
                         mSet <- MetaboShiny::change.mSet(mSet, 
                                                          stats_var = mSet.settings$exp.var, 
                                                          time_var =  mSet.settings$time.var,
                                                          stats_type = mSet.settings$exp.type)
                         mSet$dataSet$paired <- mSet.settings$paired
                         mSet
                       }) 
        
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
        new.name = if(mSetter$do == "load") input$storage_choice else name.mSet(mSet)
        
        if(new.name %in% names(mSet$storage)){
          mSet <- MetaboShiny::load.mSet(mSet, new.name)
        }else{
          mSet$analSet <- list(type = "stat")
        }
        
        mSet$dataSet$cls.name <- new.name
        mSet$settings$orig.count <- orig.count
        shiny::updateCheckboxInput(session, 
                                   "paired", 
                                   value = mSet$dataSet$paired) 
        lcl$last_mset <- mSet$dataSet$cls.name
        success = T
        shiny::showNotification("Changed mSet...")
      })
      
      if(success){
        # filtering?
        if(mSet$settings$filt.type != "none"){
          shiny::showNotification("Filtering dataset...")
          # TODO; add option to only keep columns that are also in QC ('qcfilter'?)
          keep.mz <- colnames(MetaboAnalystR::FilterVariable(mSet,
                                                             filter = mSet$settings$filt.type,
                                                             qcFilter = "F",
                                                             rsd = 25)$dataSet$filt)
          mSet$dataSet$norm <- mSet$dataSet$norm[,keep.mz]
          mSet$dataSet$proc <- mSet$dataSet$proc[,keep.mz]
        }
        
        # if(mSet$metshiParams$prematched & input$show_prematched_mz_only){
        #   orig.mz <- ncol(mSet$dataSet$norm)
        #   conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb)
        #   
        #   showadd = input$prematch_show_add
        #   showiso = input$prematch_show_iso
        #   
        #   showadd <- if(length(showadd) > 0) paste0(showadd, collapse=" OR adduct = '", "'") 
        #   showiso <- if(length(showiso) > 0) if(length(showiso) == 2) c() else showiso
        #   
        #   firstpart = "SELECT DISTINCT query_mz FROM match_mapper"
        #   addfrag = if(length(showadd) > 0) gsubfn::fn$paste("AND (adduct = '$showadd)") else ""
        #   isofrag = if(length(showiso) > 0) switch(showiso, 
        #                                          main = "AND `%iso` > 99.9999", 
        #                                          minor = "AND `%iso` < 99.9999") else ""
        #   
        #   query = gsubfn::fn$paste("$firstpart WHERE query_mz IS NOT NULL $addfrag $isofrag")
        # 
        #   matched_mz = DBI::dbGetQuery(conn, query)[,1]
        #   DBI::dbDisconnect(conn)
        #   mSet <- subset_mSet_mz(mSet, matched_mz)
        #   
        #   new.mz <- ncol(mSet$dataSet$norm)
        #   # update on screen
        #   output$mz_count <- shiny::renderText({
        #     paste0("m/z: ", 
        #            as.character(new.mz), if(orig.mz == new.mz) "" else paste0("/",as.character(orig.mz)))
        #   })
        # }else{
        #   output$mz_count <- shiny::renderText("")
        # }
        
        # =========
        
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
        lcl$hasChanged <<- FALSE
      }else{
        MetaboShiny::metshiAlert("Failed! Restoring old mSet...")
      }
      mSetter$do <- NULL
      }
  }
})