# create listener
statsmanager <- shiny::reactiveValues()

shiny::observe({
  
  if(is.null(statsmanager$calculate)){
    
    NULL # if not reloading anything, nevermind
    
  }else{
    
    if(!is.null(mSet)){
      
      mSet.old <- mSet
      success = F
      try({
        results <- runStats(mSet,
                            input = input,
                            lcl = lcl,
                            analysis = statsmanager$calculate,
                            ml_queue = ml_queue)
        
        mSet <- results$mSet
        lcl <- results$lcl
        
        if(typeof(mSet) != "double"){
          success <- T
        }
      })
      
      if(success){
        mSet <<- mSet
        save_info$has_changed <- TRUE
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_plots]"))
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_tables]"))
        shinyBS::updateCollapse(session, paste0("collapse_",input$statistics),open = paste0("collapse_", 
                                                                                            statsmanager$calculate, 
                                                                                            c("_tables","_plots")))
        if(lcl$beep){
          beepr::beep(sound = lcl$aes$which_beep)
          Sys.sleep(0.6)
          beepr::beep(sound = lcl$aes$which_beep)
          Sys.sleep(0.6)
          beepr::beep(sound = lcl$aes$which_beep)
        }
        
      }else{
        MetaboShiny::metshiAlert("Analysis failed!")
        #shiny::showNotification(msg.vec)
        mSet <<- mSet.old
      }
      lcl <<- lcl
    }
    # - - - -
    statsmanager$calculate <- NULL # set reloading to 'off'
  }
})