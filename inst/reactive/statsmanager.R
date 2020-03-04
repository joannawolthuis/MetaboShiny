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
        switch(statsmanager$calculate,
               vennrich = {
                 if("storage" %not in% names(mSet)){
                   mSet$storage <- list()
                 }
                 # TODO: use this in venn diagram creation
                 mSet <- MetaboShiny::store.mSet(mSet)
               },
               pattern = {
                 # pearson kendall spearman
                 NULL
               },
               pca = {
                 NULL
               },
               meba = {
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::performMB(mSet, 10) # perform MEBA analysis
                 })
               },
               asca = {
                 # perform asca analysis
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::Perform.ASCA(mSet, 1, 1, 2, 2)
                   mSet <- MetaboAnalystR::CalculateImpVarCutoff(mSet, 0.05, 0.9)
                 })
               },
               heatmap = {
                 # reset
                 mSet <- MetaboShiny::calcHeatMap(mSet, 
                                                  signif.only = input$heatsign,
                                                  source.anal = input$heattable,
                                                  top.hits = input$heatmap_topn,
                                                  cols = lcl$aes$mycols)
               },
               tt = {
                 withProgress({
                   mSet <- MetaboAnalystR::Ttests.Anal(mSet,
                                                       nonpar = input$tt_nonpar,
                                                       threshp = 0.1, # TODO: make the threshold user defined...
                                                       paired = mSet$dataSet$paired,
                                                       equal.var = input$tt_eqvar
                   )
                 })
               },
               fc = {
                 withProgress({
                   if(mSet$dataSet$paired){
                     mSet <- MetaboAnalystR::FC.Anal.paired(mSet,
                                                             1.5, # TODO: make this threshold user defined
                                                             1)  
                   }else{
                     mSet <- MetaboAnalystR::FC.Anal.unpaired(mSet,
                                                               1.5, # TODO: make this threshold user defined
                                                               1) 
                   }
                   if(!is.null(mSet$analSet$fc$sig.mat)){
                     rownames(mSet$analSet$fc$sig.mat) <- gsub(rownames(mSet$analSet$fc$sig.mat), 
                                                               pattern = "^X",
                                                               replacement = "")
                   }else{
                     mSet$analSet$fc$sig.mat <- data.frame()
                   }
                 })
               },
               aov = {
                 aovtype = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
                 redo = aovtype %not in% names(mSet$analSet)
                 if(redo){ # if done, don't redo
                   shiny::withProgress({
                     mSet <- switch(mSet$dataSet$exp.type,
                                    "1fm"=MetaboAnalystR::ANOVA.Anal(mSet, thresh=0.1,post.hoc = "fdr",nonpar = F),
                                    "2f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "", 1, 1),
                                    "t"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "time0", 1, 1),
                                    "t1f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.1, "fdr", "time", 1, 1))
                   })
                 }
               },
               volc = {
                 shiny::withProgress({
                   mSet <-  MetaboAnalystR::Volcano.Anal(mSet,
                                                         paired = mSet$dataSet$paired, 
                                                         1.5, 0,
                                                         0.75, F, 0.1,
                                                         TRUE, "raw") # TODO: make thresholds user-defined
                 })
               },
               power = {
                 NULL
               },
               wordcloud = {
                 if(input$wordcloud_manual){
                   shiny::withProgress({ # progress bar
                     abstracts = MetaboShiny::getAbstracts(searchTerms = input$wordcloud_searchTerm,
                                                           mindate = input$wordcloud_dateRange[1],
                                                           maxdate = input$wordcloud_dateRange[2],
                                                           retmax = input$wordcloud_absFreq)
                     lcl$tables$abstracts <- abstracts
                     shiny::setProgress(0.5)
                     lcl$tables$wordcloud_orig <- MetaboShiny::getWordFrequency(abstracts$abstract)
                     lcl$tables$wordcloud_filt <- lcl$tables$wordcloud_orig
                   }, message = "Searching...", max = 1)
                 }else{
                   if(nrow(shown_matches$forward_full) > 0){
                     lcl$tables$wordcloud_orig <- MetaboShiny::getWordFrequency(shown_matches$forward_full$description)
                     lcl$tables$wordcloud_filt <- lcl$tables$wordcloud_orig
                   }
                 }
               }) 
        if(typeof(mSet) != "double"){
          success <- T
        }
      })
      
      if(success){
        mSet <<- mSet
        lcl$hasChanged <<- TRUE
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_plots]"))
        shinyjs::show(selector = paste0("div.panel[value=collapse_", statsmanager$calculate, "_tables]"))
        shinyBS::updateCollapse(session, paste0("collapse_",input$statistics),open = paste0("collapse_", 
                                                                                            statsmanager$calculate, 
                                                                                            c("_tables","_plots")))
      }else{
        MetaboShiny::metshiAlert("Analysis failed!")
        mSet <<- mSet.old
      }
      lcl <<- lcl
    }
    # - - - -
    statsmanager$calculate <- NULL # set reloading to 'off'
  }
})
