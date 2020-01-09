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
               venn = {
                 # save previous mset
                 mset_name = mSet$dataSet$cls.name
                 # TODO: use this in venn diagram creation
                 mSet$storage[[mset_name]] <-  list(analysis = mSet$analSet)
               },
               pattern = {
                 # pearson kendall spearman
                 NULL
               },
               pca = {
                 shiny::withProgress({
                   mSet <- MetaboAnalystR::PCA.Anal(mSet) # perform PCA analysis
                 })
               },
               meba = {
                 shiny::withProgress({
                   mSet <-  MetaboAnalystR::performMB(mSet, 10) # perform MEBA analysis
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
                 mSet$analSet$heatmap <-  NULL
                
                 mSet <- MetaboShiny::calcHeatMap(mSet, 
                                                  signif.only = input$heatsign,
                                                  source.tbl = input$heattable,
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
                     mSet <-  MetaboAnalystR::FC.Anal.paired(mSet,
                                                             1.5, # TODO: make this threshold user defined
                                                             1)  
                   }else{
                     mSet <-  MetaboAnalystR::FC.Anal.unpaired(mSet,
                                                               1.5, # TODO: make this threshold user defined
                                                               1) 
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
               match_wordcloud_pm = {
                 shiny::withProgress({
                   abstr <- RISmed::EUtilsSummary(
                     as.character(input$pm_query),
                     type="esearch",
                     db="pubmed",
                     datetype='pdat',
                     mindate=input$pm_year[1],
                     maxdate=input$pm_year[2],
                     retmax=as.numeric(input$pm_max))
                   
                   shiny::setProgress(0.2)
                   
                   fetch <- RISmed::EUtilsGet(abstr)
                   wcdata = data.frame()
                   if(length(fetch@PMID) > 0){
                     res <- MetaboShiny::abstracts2wordcloud(abstracts = fetch,
                                                             top = input$wc_topn_pm, 
                                                             queryword = input$pm_query)
                     
                     shiny::setProgress(0.4)
                     
                     tbl <- data.frame(
                       pmids = RISmed::PMID(fetch),
                       titles = RISmed::ArticleTitle(fetch)
                       #abstracts = RISmed::AbstractText(fetch)
                     )
                     
                     lcl$tables$pm_absdata <-  tbl
                     
                     shiny::setProgress(0.6)
                     
                     wcdata <- res
                     colnames(wcdata) <- c("name", "value")
                   }
                   
                   lcl$tables$word_freq_pm <-  wcdata
                   
                 })
               },
               match_wordcloud = {
                 if(nrow(shown_matches$forward_full) > 0){
                   
                   # remove unwanted words (defined in global) from description
                   filtered_descriptions <- sapply(1:length(shown_matches$forward_full$description),
                                                   function(i){
                                                     # get description
                                                     desc <- shown_matches$forward_full$description[[i]]
                                                     # return
                                                     desc
                                                   })
                   
                   require(tm)
                   
                   docs <- tm::VCorpus(tm::VectorSource(filtered_descriptions))
                   
                   # Convert all text to lower case
                   docs <- tm::tm_map(docs, tm::content_transformer(tolower))
                   
                   # Remove punctuations
                   docs <- tm::tm_map(docs, tm::removePunctuation)
                   
                   # Remove numbers
                   docs <- tm::tm_map(docs, tm::removeNumbers)
                   
                   # # stem words
                   # docs <- tm::tm_map(docs, tm::stemDocument)
                   
                   # Remove english common stopwords
                   docs <- tm::tm_map(docs, 
                                      tm::removeWords, 
                                      gbl$vectors$wordcloud$skip)
                   
                   # Remove whitespace
                   docs <- tm::tm_map(docs, tm::stripWhitespace)
                   
                   # # Text stemming
                   
                   doc_mat <- tm::TermDocumentMatrix(docs)
                   
                   m <- as.matrix(doc_mat)
                   
                   v <- sort(rowSums(m), decreasing = TRUE)
                   
                   lcl$tables$word_freq <-  data.frame(name = names(v), value = v)
                   
                 }
               }) 
        success <- T
      })
      if(success){
        mSet <<- mSet
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
