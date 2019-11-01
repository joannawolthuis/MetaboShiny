# create listener
statsmanager <- shiny::reactiveValues()

shiny::observe({

  if(is.null(statsmanager$calculate)){

    NULL # if not reloading anything, nevermind

  }else{
    
    print(statsmanager$calculate)

    if(!is.null(mSet)){

      switch(statsmanager$calculate,
             venn = {
               # save previous mset
               mset_name = mSet$dataSet$cls.name
               # TODO: use this in venn diagram creation
               mSet$storage[[mset_name]] <<- list(analysis = mSet$analSet)
             },
             pca = {
               shiny::withProgress({
                 mSet <<- MetaboAnalystR::PCA.Anal(mSet) # perform PCA analysis
               })
             },
             meba = {
               shiny::withProgress({
                 mSet <<- MetaboAnalystR::performMB(mSet, 10) # perform MEBA analysis
               })

             },
             asca = {
               try({
                 # perform asca analysis
                 shiny::withProgress({
                   mSet <<- MetaboAnalystR::Perform.ASCA(mSet, 1, 1, 2, 2)
                   mSet <<- MetaboAnalystR::CalculateImpVarCutoff(mSet, 0.05, 0.9)
                 })
               })
             },
             heatmap = {
               # reset

               mSet$analSet$heatmap <<- NULL

               if(!(input$heattable %in% names(mSet$analSet))) return(NULL)
               switch(input$heattable,
                             tt={
                               decreasing <- F
                               if(input$heatsign){
                                 tbl=as.data.frame(mSet$analSet$tt$sig.mat)
                                 sigvals = (tbl$p.value)
                               }else{
                                 tbl = as.data.frame(mSet$analSet$tt$p.value)
                                 sigvals = tbl[,1]
                               }
                             },
                             fc={
                               decreasing <- T
                               
                               if(input$heatsign){
                                 tbl <- as.data.frame(mSet$analSet$fc$sig.mat)
                                 sigvals=abs(tbl$`log2(FC)`)
                               }else{
                                 tbl = as.data.frame(abs(mSet$analSet$fc$fc.log))
                                 sigvals = tbl[,1]
                               }
                             },
                             aov= {
                               decreasing=F
                               tbl=as.data.frame(mSet$analSet$aov$sig.mat)
                               sigvals = tbl$p.value
                             },
                             aov2={
                               decreasing=F
                               tbl = as.data.frame(mSet$analSet$aov2$sig.mat)
                               sigvals = tbl$`Interaction(adj.p)`
                             },
                             asca={
                               decreasing=T
                               tbl = as.data.frame(mSet$analSet$asca$sig.list$Model.ab)
                               sigvals = tbl$Leverage
                             },
                             meba={
                               decreasing=T
                               tbl = as.data.frame(mSet$analSet$MB$stats)
                               sigvals = tbl$`Hotelling-T2`
                             })
               # change top hits used in heatmap depending on time series / bivariate / multivariate mode
               # reordering of hits according to most significant at the top
            
               if(!exists("sigvals")) return(NULL)
               if(is.null(sigvals)) return(NULL)
               if(length(sigvals) == 0 ) return(NULL)

               #check top x used (slider bar in UI), if more than total matches use total matches
               topn = if(length(sigvals) < input$heatmap_topn) length(sigvals) else input$heatmap_topn
               mzorder <- order(sigvals, decreasing = decreasing)
               mzsel <- rownames(tbl)[mzorder]#[1:topn]

               # reorder matrix used
               x <- mSet$dataSet$norm[,mzsel]
               final_matrix <- t(x) # transpose so samples are in columns

               # check if the sample order is correct - mSet$..$ norm needs to match the matrix
               sample_order <- match(colnames(final_matrix), rownames(mSet$dataSet$norm))

               if(mSet$dataSet$exp.type %in% c("2f", "t1f", "t")){

                   # create convenient table with the ncessary info
                   translator <- data.table::data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],GroupA=mSet$dataSet$facA[sample_order], GroupB=mSet$dataSet$facB[sample_order])
                   hmap.lvls <- c(levels(mSet$dataSet$facA), levels(mSet$dataSet$facB))

                   # reorder first by time, then by sample
                   split.translator <- split(translator, by = c("GroupB"))
                   split.translator.ordered <- lapply(split.translator, function(tbl) tbl[order(tbl$GroupA)])
                   translator <- data.table::rbindlist(split.translator.ordered)

                   # ensure correct sample order
                   final_matrix <- final_matrix[,match(translator$Sample, colnames(final_matrix))]

                   # disable automatic ordering of samples through clustering
                   my_order=F

               }else{
                 # no complicated reordering necessary
                 translator <- data.table::data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],
                                                      Group=mSet$dataSet$cls[sample_order])
                 hmap.lvls <- levels(mSet$dataSet$cls)
                 my_order = T # enable sorting through dendrogram
               }

               # create name - to - color mapping vector for the plotting functions
               color.mapper <- {
                 classes <- hmap.lvls
                 cols <- sapply(1:length(classes), function(i) lcl$aes$mycols[i]) # use user-defined colours
                 names(cols) <- classes
                 # - - -
                 cols
               }

               # write to mSet
               mSet$analSet$heatmap <<- list(
                 matrix = final_matrix,
                 translator = translator,
                 colors = color.mapper,
                 my_order = my_order)
             },
             tt = {
               withProgress({
                 mSet <<- MetaboAnalystR::Ttests.Anal(mSet,
                                      nonpar = input$tt_nonpar,
                                      threshp = 0.05, # TODO: make the threshold user defined...
                                      paired = FALSE,
                                      equal.var = input$tt_eqvar
                 )
               })
             },
             fc = {
               withProgress({
                 mSet <<- MetaboAnalystR::FC.Anal.unpaired(mSet,
                                           2.0, # TODO: make this threshold user defined
                                           1)
               })
             },
             aov = {
               
              aovtype = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
              redo = aovtype %not in% names(mSet$analSet)
               
               if(redo){ # if done, don't redo
                 shiny::withProgress({
                   mSet <<- switch(mSet$dataSet$exp.type,
                                   "1f"=MetaboAnalystR::ANOVA.Anal(mSet, thresh=0.05,nonpar = F),
                                   "2f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.05, "fdr", "", 1, 1),
                                   "t"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.05, "fdr", "time0", 1, 1),
                                   "t1f"=MetaboAnalystR::ANOVA2.Anal(mSet, 0.05, "fdr", "time", 3, 1))
                 })
               }
             },
             volc = {
               shiny::withProgress({
                 mSet <<- MetaboAnalystR::Volcano.Anal(mSet,
                                       FALSE, 2.0, 0,
                                       0.75,F, 0.1,
                                       TRUE, "raw") # TODO: make thresholds user-defined
               })
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
                   
                   lcl$tables$pm_absdata <<- tbl
                   
                   shiny::setProgress(0.6)
                   
                   wcdata <- res
                   colnames(wcdata) <- c("name", "value")
                 }
                 
                 lcl$tables$word_freq_pm <<- wcdata
                 
               })
             },
             match_wordcloud = {
               
               if(nrow(shown_matches$forward_full) > 0){
                
                  try({
                   
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
                   
                   lcl$tables$word_freq <<- data.frame(name = names(v), value = v)
                   
                 })
               }
             })
    }
    # - - - -
    statsmanager$calculate <- NULL # set reloading to 'off'
  }
})
