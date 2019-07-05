# create listener
statsmanager <- reactiveValues()

observe({

  if(is.null(statsmanager$calculate)){

    NULL # if not reloading anything, nevermind

  }else{

    if(!is.null(mSet)){

      switch(statsmanager$calculate,
             venn = {
               # save previous mset
               mset_name = mSet$dataSet$cls.name
               # TODO: use this in venn diagram creation
               mSet$storage[[mset_name]] <<- list(analysis = mSet$analSet)
             },
             pca = {
               withProgress({
                 mSet <<- PCA.Anal(mSet) # perform PCA analysis
               })
             },
             meba = {
               withProgress({
                 mSet <<- performMB(mSet, 10) # perform MEBA analysis
               })

             },
             asca = {
               # perform asca analysis
               withProgress({
                 mSet <<- Perform.ASCA(mSet, 1, 1, 2, 2)
                 mSet <<- CalculateImpVarCutoff(mSet, 0.05, 0.9)
               })

             },
             heatmap = {
               # reset

               mSet$analSet$heatmap <<- NULL

               # change top hits used in heatmap depending on time series / bivariate / multivariate mode
               # reordering of hits according to most significant at the top
               if(interface$mode == "bivar"){

                 if(!is.null(input$heatmode)){

                   if(input$heatmode){

                     if("tt" %in% names(mSet$analSet)){
                       if(input$heatsign){
                         tbl <- req(as.data.frame(mSet$analSet$tt$sig.mat))
                       }else{
                         tbl <- data.frame('p.value' = mSet$analSet$tt$p.value)
                       }

                       used.values <- "p.value"

                       decreasing <- F

                     }else{

                       NULL

                     }


                   }else{

                     if("fc" %in% names(mSet$analSet)){
                       
                       if(input$heatsign){
                         tbl <- req(as.data.frame(mSet$analSet$fc$sig.mat))
                         tbl$abs_log2 <- abs(tbl$`log2(FC)`)
                       }else{
                         tbl <- data.frame("abs_log2" = abs(mSet$analSet$fc$fc.log))
                       }

                       used.values <- "abs_log2"
                       decreasing <- T
                     }else{
                       NULL
                     }

                   }
                 }
               }else if(interface$mode == "multivar"){

                 print("multivar mode")

                 tbl <- as.data.frame(mSet$analSet$aov$sig.mat)
                 used.values <- "p.value"
                 decreasing <- F

               }else{
                 if(!is.null(input$heatmode)){
                   
                 if(input$heatmode){
                   tbl <- as.data.frame(mSet$analSet$asca$sig.list$Model.ab)
                   used.values <- "Leverage"
                 }else{
                   tbl <- as.data.frame(mSet$analSet$MB$stats)
                   used.values <- "Hotelling-T2"
                 }
                 decreasing = T

               if(!exists("tbl")) return(NULL)
               if(is.null(tbl)) return(NULL)
               if(nrow(tbl) == 0 ) return(NULL)

               # check top x used (slider bar in UI), if more than total matches use total matches
               topn = if(length(tbl[[used.values]]) < input$heatmap_topn) length(tbl[[used.values]]) else input$heatmap_topn
               mzorder <- order(tbl[[used.values]], decreasing = decreasing)
               mzsel <- rownames(tbl)[mzorder]#[1:topn]

               # reorder matrix used
               x <- mSet$dataSet$norm[,mzsel]
               final_matrix <- t(x) # transpose so samples are in columns

               # check if the sample order is correct - mSet$..$ norm needs to match the matrix
               sample_order <- match(colnames(final_matrix), rownames(mSet$dataSet$norm))

               if(timebutton$status == "on"){ # check if time series
                 if(input$timecourse_trigger){

                   # create convenient table with the ncessary info
                   translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$exp.fac[sample_order], Time=mSet$dataSet$time.fac[sample_order])
                   hmap.lvls <- c(levels(mSet$dataSet$exp.fac), levels(mSet$dataSet$time.fac))

                   # reorder first by time, then by sample
                   split.translator <- split(translator, by = c("Time"))
                   split.translator.ordered <- lapply(split.translator, function(tbl) tbl[order(tbl$Group)])
                   translator <- rbindlist(split.translator.ordered)

                   # ensure correct sample order
                   final_matrix <- final_matrix[,match(translator$Sample, colnames(final_matrix))]

                   # disable automatic ordering of samples through clustering
                   my_order=F

                 }else{
                   # no complicated reordering necessary
                   translator <- data.table(Sample=as.character(rownames(mSet$dataSet$norm))[sample_order],Group=mSet$dataSet$cls[sample_order])
                   hmap.lvls <- levels(mSet$dataSet$cls)
                   my_order = T # enable sorting through dendrogram
                 }
               }else{
                 # no complicated reordering necessary
                 translator <- data.table(Sample=rownames(mSet$dataSet$norm)[sample_order],Group=mSet$dataSet$cls[sample_order])
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
               }
               }
             },
             tt = {
               withProgress({
                 mSet <<- Ttests.Anal(mSet,
                                      nonpar = input$tt_nonpar,
                                      threshp = 0.05, # TODO: make the threshold user defined...
                                      paired = FALSE,
                                      equal.var = input$tt_eqvar
                 )
               })
             },
             fc = {
               withProgress({
                 mSet <<- FC.Anal.unpaired(mSet,
                                           2.0, # TODO: make this threshold user defined
                                           1)
               })
             },
             aov = {
               
               if(!is.null(input$timecourse_trigger)){
                 redo = switch(input$timecourse_trigger,
                               {"aov2" %not in% names(mSet$analSet)},
                               {"aov" %not in% names(mSet$analSet)})
               }else{
                 redo = "aov" %not in% names(mSet$analSet)
               }

               print(redo)
               
               if(redo){ # if done, don't redo
                 withProgress({
                   if(!is.null(input$timecourse_trigger)){
                     mSet <<- if(input$timecourse_trigger){
                       ANOVA2.Anal(mSet, 0.05, "fdr", "time", 3, 1)
                     }else{
                       ANOVA.Anal(mSet, thresh=0.05,nonpar = F)
                     }
                   }else{
                     mSet <<- ANOVA.Anal(mSet, thresh=0.05,nonpar = F)
                   }
                 })
               }
             },
             volc = {
               withProgress({
                 mSet <<- Volcano.Anal(mSet,
                                       FALSE, 2.0, 0,
                                       0.75,F, 0.1,
                                       TRUE, "raw") # TODO: make thresholds user-defined
               })
             },
             match_wordcloud = {
               if(nrow(lcl$tables$last_matches) > 0){
                 try({
                   
                   # remove unwanted words (defined in global) from description
                   filtered_descriptions <- sapply(1:length(lcl$tables$last_matches$description),
                                                   function(i){
                                                     # get description
                                                     desc <- lcl$tables$last_matches$description[[i]]
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
             },
             match_pie = {
               
               if(nrow(lcl$tables$last_matches) > 0){
               adduct_dist <- melt(table(lcl$tables$last_matches$adduct))
               db_dist <- melt(table(lcl$tables$last_matches$source))
               
               lcl$vectors$pie_add <<- adduct_dist
               lcl$vectors$pie_db <<- db_dist
               
               }
             })
    }
    # - - - -
    statsmanager$calculate <- NULL # set reloading to 'off'
  }
})
