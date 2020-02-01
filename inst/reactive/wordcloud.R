# normal wordcloud

observeEvent(input$wordcloud_go, {
  statsmanager$calculate <- "wordcloud"
  #datamanager$reload <- "wordcloud"
})

observeEvent(input$wordcloud_topWords, {
  #statsmanager$calculate <- "wordcloud"
  datamanager$reload <- "wordcloud"
})

observeEvent(input$wordcloud_filter, {
  if("wordcloud_orig" %in% names(lcl$tables)){
    # get all the lists of filter words that the user wants, and join them together into a big list
    filterList = c()
    try({
      filterList <- unlist(gbl$vectors$wordcloud$filters[input$wordcloud_filter])
    }, silent = T)
    # remove single character words
    # start= ^ , any character = . , end = $"
    singleChar <- grep(lcl$tables$wordcloud_orig$word, pattern = "^.{1,3}$", value = T) 
    # remove verbs ending on -es and -ed (differentiated, etc.)
    #verbs <- grep(lcl$tables$wordcloud_orig$word, pattern = ".*[ed|es]$", value = T) 
    # remove numbers (p-values and the like)
    suppressWarnings({
      numericals = lcl$tables$wordcloud_orig$word[which(!is.na(as.numeric(
        gsub(x = lcl$tables$wordcloud_orig$word, 
             pattern = ",", 
             replacement = ".")
      )))]  
    })
    # make an extra filter list for the stuff that comes specifically from this search term
    additionalFilters <- c(strsplit(input$wordcloud_searchTerm, 
                                    # dont include the words themselves
                                    split = " ")[[1]],
                           singleChar,
                           numericals
                           )
    # merge into final filter list
    filterList <- c(filterList, additionalFilters, fill = T)
    without_stopwords <- data.table::as.data.table(MetaboShiny::getFilteredWordFreqency(lcl$tables$wordcloud_orig, filterList))
    without_stopwords <- without_stopwords[without_stopwords$word != ""]
    lcl$tables$wordcloud_filt <<- without_stopwords
    datamanager$reload <- "wordcloud"    
  }
})

observeEvent(input$wordcloud_make_filter, {
  withProgress({
    filterAbstracts = MetaboShiny::getAbstracts(searchTerms = input$wordcloud_filterTerm,
                                   mindate = input$wordcloud_filterDateRange[1],
                                   maxdate = input$wordcloud_filterDateRange[2],
                                   retmax = input$wordcloud_filterAbsFreq)
    shiny::setProgress(0.5)
    topWords = MetaboShiny::getWordFrequency(filterAbstracts$abstract)
    filterFolder = file.path(lcl$paths$work_dir, "wordcloud")
    if(!dir.exists(filterFolder)){
      dir.create(filterFolder)
    }
    shiny::setProgress(0.7)
    filterList <- topWords[order(topWords$n, decreasing = TRUE)[1:input$wordcloud_filterTopN],]
    data.table::fwrite(x = filterList,
                       file = file.path(filterFolder, paste0(input$wordcloud_filterTerm, ".csv")))    
    datamanager$reload <- "general"
  })
})

observe({
  if(!is.null(input$wordcloud_selected_word)){
    selected_word_clean = gsub(pattern= ":\\d+",replacement = "",isolate(input$wordcloud_selected_word))
    abstracts <- lcl$tables$abstracts
    JustASelectedbstracts <- abstracts$abstract
    MachedAbstracts = sapply(JustASelectedbstracts, function(x) grepl(pattern = selected_word_clean, x = x))
    isWhere=which(MachedAbstracts)
    # ==========================================
    ShowAbstracts<-JustASelectedbstracts[isWhere]
    absTable = lcl$tables$abstracts[isWhere,]
    # ShowAbstractsandSummaries<-paste(ShowAbstractsandSummaries,
    #                                  collapse = "\n\n")
    # ShowAbstractsandSummaries<-sapply(ShowAbstractsandSummaries,FUN = as.character,USE.NAMES = T)
    # ==========================================
    output$wordcloud_abstracts <- renderUI({ 
      lapply(1:nrow(absTable), function(i){
        row = absTable[i, ]
        title = row$title#row[1]
        pmid = row$DOI#row[4]
        abstract = row[2]
        shiny::fluidRow(align="center", 
                 shiny::tags$b(title), 
                 shiny::a(gsubfn::fn$paste("($pmid)"), href=gsubfn::fn$paste("https://pubmed.ncbi.nlm.nih.gov/$pmid")),
                 br(),
                 shiny::tags$i(abstract),
                 hr())
        })
      }) 
  }
})

