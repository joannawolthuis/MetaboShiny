tabPanel(title="settings",
         fluidRow(align = "center",
                  shinyWidgets::switchInput(
                    inputId = "wordcloud_manual",
                    label = "Manual search", 
                    labelWidth = "80px"
                  ),
                  shiny::conditionalPanel("input.wordcloud_manual == true",
                                          helpText("Type a metabolite below and search PubMed to find documents that contain that word in the text."),
                                          textInput("wordcloud_searchTerm", label = h3("Enter your search terms"), placeholder = "enter your search terms"),
                                          helpText("You can specify the start and end years of your search, use the format YYYY"),
                                          sliderInput(inputId = "wordcloud_dateRange", label = "Date range input:",min = 2000, max = 2019,value = c(2010, 2020), sep = ""),
                                          sliderInput("wordcloud_absFreq", "Amount of abstracts to use:", min = 1,  max = 4000, value = 500)
                  ),
                  actionButton("wordcloud_go", "Go! :)")
         ))
tabPanel(title = "filters",
         helpText("Type a search term: a word filter will be created with the resulting top words."),
         textInput("wordcloud_filterTerm", label = h3("Enter your search terms"), placeholder = "enter your search terms"),
         helpText("You can specify the start and end years of your search, use the format YYYY"),
         sliderInput(inputId = "wordcloud_filterDateRange", label = "Date range input:",min = 2000, max = 2019,value = c(2010, 2020), sep = ""),
         sliderInput("wordcloud_filterAbsFreq", "Amount of abstracts to use:", min = 10,  max = 500000, value = 5000),
         actionButton("wordcloud_make_filter", "Go! :)")
)
tabPanel(title = "plot",
         wordcloud2::wordcloud2Output(outputId = "cloud"),
         sliderInput("topWords", "Top words used for plot:", min = 1,  max = 10000, value = 100),
         selectizeInput("wordcloud_filter",
                        "Filter",
                        choices = c("medical","stopwords","metabolomics"), 
                        multiple = T)
)

# =======================

# normal wordcloud

observeEvent(input$wordcloud_go, {
  statsmanager$calculate <- "wordcloud"
  #datamanager$reload <- "wordcloud"
})

observeEvent(input$wordcloud_filter, {
  if(nrow(lcl$tables$wordcloud_orig) > 0){
    shiny::showNotification("this is working!")
    # get all the lists of filter words that the user wants, and join them together into a big list
    filterList <- unlist(gbl$vectors$wordcloud$filters[input$wordcloud_filter])
    # remove single character words
    # start= ^ , any character = . , end = $"
    singleChar <- grep(lcl$tables$wordcloud_orig$word, pattern = "^.{1,3}$", value = T) 
    # remove verbs ending on -es and -ed (differentiated, etc.)
    verbs <- grep(lcl$tables$wordcloud_orig$word, pattern = ".*[ed|es]$", value = T) 
    # remove numbers (p-values and the like)
    suppressWarnings({
      numericals = lcl$tables$wordcloud_orig$word[which(!is.na(as.numeric(
        gsub(x = lcl$tables$wordcloud_orig$word, pattern = ",", replacement = ".")
      )))]  
    })
    # make an extra filter list for the stuff that comes specifically from this search term
    additionalFilters <- c(strsplit(input$wordcloud_searchTerm, 
                                    # dont include the words themselves
                                    split = " ")[[1]],
                           singleChar,
                           numericals,
                           verbs)
    # merge into final filter list
    filterList <- c(filterList, additionalFilters, fill = T)
    without_stopwords <- data.table::as.data.table(MetaboShiny::getFilteredWordFreqency(lcl$tables$wordcloud_orig, filterList))
    without_stopwords <- without_stopwords[without_stopwords$word != ""]
    lcl$tables$wordcloud_filt <- without_stopwords
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



