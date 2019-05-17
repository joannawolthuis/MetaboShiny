# triggers on clicking the 'search' button in sidebar
observeEvent(input$search_mz, {

  if(length(lcl$vectors$db_search_list) > 0){ # go through selected databases

    lcl$tables$last_matches <<- unique(multimatch(lcl$curr_mz,
                                                    lcl$vectors$db_search_list,
                                                    inshiny = T,
                                                    search_pubchem = input$magicball_pubchem_cids,
                                                    pubchem_detailed = input$magicball_pubchem_details,
                                                    calc_adducts = lcl$vectors$add_list,
                                                    patdb = normalizePath(lcl$paths$patdb),
                                                    db_dir = normalizePath(lcl$paths$db_dir)
                                                    )) # match with all

    adduct_dist <- melt(table(lcl$tables$last_matches$adduct))
    db_dist <- melt(table(lcl$tables$last_matches$source))

    if(nrow(lcl$tables$last_matches) > 0){
      # render word cloud

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

        # Remove english common stopwords
        docs <- tm::tm_map(docs, tm::removeWords, tm::stopwords("english"))
        # Remove your own stop word
        # ADD stopwords as a character vector
        docs <- tm::tm_map(docs, tm::removeWords, gbl$vectors$wordcloud$skip)

        # Remove whitespace
        docs <- tm::tm_map(docs, tm::stripWhitespace)

        # # Text stemming
        # docs <- tm_map(docs, stemDocument)

        doc_mat <- tm::TermDocumentMatrix(docs)

        m <- as.matrix(doc_mat)

        v <- sort(rowSums(m), decreasing = TRUE)

        lcl$tables$word_freq <<- data.frame(name = names(v), value = v)

        # - - return - -

        lcl$vectors$pie_add <<- adduct_dist
        lcl$vectors$pie_db <<- db_dist
        
        datamanager$reload <<- "matchplots"
        
      })
    }

    shown_matches$table <- if(nrow(lcl$tables$last_matches) > 0){
      lcl$tables$last_matches
    }else{
      data.table('name' = "Didn't find anything ( •́ .̫ •̀ )")
    }
  }
})

# triggers if isotope scoring is clicked after finding db matches
observeEvent(input$score_iso, {

  # check if the matches table even exists
  if(!data.table::is.data.table(lcl$tables$last_matches)) return(NULL)

  # check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
  if("score" %in% colnames(lcl$tables$last_matches)){
    lcl$tables$last_matches <<- lcl$tables$last_matches[,-"score"]
  }

  intprec = as.numeric(input$int_prec)/100.00

  # get table including isotope scores
  # as input, takes user method for doing this scoring
  withProgress({
    score_table <- score.isos(mSet = mSet, lcl$paths$patdb, method=input$iso_score_method, inshiny=T, intprec = intprec)
    })

  # update the match table available to the rest of metaboshiny
  lcl$tables$last_matches <<- lcl$tables$last_matches[score_table, on = c("baseformula", "adduct")]

  shown_matches$table <<- lcl$tables$last_matches

})

observeEvent(input$search_pubmed, {

  withProgress({

    abstr <- RISmed::EUtilsSummary(
      as.character(input$pm_query),
      type="esearch",
      db="pubmed",
      datetype='pdat',
      mindate=input$pm_year[1],
      maxdate=input$pm_year[2],
      retmax=as.numeric(input$pm_max))

    setProgress(0.2)

    fetch <- RISmed::EUtilsGet(abstr)

    if(length(fetch@PMID) > 0){
      res <- abstracts2wordcloud(abstracts = fetch,
                                 top = input$wc_topn_pm)

      setProgress(0.4)

      tbl <- data.frame(
        pmids = RISmed::PMID(fetch),
        titles = RISmed::ArticleTitle(fetch)
        #abstracts = RISmed::AbstractText(fetch)
      )

      setProgress(0.6)

      wcdata <- res
      colnames(wcdata) <- c("word", "freq")
      
      output$wordcloud_desc_pm  <- wordcloud2::renderWordcloud2({
        wordcloud2::wordcloud2(wcdata,
                               size = 0.7,
                               shuffle = FALSE,
                               fontFamily = getOptions(lcl$paths$opt.loc)$font4,
                               ellipticity = 1,
                               minRotation = -pi/8,
                               maxRotation = pi/8,
                               shape = 'heart')
      })
      
      output$wordbar_desc_pm <- plotly::renderPlotly({
        ggPlotWordBar(wcdata = wcdata,
                      cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                      plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                      plotlyfy = TRUE,
                      font = lcl$aes$font)})
      
      setProgress(0.8)
    }else{
      tbl <- data.table("no papers found" = "Please try another term!	(｡•́︿•̀｡)")
    }

    output$pm_tab <- DT::renderDataTable({
      DT::datatable(tbl,
                    selection = "single",
                    options = list(lengthMenu = c(5, 10, 15),
                                   pageLength = 5)
      )
    })

  })

})
