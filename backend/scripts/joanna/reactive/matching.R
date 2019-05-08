# triggers on clicking the 'search' button in sidebar
observeEvent(input$search_cpd, {

  if(length(local$vectors$db_search_list) > 0){ # go through selected databases
    
    local$tables$last_matches <<- unique(multimatch(local$curr_mz, 
                                                     local$vectors$db_search_list,
                                                     inshiny = F,
                                                     search_pubchem = input$magicball_pubchem_cids,
                                                     pubchem_detailed = input$magicball_pubchem_details,
                                                     calc_adducts = local$vectors$add_list)) # match with all
    # - - -
    adduct_dist <- melt(table(local$tables$last_matches$adduct))
    db_dist <- melt(table(local$tables$last_matches$source))
    
    if(nrow(local$tables$last_matches) > 0){
      # render word cloud
      
      try({
        
        # remove unwanted words (defined in global) from description
        filtered_descriptions <- sapply(1:length(local$tables$last_matches$description), 
                                        function(i){
                                          # get description
                                          desc <- local$tables$last_matches$description[[i]]
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
        docs <- tm::tm_map(docs, tm::removeWords, global$vectors$wordcloud$skip) 
        
        # Remove whitespace
        docs <- tm::tm_map(docs, tm::stripWhitespace)
        
        # # Text stemming
        # docs <- tm_map(docs, stemDocument)
        
        doc_mat <- tm::TermDocumentMatrix(docs)
        
        m <- as.matrix(doc_mat)
        
        v <- sort(rowSums(m), decreasing = TRUE)
        
        local$tables$word_freq <<- data.frame(name = names(v), value = v)
        
        # - - return - - 
        
        res = head(local$tables$word_freq, 30)
        
        wordcloud_plot <- wordcloud2::wordcloud2(data = data.frame(word = res$name,
                                                                   freq = res$value), 
                                                 size = 0.7, 
                                                 shuffle = FALSE,
                                                 fontFamily = options$font4,
                                                 ellipticity = 1, 
                                                 minRotation = -pi/8, 
                                                 maxRotation = pi/8,
                                                 shape = 'heart')
        output$wordcloud_desc  <- wordcloud2::renderWordcloud2(wordcloud_plot)
        
        local$vectors$pie_add <<- adduct_dist$Var1
        
        output$match_pie_add <- plotly::renderPlotly({
          plot_ly(adduct_dist, labels = ~Var1, values = ~value, size=~value*10, type = 'pie',
                  textposition = 'inside',
                  textinfo = 'label+percent',
                  insidetextfont = list(color = '#FFFFFF'),
                  hoverinfo = 'text',
                  text = ~paste0(Var1, ": ", value, ' matches'),
                  marker = list(colors = colors,
                                line = list(color = '#FFFFFF', width = 1)),
                  #The 'pull' attribute can also be used to create space between the sectors
                  showlegend = FALSE) %>%
            layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        })
        
        local$vectors$pie_db <<- db_dist$Var1
        
        output$match_pie_db <- plotly::renderPlotly({
          plot_ly(db_dist, labels = ~Var1, values = ~value, size=~value*10, type = 'pie',
                  textposition = 'inside',
                  textinfo = 'label+percent',
                  insidetextfont = list(color = '#FFFFFF'),
                  hoverinfo = 'text',
                  text = ~paste0(Var1, ": ", value, ' matches'),
                  marker = list(colors = colors,
                                line = list(color = '#FFFFFF', width = 1)),
                  #The 'pull' attribute can also be used to create space between the sectors
                  showlegend = FALSE) %>%
            layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        })
      })
    }
    
    shown_matches$table <- 
    if(nrow(local$tables$last_matches)>0){
      local$tables$last_matches 
    }else{
      data.table('name' = "Didn't find anything ( •́ .̫ •̀ )")
    }
  }
})

# triggers if isotope scoring is clicked after finding db matches
observeEvent(input$score_iso, {
  
  # check if the matches table even exists
  if(!data.table::is.data.table(local$tables$last_matches)) return(NULL)
  
  # check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
  if("score" %in% colnames(local$tables$last_matches)){
    local$tables$last_matches <<- local$tables$last_matches[,-"score"]
  }
  
  intprec = as.numeric(input$int_prec)/100.00
  
  # get table including isotope scores
  # as input, takes user method for doing this scoring
  withProgress({
    score_table <- score.isos(local$paths$patdb, method=input$iso_score_method, inshiny=T, intprec = intprec)
    })
  
  # update the match table available to the rest of metaboshiny
  local$tables$last_matches <<- local$tables$last_matches[score_table, on = c("baseformula", "adduct")]
  
  shown_matches$table <<- local$tables$last_matches
 
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
                                 top = 30)
      
      setProgress(0.4)
      
      tbl <- data.frame(
        pmids = RISmed::PMID(fetch),
        titles = RISmed::ArticleTitle(fetch)
        #abstracts = RISmed::AbstractText(fetch)
      )
      
      setProgress(0.6)
      
      wordcloud_plot <- wordcloud2::wordcloud2(data = data.frame(word = res$stsp,
                                                                 freq = res$Freq), 
                                               size = 0.7, 
                                               shuffle = FALSE, 
                                               ellipticity = 1,
                                               fontFamily = options$font4,
                                               minRotation = -pi/8, 
                                               maxRotation = pi/8,
                                               shape = 'heart')
      
      setProgress(0.8)
      
      output$wordcloud_pubmed  <- wordcloud2::renderWordcloud2(wordcloud_plot)
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