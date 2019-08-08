observeEvent(input$prematch,{
  print("trigger")
  if(length(lcl$vectors$db_prematch_list) > 0){ # go through selected databases
    blocksize=100
    blocks = split(colnames(mSet$dataSet$norm), ceiling(seq_along(1:ncol(mSet$dataSet$norm))/blocksize))
    withProgress({
      i = 0
      match_rows = lapply(blocks, function(mzs){
        MetaDBparse::searchMZ(mzs = mzs,
                              ionmodes = getIonMode(mzs, lcl$paths$patdb),
                              base.dbname = gsub(x=basename(unlist(lcl$vectors$db_prematch_list)),
                                                 pattern="\\.db",
                                                 replacement = "", perl=T),
                              ppm = 3,
                              append = F,
                              outfolder = normalizePath(lcl$paths$db_dir))
        i = i + 1
        setProgress(value = i/length(blocks))
      })
    }, min=0, max=length(blocks))
    lcl$tables$pre_matches <- unique(data.table::rbindlist(match_rows))
  }
})

# triggers on clicking the 'search' button in sidebar
observeEvent(input$search_mz, {

  if(length(lcl$vectors$db_search_list) > 0){ # go through selected databases

    match_rows <- pbapply::pblapply(lcl$vectors$db_search_list, function(db){
      # get ion modes
      
      # - - - - - - - 
      matches = MetaDBparse::searchMZ(mzs = lcl$curr_mz, 
                                      ionmodes = getIonMode(lcl$curr_mz, lcl$paths$patdb),
                                      base.dbname = gsub(basename(db), pattern="\\.db", replacement=""),
                                      ppm=3,
                                      append = F, 
                                      outfolder=normalizePath(lcl$paths$db_dir))
    })
    
    lcl$tables$last_matches <- unique(data.table::rbindlist(match_rows))

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
  if(!data.table::is.data.table(shown_matches$table)) return(NULL)

  # check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
  if("score" %in% colnames(shown_matches$table)){
    shown_matches$table <<- shown_matches$table[,-"score"]
  }

  intprec = as.numeric(input$int_prec)/100.00

  # get table including isotope scores
  # as input, takes user method for doing this scoring
  withProgress({
    score_table <- score.isos(table = shown_matches$table, mSet = mSet, lcl$paths$patdb, method=input$iso_score_method, inshiny=T, intprec = intprec)
    })

  # update the match table available to the rest of metaboshiny
  shown_matches$table <<- shown_matches$table[score_table, on = c("baseformula", "adduct")]
})

observeEvent(input$search_pubmed, {

  withProgress({

    #input <- list(pm_query = "glucose", pm_year = c(2010, 2019), pm_max = 300)
    
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
                                 top = input$wc_topn_pm, 
                                 queryword = input$pm_query)

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
