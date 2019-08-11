observeEvent(input$clear_prematch,{
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS prematch_content")
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS prematch_mapper")
  RSQLite::dbExecute(conn, "VACUUM")
  mSet$metshiParams$prematched <<- FALSE
  RSQLite::dbDisconnect(conn)
})
  
observeEvent(input$prematch,{
  
  if(is.null(mSet)){
    print("Requires mSet")
    return(NULL)
  }
  
  if(length(lcl$vectors$db_prematch_list) > 0){ # go through selected databases
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
    RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS prematch_content")
    RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS prematch_mapper")
    RSQLite::dbExecute(conn, "CREATE TABLE prematch_mapper(query_mz decimal(30,13),
                                                           structure TEXT,
                                                           `%iso` decimal(30,13),
                                                           adduct TEXT,
                                                           dppm decimal(30,13))")    
    RSQLite::dbExecute(conn, "CREATE TABLE prematch_content(name TEXT,
                                                            baseformula TEXT,
                                                            identifier TEXT,
                                                            description VARCHAR(255),
                                                            structure TEXT,
                                                            source TEXT)")
    RSQLite::dbExecute(conn, "CREATE INDEX map_mz ON prematch_mapper(query_mz)")
    RSQLite::dbExecute(conn, "CREATE INDEX cont_struc ON prematch_content(structure)")
    
    blocksize=100
    blocks = split(colnames(mSet$dataSet$norm), ceiling(seq_along(1:ncol(mSet$dataSet$norm))/blocksize))
    withProgress({
      i = 0
      matches = pbapply::pblapply(blocks, function(mzs){
        res = MetaDBparse::searchMZ(mzs = mzs,
                              ionmodes = getIonMode(mzs, lcl$paths$patdb),
                              base.dbname = gsub(x=gsub(basename(unlist(lcl$vectors$db_prematch_list)), pattern="\\.db", replacement=""),
                                                 pattern="\\.db",
                                                 replacement = "", perl=T),
                              ppm = as.numeric(mSet$ppm),
                              append = F,
                              outfolder = normalizePath(lcl$paths$db_dir))
        i <<- i + 1
        setProgress(value = i)
        list(mapper = unique(res[,c("query_mz", "structure", "%iso", "adduct", "dppm")]), 
             content = unique(res[,-c("query_mz", "%iso", "adduct", "dppm")]))
      })
    }, min=0, max=length(blocks))
    
    RSQLite::dbWriteTable(conn, "prematch_mapper", unique(data.table::rbindlist(lapply(matches, function(x) x$mapper))), append=T)
    RSQLite::dbWriteTable(conn, "prematch_content", unique(data.table::rbindlist(lapply(matches, function(x) x$content))), append=T)
    
    mSet$metshiParams$prematched<<-T
    
    RSQLite::dbDisconnect(conn)
  }
})

# triggers on clicking the 'search' button in sidebar
observeEvent(input$search_mz, {
  if(length(lcl$vectors$db_search_list) > 0){ # go through selected databases
      # get ion modes
    withProgress({
      lcl$tables$last_matches <<- MetaDBparse::searchMZ(mzs = lcl$curr_mz, 
                                                        ionmodes = getIonMode(lcl$curr_mz, 
                                                                              lcl$paths$patdb),
                                                        base.dbname = gsub(basename(unlist(lcl$vectors$db_search_list)), 
                                                                           pattern="\\.db", 
                                                                           replacement=""),
                                                        ppm=as.numeric(mSet$ppm),
                                                        append = F, 
                                                        outfolder=normalizePath(lcl$paths$db_dir))  
    })
    
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
