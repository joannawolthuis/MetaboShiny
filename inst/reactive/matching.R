shiny::observeEvent(input$clear_prematch,{
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
  
  RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_mz")
  RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_struc")
  RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS cont_struc")
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_content")
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_mapper")
  RSQLite::dbExecute(conn, "VACUUM")
  mSet$metshiParams$prematched <<- FALSE
  search_button$go <- TRUE
  RSQLite::dbDisconnect(conn)
})
  
shiny::observeEvent(input$prematch,{
  if(is.null(mSet)){
    MetaboShiny::metshiAlert("Requires mSet!")
    return(NULL)
  }
  if(length(lcl$vectors$db_prematch_list) > 0){ # go through selected databases
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
    RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_mz")
    RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_struc")
    RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS cont_struc")
    RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_content")
    RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_mapper")
    RSQLite::dbExecute(conn, "CREATE TABLE IF NOT EXISTS match_mapper(query_mz decimal(30,13),
                                                                       structure TEXT,
                                                                       `%iso` decimal(30,13),
                                                                       adduct TEXT,
                                                                       dppm decimal(30,13))")    
    RSQLite::dbExecute(conn, "CREATE TABLE IF NOT EXISTS match_content(query_mz decimal(30,13),
                                                                      name TEXT,
                                                                      baseformula TEXT,
                                                                      fullformula TEXT,
                                                                      finalcharge INT,
                                                                      identifier TEXT,
                                                                      description VARCHAR(255),
                                                                      structure TEXT,
                                                                      source TEXT)")
    
    blocksize=100
    blocks = split(colnames(mSet$dataSet$norm), ceiling(seq_along(1:ncol(mSet$dataSet$norm))/blocksize))
    shiny::withProgress({
      i = 0
      matches = pbapply::pblapply(blocks, function(mzs){
        res = MetaDBparse::searchMZ(mzs = mzs,
                              ionmodes = getIonMode(mzs, lcl$paths$patdb),
                              base.dbname = gsub(x=gsub(basename(unlist(lcl$vectors$db_prematch_list)), 
                                                        pattern="\\.db", 
                                                        replacement=""),
                                                 pattern="\\.db",
                                                 replacement = "", 
                                                 perl=T),
                              ppm = as.numeric(mSet$ppm),
                              append = F,
                              outfolder = normalizePath(lcl$paths$db_dir))
        i <<- i + 1
        shiny::setProgress(value = i)
        list(mapper = unique(res[,c("query_mz", "structure", 
                                    "%iso", "adduct", "dppm")]),
             content = unique(res[,-c("%iso", "adduct", "dppm")]))
      })
    }, min=0, max=length(blocks))
    
    RSQLite::dbWriteTable(conn, 
                          "match_mapper", 
                          unique(data.table::rbindlist(lapply(matches, function(x) x$mapper))), append=T)
    RSQLite::dbWriteTable(conn, 
                          "match_content", 
                          unique(data.table::rbindlist(lapply(matches, function(x) x$content))), append=T)
    
    RSQLite::dbExecute(conn, "CREATE INDEX map_mz ON match_mapper(query_mz)")
    RSQLite::dbExecute(conn, "CREATE INDEX map_struc ON match_mapper(structure)")
    RSQLite::dbExecute(conn, "CREATE INDEX cont_struc ON match_content(structure)")
    
    mSet$metshiParams$prematched<<-T
    search_button$on <<- FALSE
    
    RSQLite::dbDisconnect(conn)
  }else{
    MetaboShiny::metshiAlert("Please build at least one database to enable this feature!")
    return(NULL)
  }
})

# triggers on clicking the 'search' button in sidebar
shiny::observeEvent(input$search_mz, {
  if(length(lcl$vectors$db_search_list) > 0 & my_selection$mz != ""){ # go through selected databases
      # get ion modes
    shiny::withProgress({
      conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_mz")
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_struc")
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS cont_struc")
      RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_mapper") 
      RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_content")   
      RSQLite::dbExecute(conn, "CREATE TABLE IF NOT EXISTS match_mapper(query_mz decimal(30,13),
                                                                       structure TEXT,
                                                                       `%iso` decimal(30,13),
                                                                       adduct TEXT,
                                                                       dppm decimal(30,13))")    
      RSQLite::dbExecute(conn, "CREATE TABLE IF NOT EXISTS match_content(query_mz decimal(30,13),
                                                                        name TEXT,
                                                                        baseformula TEXT,
                                                                        fullformula TEXT,
                                                                        finalcharge INT,
                                                                        identifier TEXT,
                                                                        description VARCHAR(255),
                                                                        structure TEXT,
                                                                        source TEXT)")
      
      res <- MetaDBparse::searchMZ(mzs = my_selection$mz, 
                                   ionmodes = getIonMode(my_selection$mz, 
                                                         lcl$paths$patdb),
                                   base.dbname = gsub(basename(unlist(lcl$vectors$db_search_list)), 
                                                      pattern="\\.db", 
                                                      replacement=""),
                                   ppm=as.numeric(mSet$ppm),
                                   append = F, 
                                   outfolder=normalizePath(lcl$paths$db_dir))
      
      if(nrow(res)>0){
        mapper = unique(res[,c("query_mz", "structure", 
                               "%iso", "adduct", "dppm")]) 
        content = unique(res[,-c("%iso", "adduct", "dppm")])
        RSQLite::dbWriteTable(conn, "match_mapper", mapper, append=T)
        RSQLite::dbWriteTable(conn, "match_content", content, append=T)
        RSQLite::dbExecute(conn, "CREATE INDEX map_mz ON match_mapper(query_mz)")
        RSQLite::dbExecute(conn, "CREATE INDEX map_struc ON match_mapper(structure)")
        RSQLite::dbExecute(conn, "CREATE INDEX cont_struc ON match_content(structure)")
        RSQLite::dbDisconnect(conn)
        search$go <<- TRUE
      }
    })
    
  }
})

# triggers if isotope scoring is clicked after finding db matches
shiny::observeEvent(input$score_iso, {

  # check if the matches table even exists
  if(!data.table::is.data.table(shown_matches$forward_unique)) return(NULL)

  # check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
  if("score" %in% colnames(shown_matches$forward)){
    shown_matches$forward_unique <<- shown_matches$forward_unique[,-"score"]
  }

  intprec = as.numeric(input$int_prec)/100.00

  # get table including isotope scores
  # as input, takes user method for doing this scoring
  shiny::withProgress({
    score_table <- MetaboShiny::score.isos(table = shown_matches$forward_unique, 
                              mSet = mSet, 
                              ppm = as.numeric(mSet$ppm),
                              patdb = lcl$paths$patdb, 
                              method = input$iso_score_method, 
                              inshiny = T, intprec = intprec)
    })
  shown_matches$forward_unique <- shown_matches$forward_unique[score_table, on = c("fullformula")]
})

shiny::observeEvent(input$search_pubmed,{
  statsmanager$calculate <- "match_wordcloud_pm"
  datamanager$reload <- "match_wordcloud_pm"
})
