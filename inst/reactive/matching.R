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
  
lapply(c("prematch","search_mz"), function(search_type){
  shiny::observeEvent(input[[search_type]],{
    if(is.null(mSet)){
      MetaboShiny::metshiAlert("Requires mSet!")
      return(NULL)
    }
    
    continue = switch(search_type, 
                      prematch = length(lcl$vectors$db_prematch_list) > 0,
                      search_mz = length(lcl$vectors$db_search_list) > 0 & my_selection$mz != "")
    
    db_list = switch(search_type,
                     prematch = lcl$vectors$db_prematch_list,
                     search_mz = lcl$vectors$db_search_list)
    
    if(continue){ # go through selected databases
      
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
      blocks = switch(search_type,
                      prematch = split(colnames(mSet$dataSet$norm), ceiling(seq_along(1:ncol(mSet$dataSet$norm))/blocksize)),
                      search_mz = list(my_selection$mz))
      shiny::withProgress({
        i = 0
        matches = pbapply::pblapply(blocks, function(mzs){
          if("cmmmediator" %in% db_list){
            res.online <- MetaDBparse::searchMZonline(mz = mzs, 
                                                      mode = MetaboShiny::getIonMode(mzs, 
                                                                        lcl$paths$patdb),
                                                      ppm = as.numeric(mSet$ppm),
                                                      which_db = "cmmr")
            if(nrow(res.online) > 0 ){
              if(!("structure" %in% colnames(res.online))){
                res.online$structure <- c("")
              }
              res.online <- data.table::rbindlist(lapply(1:nrow(res.online), function(i){
                row = res.online[i,]
                row$finalcharge = "?"
                row$fullformula = "?"
                try({
                  add = grep(adducts$Name, pattern = paste0("\\[\\Q", row$adduct, "\\E\\]"), value=T)
                  adj = MetaDBparse::doAdduct(structure = "", 
                                              formula = row$baseformula, 
                                              charge = 0, 
                                              query_adduct = add, 
                                              adduct_table = adducts)
                  row$finalcharge = adj$final.charge
                  row$fullformula = adj$final  
                })
                row
              }))
              colnames(res.online)[colnames(res.online) == "perciso"] <- "%iso"
              res.online$source <- c("cmmmediator")  
            }else{
              res.online = data.table::data.table()  
            }
          }else{
            res.online = data.table::data.table()
          }
          predict = length(intersect(c("magicball",
                                       "pubchem",
                                       "chemspider"), db_list)) > 0
          if(predict){
            res.rows.predict = pbapply::pblapply(mzs, function(mz){
              res.predict = MetaDBparse::getPredicted(mz = as.numeric(mz), 
                                                      ppm = as.numeric(mSet$ppm),
                                                      mode = MetaboShiny::getIonMode(mz, 
                                                                                     lcl$paths$patdb))
              if(("pubchem" %in% db_list | "chemspider" %in% db_list) & nrow(res.predict) > 0){
                if(lcl$apikey == "") shiny::shinyNotification("Skipping ChemSpider, you haven't entered an API key!")
                res.big.db = MetaDBparse::searchFormulaWeb(unique(res.predict$baseformula),
                                                           search = intersect(db_list, 
                                                                              if(lcl$apikey != "") c("pubchem","chemspider") else "pubchem"),
                                                           detailed = input$predict_details,
                                                           apikey = lcl$apikey)
                
                res.big.db$`%iso` <- c(100)
                form_add_only <- res.predict[,c("baseformula",
                                                "adduct")]
                results_full <- merge(res.big.db, form_add_only)
                withSmi = which(results_full$structure != "")
                results_nosmi <- results_full[ -withSmi ]
                results_withsmi <- results_full[ withSmi ]
                
                if(input$predict_structure_check){
                  mols = MetaDBparse::smiles.to.iatom(results_withsmi$structure)
                  new.smi = MetaDBparse::iatom.to.smiles(mols)
                  results_withsmi$structure <- new.smi
                }
                if(input$predict_adduct_rules){
                  rulematch = MetaDBparse::countAdductRuleMatches(mols, adduct_rules)
                  structure.adducts.possible =  MetaDBparse::checkAdductRule(rulematch,
                                                                             adduct_table)
                  keep <- sapply(1:nrow(results_withsmi), function(i){
                    adduct = results_withsmi[i, "adduct"][[1]]
                    if(!is.na(adduct)){
                      structure.adducts.possible[i, ..adduct][[1]]
                    }
                  })
                  results_withsmi <- results_withsmi[keep,]
                }
                return(rbind(results_nosmi, results_withsmi))
              }else{
                return(data.table::data.table())
              }  
            })
            res.predict = data.table::rbindlist(res.rows.predict)
          }else{
            res.predict = data.table::data.table()
          }
          dbs.local = setdiff(gsub(x=gsub(basename(unlist(db_list)), 
                                          pattern="\\.db", 
                                          replacement=""),
                                   pattern="\\.db",
                                   replacement = "", 
                                   perl=T), gbl$vectors$db_no_build)
          if(length(dbs.local)>1){
            res.local = MetaDBparse::searchMZ(mzs = mzs,
                                              ionmodes = MetaboShiny::getIonMode(mzs, lcl$paths$patdb),
                                              base.dbname = dbs.local,
                                              ppm = as.numeric(mSet$ppm),
                                              append = F,
                                              outfolder = normalizePath(lcl$paths$db_dir))
          }else{
            res.local <- data.table::data.table()
          }
          res <- data.table::rbindlist(list(res.local, res.online, res.predict),use.names = T)
          if(nrow(res) > 0){
            list(mapper = unique(res[,c("query_mz", 
                                        "structure", 
                                        "%iso", 
                                        "adduct", 
                                        "dppm")]),
                 content = unique(res[,c("query_mz",
                                         "name",
                                         "baseformula",
                                         "fullformula",
                                         "finalcharge",
                                         "identifier",
                                         "description",
                                         "structure",
                                         "source")]))  
          }else{
            list(mapper = data.table::data.table(),
                 content = data.table::data.table())
          }
          
        })
      }, min=0, max=length(blocks))
      
      RSQLite::dbWriteTable(conn, 
                            "match_mapper", 
                            unique(data.table::rbindlist(lapply(matches, function(x) x$mapper))), append=T,use.names = T)
      RSQLite::dbWriteTable(conn, 
                            "match_content", 
                            unique(data.table::rbindlist(lapply(matches, function(x) x$content))), append=T,use.names = T)
      
      RSQLite::dbExecute(conn, "CREATE INDEX map_mz ON match_mapper(query_mz)")
      RSQLite::dbExecute(conn, "CREATE INDEX map_struc ON match_mapper(structure)")
      RSQLite::dbExecute(conn, "CREATE INDEX cont_struc ON match_content(structure)")
      
      if(search_type == "prematch"){
        mSet$metshiParams$prematched<<-T
        search_button$on <- FALSE   
      }else{
        search$go <- TRUE
      }
      RSQLite::dbDisconnect(conn)
    }else{
      MetaboShiny::metshiAlert("Please build at least one database to enable this feature!")
      return(NULL)
    }
  })
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
