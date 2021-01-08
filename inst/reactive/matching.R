output$manual_search <- renderUI({
  if(search_button$on){
    tags$button(
      id = "search_mz",
      class = "btn btn-default action-button",
      img(src = "detective.png",
          height = "50px")
    )
  }else{
    fluidRow(align="center", 
             shiny::img(src = "pawprint.png",height = "50px"),
             br(),
             tags$h2("pre-matched")
    )
  }
})

shiny::observeEvent(input$classify_matches, {
  library(classyfireR)
  smi = shown_matches$forward_unique$structure
  names(smi) = smi
  rand = gsub("file","",basename(tempfile()))
  cfrid = gsubfn::fn$paste('metaboshiny_$rand')
  classif = classyfireR::submit_query(label = cfrid, input = smi, type = 'STRUCTURE')
  classif@classification
})

shiny::observeEvent(input$clear_prematch,{
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
  RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_mz")
  RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_struc")
  RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS cont_struc")
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_content")
  RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS match_mapper")
  RSQLite::dbExecute(conn, "VACUUM")
  mSet$metshiParams$prematched <<- FALSE
  # remove ALL mSet storage that had prematched m/z
  #mSet$storage <<- mSet$storage[!grepl(pattern = "\\(prematched m/z only\\)", names(mSet$storage))]
  search_button$go <- search_button$on <- TRUE
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
    
    db_list <- switch(search_type,
                      prematch = lcl$vectors$db_prematch_list,
                      search_mz = lcl$vectors$db_search_list)
    
    if(continue){ # go through selected databases
      
      conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_mz")
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS map_ba")
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS cont_ba")
      RSQLite::dbExecute(conn, "DROP INDEX IF EXISTS cont_str")
      
      blocksize=100
      blocks = switch(search_type,
                      prematch = split(colnames(mSet$dataSet$norm), ceiling(seq_along(1:ncol(mSet$dataSet$norm))/blocksize)),
                      search_mz = list(my_selection$mz))
      
      shiny::withProgress({
        i = 0
        matches = pbapply::pblapply(blocks, function(mzs){
          
           i = i + 1
          
          try({
            shiny::setProgress(i)
          })
          
          if(any(grepl(mzs, pattern="/"))){
            eachPPM = T
            ppm = ceiling(as.numeric(gsub(mzs, pattern="^.*/", replacement="")))
            mzs = gsub(mzs, pattern="/.*$", replacement="")
          }else{
            eachPPM = F
            ppm = mSet$ppm
            mzs = gsub(mzs, pattern="RT.*$", replacement="")
          }
          
          ionmode = sapply(mzs, function(mz) if(grepl(mz, pattern="\\-")) "negative" else "positive")
          mzs = stringr::str_match(mzs, "(^\\d+\\.\\d+)")[,2]
          
          if("cmmmediator" %in% db_list){
            ionmode = sapply(mzs, function(mz) if(grepl(mz, pattern="\\-")) "negative" else "positive")
            
            res.online <- MetaDBparse::searchMZonline(mz = mzs, 
                                                      mode = ionmode,
                                                      ppm = ppm,
                                                      which_db = "cmmr",
                                                      apikey = lcl$apikey)
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
                  row$adduct <- add
                  adj = MetaDBparse::doAdduct(structure = "", 
                                              formula = row$baseformula, 
                                              charge = 0, 
                                              query_adduct = add, 
                                              adduct_table = adducts)
                  row$finalcharge = adj$final.charge
                  row$fullformula = adj$final
                })
                row
              }), fill=T)
              colnames(res.online)[colnames(res.online) == "perciso"] <- "%iso"
              res.online$source <- c("cmmmediator")  
            }else{
              res.online = data.table::data.table()  
            }
          }else{
            res.online = data.table::data.table()
          }
          
          pred_dbs = intersect(c("magicball",
                                 "pubchem",
                                 "chemspider",
                                 "knapsack",
                                 "supernatural2",
                                 "chemidplus"), db_list)
          predict = length(pred_dbs) > 0
          
          if(predict){
            
            res.rows.predict = pbapply::pblapply(1:length(mzs), function(i){
              
              mz = mzs[i]
              
              if(eachPPM){
                ppm <- ppm[i]
              }else{
                ppm <- as.numeric(mSet$ppm)
              }
              
              res.predict = MetaDBparse::getPredicted(mz = as.numeric(mz),
                                                      ppm = ppm,
                                                      mode = ionmode,
                                                      rules = input$predict_rules,
                                                      elements = input$predict_elements)
              if(length(pred_dbs) == 1){
                if(pred_dbs == "magicball"){
                  search_db = F
                }else{
                  search_db = T
                }
              }else{
                search_db = T
              }
              
              if(nrow(res.predict) > 0 & search_db){
                if(lcl$apikey == " ") shiny::showNotification("Skipping ChemSpider, you haven't entered an API key!")
                res.big.db = MetaDBparse::searchFormulaWeb(unique(res.predict$baseformula),
                                                           search = intersect(db_list, 
                                                                              if(lcl$apikey != " ") c("pubchem", "chemspider", "knapsack","supernatural2","chemidplus") else 
                                                                                                       c("pubchem", "knapsack", "supernatural2","chemidplus")),
                                                           detailed = input$predict_details,
                                                           apikey = lcl$apikey)
                
                form_add_only <- res.predict[,c("baseformula",
                                                "adduct","query_mz",
                                                "%iso","dppm",
                                                "fullformula",
                                                "finalcharge")]
                keep = form_add_only$baseformula %in% res.big.db$baseformula
                form_add_only <- form_add_only[keep,]
                
                if(nrow(res.big.db) > 0){
                  results_full <- merge(res.big.db, 
                                        form_add_only, 
                                        on = "baseformula", 
                                        all.y = ifelse("magicball" %in% db_list, TRUE, FALSE), 
                                        allow.cartesian=T)  
                }else{
                  results_full <- res.predict
                }
                
                if(!("source" %in% colnames(results_full))){
                  results_full$source = "magicball"
                }
                results_full[is.na(source),]$compoundname <- results_full[is.na(source),]$baseformula
                results_full[is.na(source),]$source <- c("magicball")
                withSmi = which(results_full$structure != "")
                if(length(withSmi) > 0){
                  results_nosmi <- results_full[ -withSmi ]
                  results_nosmi$structure = paste0("[",results_nosmi$fullformula,"]0")
                  results_withsmi <- results_full[ withSmi ]  
                  if(input$predict_structure_check){
                    mols = MetaDBparse::smiles.to.iatom(results_withsmi$structure)
                    new.smi = MetaDBparse::iatom.to.smiles(mols)
                    results_withsmi$structure <- new.smi
                  }
                  if(input$predict_adduct_rules){
                    rulematch = MetaDBparse::countAdductRuleMatches(mols,
                                                                    adduct_rules)
                    structure.adducts.possible =  MetaDBparse::checkAdductRule(rulematch,
                                                                               adducts)
                    keep <- sapply(1:nrow(results_withsmi), function(i){
                      adduct = results_withsmi[i, "adduct"][[1]]
                      if(!is.na(adduct)){
                        structure.adducts.possible[i, ..adduct][[1]]
                      }
                    })
                    results_withsmi <- results_withsmi[keep,]
                  }
                  res = rbind(results_nosmi, results_withsmi)
                }else{
                  res = results_full
                }
                return(unique(res))
              }else{
                return(res.predict)
              }  
            })
            res.predict = unique(data.table::rbindlist(res.rows.predict, fill=T))
            }else{
            res.predict = data.table::data.table()
          }
          dbs.local = setdiff(gsub(x=gsub(basename(unlist(db_list)), 
                                          pattern="\\.db", 
                                          replacement=""),
                                   pattern="\\.db",
                                   replacement = "", 
                                   perl=T), gbl$vectors$db_no_build)
          if(length(dbs.local)>0){
            res.local = MetaDBparse::searchMZ(mzs = mzs,
                                              ionmodes = ionmode,
                                              base.dbname = dbs.local,
                                              ppm = ppm,
                                              append = F,
                                              outfolder = normalizePath(lcl$paths$db_dir))
          }else{
            res.local <- data.table::data.table()
          }
          res <- data.table::rbindlist(list(res.local, res.online, res.predict), use.names = T, fill=T)
          if(nrow(res) > 0){
            res[, c("query_mz") := lapply(.SD, as.character), .SDcols="query_mz"]
            ionmapper <- rep("", nrow(res))
            isneg <- grepl(res$adduct, pattern = "\\]\\d\\-")
            ionmapper[isneg] <- "-"
            res$query_mz <- paste0(res$query_mz, ionmapper)
            isocols = c("n2H", "n13C", "n15N")
            if("n2H" %in% colnames(res)){
              extracols=isocols
            }else{
              extracols=c()
            }
            getCols = c("query_mz",
                        "compoundname",
                        "baseformula",
                        "adduct",
                        "fullformula",
                        "finalcharge",
                        "identifier",
                        "description",
                        "structure",
                        "source",
                        extracols)
            
            list(mapper = unique(res[,c("query_mz", 
                                        "baseformula",
                                        "adduct", 
                                        "%iso",
                                        "dppm")]),
                 content = unique(res[, ..getCols])) 
          }else{
            list(mapper = data.table::data.table(),
                 content = data.table::data.table())
          }
        })
      }, min=0, max=length(blocks))
      
      mapper = unique(data.table::rbindlist(lapply(matches, function(x) x$mapper)))
      content = unique(data.table::rbindlist(lapply(matches, function(x) x$content)))
      
      if(nrow(mapper)>0){
        
        RSQLite::dbWriteTable(conn, 
                              "match_mapper", 
                              mapper, overwrite=T, use.names = T)
        RSQLite::dbWriteTable(conn, 
                              "match_content", 
                              content, overwrite=T, use.names = T)
        RSQLite::dbExecute(conn, "CREATE INDEX map_mz ON match_mapper(query_mz)")
        RSQLite::dbExecute(conn, "CREATE INDEX map_ba ON match_mapper(baseformula, adduct)")
        RSQLite::dbExecute(conn, "CREATE INDEX cont_ba ON match_content(baseformula, adduct)")
        RSQLite::dbExecute(conn, "CREATE INDEX cont_str ON match_content(structure)")
        if(search_type == "prematch"){
          mSet$metshiParams$prematched<<-T
          filemanager$do <- "save"
          search_button$on <- FALSE   
        }else{
          search$go <- TRUE
        }
        RSQLite::dbDisconnect(conn)
      }else{
        shiny::showNotification("No matches found!")  
        shown_matches$forward_unique <<- data.table::data.table()
        shown_matches$forward_full <<- data.table::data.table()
      }
      
    }else{
      if(length(grep(pattern = "supernatural|pubchem|magicball|chemspider|cmmediator|knapsack|chemidplus", 
                     x = lcl$vectors$built_dbs, value = T, invert = T)) == 0){
        MetaboShiny::metshiAlert("Please build at least one database (or select an online one) to enable this feature!")
      }else{
        MetaboShiny::metshiAlert("Please select at least one database to enable this feature!")
      }
      shinyBS::updateCollapse(session = session, 
                              id = "search_panel",
                              open = "databases",
                              style = "info")
      shiny::updateTabsetPanel(session, 
                               "tab_iden_1", 
                               selected = "pick_databases")
    }
  })
})

# triggers if isotope scoring is clicked after finding db matches
shiny::observeEvent(input$score_iso, {
  # check if the matches table even exists
  if(!data.table::is.data.table(shown_matches$forward_unique)) return(NULL)
  # check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
  
  if("isoScore" %in% colnames(shown_matches$forward_unique)){
    shown_matches$forward_unique <<- shown_matches$forward_unique[,-"isoScore"]
  }

  # get table including isotope scores
  # as input, takes user method for doing this scoring
  shiny::withProgress({
    score_table <- score.isos(qmz=my_selection$mz,
                              table = shown_matches$forward_unique, 
                              mSet = mSet,
                              ppm = as.numeric(mSet$ppm),
                              dbdir = lcl$paths$db_dir,
                              method = input$iso_score_method,
                              inshiny = F,
                              intprec = as.numeric(input$int_prec),
                              rtmode = input$iso_use_rt,
                              rtperc = input$iso_rt_perc,
                              useint = input$iso_use_int)
    colnames(score_table)[ colnames(score_table) == "score"] <- "isoScore"
    })
  shown_matches$forward_unique <- shown_matches$forward_unique[unique(score_table), on = c("fullformula")]
})

shiny::observeEvent(input$score_add, {
  # check if the matches table even exists
  if(!data.table::is.data.table(shown_matches$forward_unique)) return(NULL)
  # check if a previous scoring was already done (remove that column if so, new score is generated in a bit)
  if("addScore" %in% colnames(shown_matches$forward_unique)){
    shown_matches$forward_unique <<- shown_matches$forward_unique[,-"addScore"]
  }
  
  # get table including isotope scores
  # as input, takes user method for doing this scoring
  shiny::withProgress({
    score_table <- score.add(qmz = my_selection$mz,
                            table = shown_matches$forward_unique, 
                            adducts_considered=input$score_adducts,
                            mSet = mSet,
                            rtperc = input$add_rt_perc,
                            rtmode=input$add_use_rt,
                            mzppm = mSet$ppm,
                            dbdir = lcl$paths$db_dir,
                            inshiny = T)
  })
  colnames(score_table)[ colnames(score_table) == "score"] <- "addScore"
  shown_matches$forward_unique <- shown_matches$forward_unique[score_table, on = c("structure")]
})

shiny::observeEvent(input$search_pubmed,{
  statsmanager$calculate <- "match_wordcloud_pm"
  datamanager$reload <- "match_wordcloud_pm"
})

shiny::observe({
  # - - filters - -
  if(search$go){
    shiny::withProgress({
      if(input$tab_iden_2 == "mzmol"){
        if(lcl$prev_mz != my_selection$mz & !identical(lcl$vectors$prev_dbs, lcl$vectors$db_search_list)){
          matches = data.table::as.data.table(get_prematches(who = gsub(my_selection$mz, 
                                                                        pattern="/.*$|RT.*$", 
                                                                        replacement=""),
                                                             what = "map.query_mz",
                                                             patdb = lcl$paths$patdb,
                                                             showadd = c(),
                                                             showdb = c(),
                                                             showiso = c(),
                                                             showIsolabels = input$show_iso_labels))
          if(nrow(matches) == 0){
            shown_matches$forward_unique <- data.table::data.table()
            shown_matches$forward_full <- data.table::data.table()
            return(NULL)
          }
          lcl$prev_mz <<- my_selection$mz
          lcl$vectors$prev_dbs <<- lcl$vectors$db_search_list
          pieinfo$db <- reshape::melt(table(matches$source))
          pieinfo$add <- reshape::melt(table(matches$adduct))
          pieinfo$iso <- reshape::melt(table(matches$isocat))
        }
        
        mzMode = if(grepl(my_selection$mz, pattern="\\-")) "negative" else "positive"
        
        matches = data.table::as.data.table(get_prematches(who = gsub(my_selection$mz, pattern="/.*$|RT.*$", replacement=""),
                                                           what = "map.query_mz",
                                                           patdb = lcl$paths$patdb,
                                                           showadd = result_filters$add[[mzMode]],
                                                           showdb = result_filters$db,
                                                           showiso = result_filters$iso,
                                                           showIsolabels = input$show_iso_labels))  
        
        if(nrow(matches)>0){
          
          shiny::setProgress(0.2)
          
          matches$compoundname[matches$source != "magicball"] <- tolower(matches$compoundname[matches$source != "magicball"])
          
          # =====
          
          uniques = data.table::as.data.table(unique(data.table::as.data.table(matches)[, -c("source", 
                                                                                             "description",
                                                                                             "identifier"),
                                                                                        with=F]))
          shiny::setProgress(0.6)
          
          # === aggregate ===
          
          info_only = unique(matches[,c("compoundname", 
                                        "source", 
                                        "structure", 
                                        "description",
                                        "identifier"),with=F])
          info_only$description <- paste0("Database ID: ", info_only$identifier, ". ", info_only$description)
          info_only <- unique(info_only[,-"identifier"])
          
          info_no_na <- info_only[!is.na(info_only$structure)]
          info_na <- info_only[is.na(info_only$structure)]
          
          info_aggr <- aggregate(info_no_na, by = list(info_no_na$compoundname), FUN = function(x) paste0(x, collapse = "SEPERATOR"))
          info_aggr <- aggregate(info_aggr, by = list(info_aggr$structure), FUN = function(x) paste0(x, collapse = "SEPERATOR"))
          
          info_aggr <- rbind(info_aggr, info_na, fill=T)
          
          # fix structures
          split_structs <- strsplit(info_aggr$structure, split = "SEPERATOR")
          main_structs <- unlist(lapply(split_structs, function(x) x[[1]]))
          info_aggr$structure <- main_structs
          
          # move extra names to descriptions
          split_names <- strsplit(info_aggr$compoundname, split = "SEPERATOR")
          main_names <- unlist(lapply(split_names, function(x) if(length(x) > 1) x[[1]] else x[1]))
          
          synonyms <- unlist(lapply(split_names, function(x){
            if(length(x)>1){
              paste0("SYNONYMS: ", paste0(unique(x[2:length(x)]), collapse=", "),".SEPERATOR") 
            }else NA
          }))
          
          info_aggr$compoundname <- main_names
          has.syn <- which(!is.na(synonyms))
          
          info_aggr$description[has.syn] <- paste0(synonyms[has.syn], info_aggr$description[has.syn])
          info_aggr <- data.table::as.data.table(info_aggr)
          
          # =================
          
          uniques = uniques[structure %in% info_aggr$structure]
          
          is.no.na.uniq <- which(!is.na(uniques$structure))
          is.no.na.info <- which(!is.na(info_aggr$structure))
          
          uniq.to.aggr <- match(uniques[is.no.na.uniq]$structure, 
                                info_aggr[is.no.na.info]$structure)
          
          uniques$compoundname[is.no.na.uniq] <- info_aggr$compoundname[is.no.na.info][uniq.to.aggr]
          uniques$structure[is.no.na.uniq] <- info_aggr$structure[is.no.na.info][uniq.to.aggr]
          
          uniques <- unique(uniques)
          
          shown_matches$forward_unique <- uniques[,-grepl(colnames(uniques), pattern = "Group\\.\\d"),with=F]
          shown_matches$forward_full <- info_aggr[,-grepl(colnames(info_aggr), pattern = "Group\\.\\d"),with=F]
          
        }else{
          shown_matches$forward_unique <- data.table::data.table()
          shown_matches$forward_full <- data.table::data.table()
        }
        
        my_selection$name <- ""
        my_selection$form <- ""
        my_selection$struct <- ""  
        
      }else{
        NULL
      }
      
      shiny::setProgress(0.8)
      
    })
    search$go <- FALSE #reset self
  }
})
