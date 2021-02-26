# triggers on clicking the 'browse database' function
observeEvent(input$browse_db,{
  # get all compounds in the selected databases
  cpd_list <- lapply(lcl$vectors$db_search_list, FUN=function(match.table){
    res = MetaDBparse::showAllBase(outfolder = lcl$paths$db_dir,
                                   base.dbname = gsub(basename(match.table), 
                                                      pattern = "\\.db",
                                                      replacement=""))
    res$source = match.table
    res
  })
  # join the individual result tables together
  browse_content$table <<- unique(data.table::as.data.table(data.table::rbindlist(cpd_list)))
  # render table for UI
})

# triggers on reverse searching TODO: fix this, it's broken
observeEvent(input$revsearch_mz, {
  curr_row <- input$browse_tab_rows_selected
  if(!is.null(curr_row)){
    my_selection$structure <- browse_content$table[curr_row,c('structure')][[1]]
  }
})

shiny::observe({
  if(my_selection$struct != "" & input$tab_iden_2 == "molmz"){
    if(FALSE){#!mSet$metshiParams$prematched){
      MetaboShiny::metshiAlert("Please perform pre-matching first to enable this feature!")
      return(NULL)
    }else{
      if(lcl$prev_struct != my_selection$struct & !is.null(input$browse_tab_rows_selected)){
      
        if(!mSet$metshiParams$prematched){
          
          # spaghetti
          rev_matches = metshiRevSearch(mSet, my_selection$struct, "extended", lcl$paths$db_dir)

          browse_info = browse_content$table[input$browse_tab_rows_selected, ]
          
          rev_matches = merge(browse_info,rev_matches, by="structure")
          
          # WRITE TO PREMATCH DB
          isocols = c("n2H", "n13C", "n15N")
          if("n2H" %in% colnames(rev_matches)){
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
          
          mapper = unique(rev_matches[,c("query_mz", 
                                      "baseformula",
                                      "adduct", 
                                      "%iso",
                                      "dppm")])
          
          content = unique(rev_matches[, ..getCols])
          
          conn <- RSQLite::dbConnect(RSQLite::SQLite(), lcl$paths$patdb) # change this to proper var later
          
          RSQLite::dbWriteTable(conn, 
                                "match_mapper", 
                                mapper, 
                                overwrite = !mSet$metshiParams$prematched, 
                                append = mSet$metshiParams$prematched, 
                                use.names = T)
          
          RSQLite::dbWriteTable(conn, 
                                "match_content", 
                                content, 
                                overwrite = !mSet$metshiParams$prematched,
                                append = mSet$metshiParams$prematched,
                                use.names = T)
          
          if(!mSet$metshiParams$prematched){
            RSQLite::dbExecute(conn, "CREATE INDEX map_mz ON match_mapper(query_mz)")
            RSQLite::dbExecute(conn, "CREATE INDEX map_ba ON match_mapper(baseformula, adduct)")
            RSQLite::dbExecute(conn, "CREATE INDEX cont_ba ON match_content(baseformula, adduct)")
            RSQLite::dbExecute(conn, "CREATE INDEX cont_str ON match_content(structure)")  
          }
          RSQLite::dbDisconnect(conn)
          # ===================
        }
        rev_matches = get_prematches(who = my_selection$struct,
                                     what = "con.structure",
                                     patdb = lcl$paths$patdb,
                                     showadd = c(),
                                     showiso = c(),
                                     showdb = c())
        lcl$prev_struct <<- my_selection$struct
        if(nrow(rev_matches) == 0){
          shown_matches$reverse <- data.table::data.table(result = "No matches found.")
          return(NULL)
        }else{
          pieinfo$db <- reshape::melt(table(rev_matches$source)) # TODO: remove
          pieinfo$add <- reshape::melt(table(rev_matches$adduct))
          pieinfo$iso <- reshape::melt(table(rev_matches$isocat))
        }
      }
      mzMode =if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
      rev_matches = get_prematches(who = my_selection$struct,
                                                what = "con.structure",
                                                patdb = lcl$paths$patdb,
                                                showadd = result_filters$add[[mzMode]],
                                                showiso = result_filters$iso,
                                                showdb = result_filters$db,
                                                showIsolabels = input$show_iso_labels)
      if(nrow(rev_matches)>0){
        isocols = c("n2H", "n13C", "n15N")
        if("n2H" %in% colnames(rev_matches)){
          extracols=isocols
        }else{
          extracols=c()
        }
        shown_matches$reverse <- unique(rev_matches[,c("query_mz", "adduct",
                                                       "%iso", "dppm", 
                                                       extracols)])
      }else{
        shown_matches$reverse <- data.table::data.table(result = "No matches found.")
      }
    }
  } 
})
