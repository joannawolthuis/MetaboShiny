# triggers on clicking the 'browse database' function
observeEvent(input$browse_db,{
  # get all compounds in the selected databases
  cpd_list <- lapply(lcl$vectors$db_search_list, FUN=function(match.table){
    res = MetaDBparse::showAllBase(outfolder = lcl$paths$db_dir,
                                   base.dbname = gsub(basename(match.table), 
                                                      pattern = "\\.db",
                                                      replacement=""))
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
    if(!mSet$metshiParams$prematched){
      MetaboShiny::metshiAlert("Please perform pre-matching first to enable this feature!")
      return(NULL)
    }else{
      
      if(lcl$prev_struct != my_selection$struct){
        rev_matches = MetaboShiny::get_prematches(who = my_selection$struct,
                                                  what = "con.structure",
                                                  patdb = lcl$paths$patdb,
                                                  showadd = c(),
                                                  showiso = c(),
                                                  showdb = c())  
        lcl$prev_struct <<- my_selection$struct
        if(nrow(rev_matches) == 0){
          shown_matches$reverse <- data.table::data.table()
          return(NULL)
        }else{
          pieinfo$db <- reshape::melt(table(rev_matches$source))
          pieinfo$add <- reshape::melt(table(rev_matches$adduct))
          pieinfo$iso <- reshape::melt(table(rev_matches$isocat))
        }
      }
      mzMode =if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
      rev_matches = MetaboShiny::get_prematches(who = my_selection$struct,
                                                what = "con.structure",
                                                patdb = lcl$paths$patdb,
                                                showadd = result_filters$add[[mzMode]],
                                                showiso = result_filters$iso,
                                                showdb = result_filters$db)
      if(nrow(rev_matches)>0){
        shown_matches$reverse <- unique(rev_matches[,c("query_mz", "adduct", "%iso", "dppm")])
      }else{
        shown_matches$reverse <- data.table::data.table()
      }
    }
  } 
})
