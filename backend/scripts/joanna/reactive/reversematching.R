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
  lcl$tables$browse_table <<- unique(as.data.table(rbindlist(cpd_list)))
  # render table for UI
  output$browse_tab <-DT::renderDataTable({
    remove_cols = c("description", "structure", "formula", "charge")
    remove_idx <- which(colnames(lcl$tables$browse_table) %in% remove_cols)
    # don't show some columns but keep them in the original table, so they can be used
    # for showing molecule descriptions, structure
    DT::datatable(lcl$tables$browse_table,
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 15),
                                 pageLength = 5,
                                 columnDefs = list(list(visible=FALSE, 
                                                        targets=remove_idx)))
    )
  })

})

# triggers on reverse searching TODO: fix this, it's broken
observeEvent(input$revsearch_mz, {
  curr_row <- input$browse_tab_rows_selected
  search_cmd <- lcl$tables$browse_table[curr_row,c('structure')][[1]]

  if(!mSet$metshiParams$prematched){
    print("Please perform pre-matching first to enable this feature!")
    return(NULL)
  }else{
    lcl$tables$hits_table <<- unique(get_prematches(who = search_cmd,
                              what = "map.structure", #map.mz as alternative
                              patdb = lcl$paths$patdb)[,c("query_mz", "adduct", "%iso")])
    shown_matches$reverse <- if(nrow(lcl$tables$hits_table) > 0){
      lcl$tables$hits_table
    }else{
      data.table('name' = "Didn't find anything ( •́ .̫ •̀ )")
    }
  }
})



