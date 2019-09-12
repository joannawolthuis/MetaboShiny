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
  my_selection$structure <- browse_content$table[curr_row,c('structure')][[1]]
})



