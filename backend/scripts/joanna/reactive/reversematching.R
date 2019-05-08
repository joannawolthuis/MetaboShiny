# triggers on clicking the 'browse database' function
observeEvent(input$browse_db,{
  # get all compounds in the selected databases
  cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
    browse_db(match.table)
  })
  # join the individual result tables together
  global$tables$browse_table <<- unique(as.data.table(rbindlist(cpd_list)))
  # render table for UI
  output$browse_tab <-DT::renderDataTable({
    remove_cols = c("description", "structure", "formula", "dppm", "charge")
    remove_idx <- which(colnames(local$tables$browse_table) %in% remove_cols)
    # don't show some columns but keep them in the original table, so they can be used
    # for showing molecule descriptions, structure
    DT::datatable(local$tables$browse_table,
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 15),
                                 pageLength = 5,
                                 columnDefs = list(list(visible=FALSE, targets=remove_idx)))
    )
  })

})

# triggers on reverse searching TODO: fix this, it's broken
observeEvent(input$revsearch_mz, {
  curr_row <- input$browse_tab_rows_selected
  # curr_row <- grep(local$tables$browse_table$description, pattern="Creatine riboside")
  # -------------------
  search_cmd <- local$tables$browse_table[curr_row,c('formula', 'charge')]
  # -------------------
  cpd_list <- lapply(local$vectors$db_search_list, FUN=function(match.table){
    get_mzs(search_cmd$formula, search_cmd$charge, match.table, patdb=local$paths$patdb)})
  # ------------------
  local$tables$hits_table <<- unique(as.data.table(rbindlist(cpd_list)))
  output$hits_tab <-DT::renderDataTable({
    DT::datatable(local$tables$hits_tab,
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
  })
})

observeEvent(input$browse_tab_rows_selected,{
  curr_row <- input$browse_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_def <- local$tables$browse_table[curr_row, description]
  output$browse_definition <- renderText(curr_def)
})

observeEvent(input$hits_tab_rows_selected,{
  curr_row <- input$hits_tab_rows_selected
  curr_mz <- local$tables$hits_table[curr_row, mzmed]
  if (is.null(curr_row)) return()
  # -----------------------------
  #TODO: this should be a function and not re-written
  # make miniplot for sidebar with current compound
  output$curr_plot <- plotly::renderPlotly({
    # --- ggplot ---
    ggplotSummary(mSet, curr_mz, shape.fac = input$shape_var, cols = local$aes$mycols, cf=global$functions$color.functions[[local$aes$spectrum]],
                  styles = input$ggplot_sum_style,
                  add_stats = input$ggplot_sum_stats, col.fac = input$col_var,txt.fac = input$txt_var,
                  plot.theme = global$functions$plot.themes[[local$aes$theme]],
                  font = local$aes$font)
  })
})
