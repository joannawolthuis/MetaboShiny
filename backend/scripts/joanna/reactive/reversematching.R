# triggers on clicking the 'browse database' function
observeEvent(input$browse_db,{
  # get all compounds in the selected databases
  cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
    browse_db(match.table)
  })
  # join the individual result tables together
  browse_table <<- unique(as.data.table(rbindlist(cpd_list)))
  # render table for UI
  output$browse_tab <-DT::renderDataTable({
    DT::datatable(browse_table[,-c("Description", "Charge")],
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
  })
})

# triggers on reverse searching TODO: fix this, it's broken
# observeEvent(input$revsearch_cpd, {
#   input$browse_tab_rows_selected
#   # -------------------
#   search_cmd <- browse_table[curr_row,c('Formula', 'Charge')]
#   # -------------------
#   cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
#     get_mzs(search_cmd$Formula, search_cmd$Charge, match.table)})
#   # ------------------
#   hits_table <<- unique(as.data.table(rbindlist(cpd_list)))
#   output$hits_tab <-DT::renderDataTable({
#     DT::datatable(hits_table,
#                   selection = 'single',
#                   autoHideNavigation = T,
#                   options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
#   })
# })

lapply(c("browse", "hits"), function(prefix){
  observeEvent(input$browse_tab_rows_selected,{
    curr_row <<- input[[paste0(prefix, "_tab_rows_selected")]]
    curr_cpd <- browse_table[curr_row, Formula]
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_def <<- browse_table[curr_row, switch(prefix, browse='Description', hits="mzmed.pgrp")]
    if(prefix == "browse"){
      output$browse_definition <- renderText(curr_def$Description)
    }
    #TODO: this should be a function and not re-written
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, shape.fac = input$second_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])})
  })
})