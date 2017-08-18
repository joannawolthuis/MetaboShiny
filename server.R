function(input, output, session) {

# ======================== DB CHECK ============================



# ======================== ASCA ============================

output$asca_tab <- DT::renderDataTable({
  datatable(asca.table, 
            selection = 'single',
            colnames = c("Mass/charge", "Leverage", "SPE"),
            autoHideNavigation = T,
            options = list(lengthMenu = c(10, 30, 50), pageLength = 10))
})

# check for selected mz row
observeEvent(input$asca_tab_rows_selected,{
  curr_row = input$asca_tab_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- asca.table[curr_row,'X']
  output$asca_plot <- renderPlot(PlotCmpdSummary(curr_mz))
})

# --- find matches ---

observeEvent(input$do,{
  match_list <- lapply(input$checkGroup, FUN=function(match.table){
    get_matches(curr_mz, match.table)})
  result_table <- as.data.table(rbindlist(match_list))
  output$match_tab <- DT::renderDataTable({
    datatable(unique(result_table[Compound != "",]),
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })
})

}