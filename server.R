library(ggplot2)
library(DT)
# choose columns to display

function(input, output, session) {
  
output$x1 <- DT::renderDataTable({
  datatable(asca.table, selection = 'single')
})

observeEvent(input$x1_rows_selected, {
  curr.row = input$x1_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr.row)) return()
  curr.mz <- asca.table[curr.row,'X']
  output$plot1 <- renderPlot(PlotCmpdSummary(curr.mz))
})

}