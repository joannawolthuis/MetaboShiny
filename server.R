# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).
library(datasets)

# Define a server for the Shiny app
function(input, output) {
  
  # Fill in the spot we created for a plot
  output$mzPlot <- renderPlot({
    # Render a plot
    PlotCmpdSummary(input$mz)
  })
}