# triggers when a new color spectrum is chosen
observeEvent(input$color_ramp,{
  # render preview plot
  output$ramp_plot <- plotly::renderPlotly({
    
    # change the current spectrum in global
    global$functions$color.functions[[getOptions("user_options.txt")$gspec]] <<- global$functions$color.functions[[input$color_ramp]]
    # change the current spectrum in user options file
    setOption("user_options.txt", "gspec", input$color_ramp)
    
    # create plot background (no grid, no lines, just color ;) )
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE,
      titlefont = list(size = 20)
    )
    
    # re-render preview plot with the new options (general heatmap using R standard volcano dataset)
    plotly::plot_ly(z = volcano, 
                    colors = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](100), 
                    type = "heatmap",
                    showscale=FALSE)  %>%
      layout(xaxis = ax, yaxis = ax)
  })
})

# triggers when new plot theme is picked
observeEvent(input$ggplot_theme,{
  
  # change default plot theme in user settings
  setOption("user_options.txt", "gtheme", input$ggplot_theme)
  
  # change preview plot (uses mtcars default R dataset)
  output$ggplot_theme_example <- renderPlot({
    p <- ggplot(mtcars) + geom_boxplot(aes(x = wt, y = mpg,
                                           colour = factor(gear)))
    p + global$functions$plot.themes[[getOptions("user_options.txt")$gtheme]]()
  })
})

# triggers when changes in interface aesthetics are applied
observeEvent(input$change_css, {
  
  # set default user color options
  setOption("user_options.txt", "col1", input$bar.col.1)
  setOption("user_options.txt", "col2", input$bar.col.2)
  setOption("user_options.txt", "col3", input$bar.col.3)
  setOption("user_options.txt", "col4", input$bar.col.4)
  
  # set default user font options
  setOption("user_options.txt", "font1", input$font.1)
  setOption("user_options.txt", "font2", input$font.2)
  setOption("user_options.txt", "font3", input$font.3)
  setOption("user_options.txt", "font4", input$font.4)
  
  # set default user font size options
  setOption("user_options.txt", "size1", input$size.1)
  setOption("user_options.txt", "size2", input$size.2)
  setOption("user_options.txt", "size3", input$size.3)
  setOption("user_options.txt", "size4", input$size.4)
})