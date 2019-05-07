# triggers when a new color spectrum is chosen
observeEvent(input$color_ramp,{
  # render preview plot
  output$ramp_plot <- plotly::renderPlotly({

    # change the current spectrum in global
    global$functions$color.functions[[getOptions(opt.loc)$gspec]] <<- global$functions$color.functions[[input$color_ramp]]
    # change the current spectrum in user options file
    setOption(opt.loc, key="gspec", value=input$color_ramp)

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
                    colors = global$functions$color.functions[[getOptions(opt.loc)$gspec]](100),
                    type = "heatmap",
                    showscale=FALSE)  %>%
      layout(xaxis = ax, yaxis = ax)
  })
})

# triggers when new plot theme is picked
observeEvent(input$ggplot_theme,{

  print(opt.loc)
  
  # change default plot theme in user settings
  setOption(opt.loc, key="gtheme", value=input$ggplot_theme)

  # change preview plot (uses mtcars default R dataset)
  output$ggplot_theme_example <- renderPlot({
    p <- ggplot(mtcars) + geom_boxplot(aes(x = wt, y = mpg,
                                           colour = factor(gear)))
    p + global$functions$plot.themes[[getOptions(opt.loc)$gtheme]]()
  })
})

# triggers when changes in interface aesthetics are applied
observeEvent(input$change_css, {

  # set default user color options
  setOption(opt.loc, key="col1", value=input$bar.col.1)
  setOption(opt.loc, key="col2", value=input$bar.col.2)
  setOption(opt.loc, key="col3", value=input$bar.col.3)
  setOption(opt.loc, key="col4", value=input$bar.col.4)

  # set default user font options
  setOption(opt.loc, key="font1", value=input$font.1)
  setOption(opt.loc, key="font2", value=input$font.2)
  setOption(opt.loc, key="font3", value=input$font.3)
  setOption(opt.loc, key="font4", value=input$font.4)

  # set default user font size options
  setOption(opt.loc, key="size1", value=input$size.1)
  setOption(opt.loc, key="size2", value=input$size.2)
  setOption(opt.loc, key="size3", value=input$size.3)
  setOption(opt.loc, key="size4", value=input$size.4)

})
