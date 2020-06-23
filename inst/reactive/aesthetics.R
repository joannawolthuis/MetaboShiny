output$ramp_plot <- plotly::renderPlotly({
  # change the current spectrum in global
  lcl$aes$spectrum <<- input$color_ramp
  # change the current spectrum in user options file
  MetaboShiny::setOption(lcl$paths$opt.loc, 
                         key="gspec", 
                         value=input$color_ramp)
  
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
                  colors = gbl$functions$color.functions[[lcl$aes$spectrum]](100),
                  type = "heatmap",
                  showscale=FALSE)  %>%
    plotly::layout(xaxis = ax, yaxis = ax)
})


# change preview plot (uses mtcars default R dataset)
  output$ggplot_theme_example <- plotly::renderPlotly({
    MetaboShiny::setOption(lcl$paths$opt.loc, key="gtheme", value=input$ggplot_theme)
    lcl$aes$theme <<- input$ggplot_theme
    
    p <- ggplot2::ggplot(mtcars) + 
      ggplot2::geom_boxplot(ggplot2::aes(x = wt, y = mpg,
                       colour = factor(gear))) +
      gbl$functions$plot.themes[[lcl$aes$theme]]()
    
    plotly::ggplotly(p)
  })


# triggers when changes in interface aesthetics are applied
shiny::observeEvent(input$change_css, {

  # set default user color options
  MetaboShiny::setOption(lcl$paths$opt.loc, key="col1", value=input$bar.col.1)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="col2", value=input$bar.col.2)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="col3", value=input$bar.col.3)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="col4", value=input$bar.col.4)

  # set default user font options
  MetaboShiny::setOption(lcl$paths$opt.loc, key="font1", value=input$font.1)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="font2", value=input$font.2)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="font3", value=input$font.3)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="font4", value=input$font.4)

  # set default user font size options
  MetaboShiny::setOption(lcl$paths$opt.loc, key="size1", value=input$size.1)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="size2", value=input$size.2)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="size3", value=input$size.3)
  MetaboShiny::setOption(lcl$paths$opt.loc, key="size4", value=input$size.4)

})
