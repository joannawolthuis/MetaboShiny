# make miniplot for sidebar with current compound
output$curr_plot <- plotly::renderPlotly({
  if(my_selection$mz != ""){
    MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, 
                               cols = lcl$aes$mycols, cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                               styles = input$ggplot_sum_style,
                               add_stats = input$ggplot_sum_stats, 
                               color.fac = input$col_var,
                               text.fac = input$txt_var,
                               fill.fac = input$fill_var,
                               plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                               font = lcl$aes$font)
  }
})

observeEvent(input$reload_plots, {
  if(!is.null(mSet)){
    if(!is.null(input$statistics)){
      uimanager$refresh <- input$statistics
      if(input$statistics %in% names(mSet$analSet)){
        plotmanager$make <- input$statistics
      }
    }  
  }
})

# pie charts
lapply(c("add", "iso", "db"), function(which_pie){
  output[[paste0("match_pie_", which_pie)]] <- plotly::renderPlotly({
    
    pievec = pieinfo[[which_pie]]
    
    if(nrow(pievec) > 0){
      
      m <- list(
        l = 0,
        r = 0,
        b = 30,
        t = 30)
      
      pulls = rep(0, nrow(pievec))
      lines = rep(1, nrow(pievec))
      
      if(!is.null(pievec)){
        if(which_pie == "add"){
          mzMode = if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
          targets = result_filters$add[[mzMode]]
        }else{
          targets = result_filters[[which_pie]]
        }
        pulls[which(as.character(pievec$Var.1) %in% targets)] <- 0.15  
        lines[which(as.character(pievec$Var.1) %in% targets)] <- 4
      }
      
      myCols <- gbl$functions$color.functions[[lcl$aes$spectrum]](n = nrow(pievec))
      
      if(length(pievec)>0){
        p = plotly::plot_ly(pievec, labels = ~Var.1, 
                            values = ~value, size=~value*10, type = 'pie',
                            textposition = 'inside',
                            textinfo = 'label+percent',
                            insidetextfont = list(colors = ggdark::invert_color(myCols)),
                            hoverinfo = 'text',
                            pull = pulls,
                            text = ~paste0(Var.1, ": ", value, ' matches'),
                            marker = list(colors = myCols,
                                          line = list(color = "gray", 
                                                      width = lines)),
                            #The 'pull' attribute can also be used to create space between the sectors
                            showlegend = FALSE) %>%
          plotly::layout(autosize = T, margin = m,
                         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))  
      }
      return(p)          
    }
  })
})

# ===== UI SWITCHER ====

shiny::observeEvent(input$heattable,{
  if(!is.null(input$overview)){
    if(input$overview == "heatmap"){
      statsmanager$calculate <- "heatmap"
      plotmanager$make <- "heatmap"
      uimanager$refresh <- "heatmap"
    }  
  }
})

shiny::observeEvent(input$network_style, {
  if("network" %in% names(mSet$analSet)){
    plotmanager$make <- "network"
  }
})

observeEvent(c(input$heatsign, input$heatlimits),{
  if(!is.null(mSet)){
    if("heatmap" %in% names(mSet$analSet)){
      plotmanager$make <- "heatmap"  
    }  
  }
})

