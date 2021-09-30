# make miniplot for sidebar with current compound
output$curr_plot <- plotly::renderPlotly({
  if(my_selection$mz != ""){
    plotmananager$make <- "summary"
  }
})

observeEvent(input$reload_plots, {
  if(!is.null(mSet)){
    if(!is.null(input$statistics)){
      uimanager$refresh <- input$statistics
      if(input$statistics %in% names(mSet$analSet)){
        plot.me = input$statistics
        if(my_selection$mz != ""){
          plot.me = c(plot.me, "summary")
        }
        plotmanager$make <- plot.me
      }
    }
  }
})

# pie charts
lapply(c("add", "iso", "db"), function(which_pie){
  output[[paste0("match_pie_", which_pie)]] <- plotly::renderPlotly({
    
    pievec = pieinfo[[which_pie]]
    if(!is.null(pievec)){
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
      }else{
        data = data.frame(text = "Please run \n a search!")
        p = ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
          ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
        plotly::ggplotly(p)
    }
    }else{
      data = data.frame(text = "Please run \n a search!")
      p = ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
        ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
      plotly::ggplotly(p)
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

