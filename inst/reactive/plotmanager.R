# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
plotmanager <- shiny::reactiveValues()

# preload pca/plsda
shiny::observe({
  if(is.null(plotmanager$make)){
    NULL # if not reloading anything, nevermind
  }else{
    if(!is.null(mSet)){
      success = F
      try({
        for(do in plotmanager$make){
          #suppressWarnings({
          toWrap <- MetaboShiny::getPlots(do, mSet, input, gbl, lcl, venn_yes, my_selection)
          lcl <<- toWrap$lcl
          toWrap <- toWrap$plots
          
          mapply(function(myplot, plotName){
            
            # if(grepl("aov|tt|fc|pattern|asca|volc", plotName)){
            #   whichAnal <- stringr::str_match(plotName, "aov|tt|fc|pattern|asca|volc")[,1]
            #   if(is.null(mSet$analSet[[whichAnal]]$sig.mat)){
            #     data = data.frame(text = "No significant hits!")
            #     myplot = ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10)
            #   }
            # }
            
            isSquare <- grepl("pca|plsda|tsne|roc|heatmap|var|samp|network", plotName) & !grepl("scree|cv|perm|venn", plotName)
            
            # === WRAPPER ===
            empty <- if(grepl(plotName, pattern="var|samp")) "output_empty2_width" else "output_empty3_width"
            
            output[[paste0(plotName, "_wrap")]] <- shiny::renderUI({
              list(conditionalPanel(
                condition = 'input.ggplotly == true',
                if(plotName != "network"){
                  plotly::plotlyOutput(paste0(plotName, "_interactive"), height = "100%") 
                }else{
                  visNetwork::visNetworkOutput(paste0(plotName, "_interactive"))
                }),
                conditionalPanel(
                  condition = 'input.ggplotly == false',
                  list(fluidRow(align="right",
                                downloadButton(outputId = paste0("download_", plotName),
                                               label = icon("Click to download"))),
                       plotOutput(plotName, height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)
                  )
                ))
            })
            
            # === PLOTS ===
            
            observe({
              
              canBe3D <- grepl("pca|plsda|tsne", plotName) & !grepl("scree|perm|cv", plotName)
              if(canBe3D){
                whichAnal <- stringr::str_match(plotName, "pca|plsda|tsne")[,1]
                is3D <- !input[[paste0(whichAnal, "_2d3d")]]
              }else{
                is3D <- plotName %in% c("heatmap", "network", "network_heatmap")
              }
              
              if(!is.null(session$clientData[[empty]])){
                
                try({
                  output[[plotName]] <- shiny::renderPlot({
                    # if(plotName %in% c("heatmap", "network_heatmap", "network")){
                    #   data = data.frame(text = "Currently only available in 'plotly' mode!\nPlease switch in the sidebar.")
                    #   ggplot2::ggplot(data) + ggplot2::geom_text(ggplot2::aes(label = text), x = 0.5, y = 0.5, size = 10) +
                    #     ggplot2::theme(text = ggplot2::element_text(family = lcl$aes$font$family)) + ggplot2::theme_bw()
                    # }else{
                      suppressWarnings({myplot})
                    #}
                  })  
                }, silent = F)
                
                plotFn <- paste0(c(gsub(":|,:", "_", mSet$dataSet$cls.name), 
                                   plotName), collapse="_")
                
                # emptyax <- list(
                #   title = "",
                #   zeroline = FALSE,
                #   showline = FALSE,
                #   showticklabels = FALSE,
                #   showgrid = FALSE
                # )
                
                output[[paste0(plotName, "_interactive")]] <- 
                  
                  if(plotName == "network"){
                    visNetwork::renderVisNetwork(myplot)
                  }else{
                    plotly::renderPlotly({
                      
                      if(!is3D){
                        myplot <- plotly::ggplotly(myplot,
                                                   tooltip = "text", 
                                                   height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)  
                      }
                      if(plotName != "heatmap"){
                        myplot <- if(grepl("venn", plotName)) myplot %>% plotly::layout(xaxis = emptyax,
                                                                                        yaxis = emptyax,
                                                                                        showlegend=F) else myplot %>% plotly::layout(showlegend=F) 
                      }else{
                        myplot <- myplot %>% plotly::layout(height = session$clientData[[empty]]/1.4,
                                                            width = session$clientData[[empty]])
                      }
                      suppressWarnings({
                        myplot %>%
                          plotly::config(
                            toImageButtonOptions = list(
                              format = if(input$plotsvg) "svg" else "png",
                              filename = paste0(plotFn, "_interactive")
                            ))   
                      })
                    })
                  }
                
                output[[paste0("download_", plotName)]] <- downloadHandler(
                  filename = function() paste0(plotFn, if(input$plotsvg) ".svg" else ".png"),
                  content = function(file){
                    ggplot2::ggsave(file, plot = myplot)
                  }
                )
              }
            })
          }, toWrap, names(toWrap))
        }
        success = T
      })
      if(!success){
        MetaboShiny::metshiAlert("Data plotting failed!")
      }
    }
    plotmanager$make <- NULL # set reloading to 'off'
  }
})
