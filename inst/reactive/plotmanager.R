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
        
        emptyax <- list(
          title = "",
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = FALSE
        )
        
        for(do in plotmanager$make){

          toWrap <- getPlots(do, mSet, 
                             input, gbl, 
                             lcl, venn_yes, 
                             my_selection)
          
          lcl <<- toWrap$lcl
          toWrap <- toWrap$plots
          
          isHeatmap = grepl("heatmap_", names(toWrap))
          toWrap$heatmap_plot = toWrap[isHeatmap]
          toWrap[which(isHeatmap)] <- NULL
          
          mapply(function(myplot, plotName){
            
            isSquare <- grepl("pca|plsda|tsne|roc|heatmap|var|samp|network|umap|ica", plotName) & !grepl("scree|cv|perm|venn", plotName)
            
            # === WRAPPER ===
            empty <- if(grepl(plotName, pattern="var|samp")) "output_empty2_width" else "output_empty3_width"
            
            output[[paste0(plotName, "_wrap")]] <- shiny::renderUI({
              if(plotName != "network"){
                list(conditionalPanel(
                  condition = 'input.ggplotly == true',
                  plotly::plotlyOutput(paste0(plotName, "_interactive"), height = "100%")) ,
                  conditionalPanel(
                    condition = 'input.ggplotly == false',
                    list(fluidRow(align="right",
                                  downloadButton(outputId = paste0("download_", plotName),
                                                 label = icon("download"))),
                         plotOutput(plotName, height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)
                    )
                  )) 
              }else{
                visNetwork::visNetworkOutput(paste0(plotName, "_interactive"),
                                             height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)
              }
            })
            
            # === PLOTS ===
            
            observe({
              
              canBe3D <- grepl("pca|plsda|tsne|umap|ica", plotName) & !grepl("scree|perm|cv", plotName)
              if(canBe3D){
                whichAnal <- stringr::str_match(plotName, "pca|plsda|tsne|umap|ica")[,1]
                is3D <- !input[[paste0(whichAnal, "_2d3d")]]
              }else{
                is3D <- plotName %in% c("network", 
                                        "network_heatmap")
              }
              
              if(!is.null(session$clientData[[empty]])){
                if(!(plotName %in% c("network",
                                     "wordcloud"))){
                  try({
                    if(length(myplot) > 1){
                      if(plotName == "heatmap_plot") myplot$heatmap_static() else{
                        if(length(myplot$layers[[1]]$data) > 0){
                          myplot$data = myplot$layers[[1]]$data
                        }
                        if(input$plot_mzlabels & (
                          any(grepl("mz|m/z", names(myplot$data)))
                        )){
                          if(length(myplot$layers[[1]]$mapping) > 0){
                            myplot$mapping = myplot$layers[[1]]$mapping
                          }
                          myX = rlang::quo_get_expr(myplot$mapping[['x']])
                          myY = rlang::quo_get_expr(myplot$mapping[['y']])
                          myText = rlang::quo_get_expr(myplot$mapping[['text']])
                          myCol = rlang::quo_get_expr(myplot$mapping[['colour']])
                          flip = grepl("tt|fc|aov|var|samp|corr", plotName)
                          
                          if(length(myplot$data) == 0){
                            myplot$data = myplot$layers[[1]]$data  
                          }
                          
                          if("significant" %in% colnames(myplot$data)){
                            plotdata = myplot$data[significant == "YES"]
                          }else{
                            plotdata = myplot$data
                          }
                          
                          myplot = myplot + ggrepel::geom_label_repel(data = plotdata,
                                                                      aes_string(y = myY,
                                                                                 x = myX,
                                                                                 label = myText),
                                                                      color="black",
                                                                      size = 5)
                        }
                      }
                      
                      output[[plotName]] <- shiny::renderPlot({
                        suppressWarnings(myplot)
                      })  
                    }
                  }, silent = F)
                  plotFn <- paste0(c(gsub(":|,:", "_", mSet$settings$cls.name), 
                                     plotName), collapse="_") 
                  if(grepl(x=plotFn, "ml")){
                    plotFn <- paste(plotFn, mSet$analSet$ml$last$method, mSet$analSet$ml$last$name, sep = "_")
                  }
                }
                
                output[[paste0(plotName, "_interactive")]] <- 
                  
                  if(plotName == "network"){
                    visNetwork::renderVisNetwork(myplot)
                  }else if(plotName == "wordcloud"){
                    wordcloud2::renderWordcloud2(myplot)
                  }else{
                    plotly::renderPlotly({
                      if(!is3D & plotName != "heatmap_plot"){
                        myplot <- plotly::ggplotly(suppressWarnings(myplot),
                                                   tooltip = "text", 
                                                   height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)  
                      }
                      if(plotName != "heatmap_plot"){
                        myplot <- if(grepl("venn", plotName)) plotly::ggplotly(suppressWarnings(myplot)) %>% plotly::layout(xaxis = emptyax,
                                                                                           yaxis = emptyax,
                                                                                           showlegend=input$legend) else myplot %>% plotly::layout(showlegend=input$legend) 
                      }else{
                        myplot <- if(plotName == "heatmap_plot") suppressWarnings(myplot$heatmap_interactive) else plotly::ggplotly(suppressWarnings(myplot)) %>% plotly::layout(height = session$clientData[[empty]]/1.4,
                                                               width = session$clientData[[empty]])
                      }
                      
                      if(canBe3D){
                        try({
                          if(length(myplot$x$data) > 0){
                            for(i in 1:length(myplot$x$data)){
                              if(myplot$x$data[[i]]$hoveron == "fills"){
                                myplot$x$data[[i]]$hoverinfo <- "skip"
                              }
                            }  
                          }
                        }, silent = T)
                      }
                      
                      suppressWarnings({
                        suppressWarnings(myplot %>%
                          plotly::config(
                            toImageButtonOptions = list(
                              format = if(input$plotsvg) "svg" else "png",
                              filename = paste0(plotFn, "_interactive")
                            ))
                          )
                      })
                    })
                  }
                
                output[[paste0("download_", plotName)]] <- downloadHandler(
                  filename = function() paste0(plotFn, if(input$plotsvg) ".svg" else ".png"),
                  content = function(file){
                    if(plotName == "heatmap_plot"){
                      saveFun = if(input$plotsvg) svg else png
                      saveFun(file=file)
                      suppressWarnings(myplot$heatmap_static())
                      dev.off()  
                    }else{
                      ggplot2::ggsave(file, plot = suppressWarnings(myplot))
                    }
                  }
                )
              }
            })
          }, toWrap, names(toWrap))
        }
        success = T
      })
      if(!success){
        metshiAlert("Data plotting failed!")
      }
    }
    plotmanager$make <- NULL # set reloading to 'off'
  }
})
