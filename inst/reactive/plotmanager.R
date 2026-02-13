# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
plotmanager <- shiny::reactiveValues()
plotDims <- shiny::reactiveValues()

dbg_plots <- function(...){
  if (!isTRUE(getOption("metaboshiny.debug_plots", FALSE)) &&
      Sys.getenv("METABOSHINY_DEBUG_PLOTS") != "1") {
    return(invisible(NULL))
  }
  msg <- paste0("[metaboshiny][plots] ", paste(..., collapse = ""))
  message(msg)
  invisible(NULL)
}

get_plot_dims_px <- function(session, plotName){
  width_px <- session$clientData[[paste0("output_", plotName, "_width")]]
  height_px <- session$clientData[[paste0("output_", plotName, "_height")]]
  if (is.null(width_px) || is.null(height_px) || width_px <= 0 || height_px <= 0) {
    return(list(width_px = 900, height_px = 650))
  }
  list(width_px = width_px, height_px = height_px)
}

is_diag_plots <- function(){
  isTRUE(getOption("metaboshiny.diag_plots", FALSE)) ||
    identical(Sys.getenv("METABOSHINY_DIAG_PLOTS"), "1") ||
    isTRUE(get0("diag_plots", ifnotfound = FALSE, inherits = TRUE))
}

screen_dpi <- 96
download_dpi <- 300

# preload pca/plsda
shiny::observe({
  if(is.null(plotmanager$make)){
    NULL # if not reloading anything, nevermin
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
              plot_h <- "500px"
              dbg_plots("renderUI wrap plotName=", plotName,
                        " ggplotly=", isTRUE(input$ggplotly),
                        " plot_h=", plot_h)
              if(plotName != "network"){
                if (isTRUE(input$ggplotly)) {
                  list(
                    fluidRow(
                      align = "right",
                      downloadButton(outputId = paste0("download_", plotName), label = "")
                    ),
                    div(
                      style = paste0("height:", plot_h, "; width: 100%;"),
                      plotly::plotlyOutput(paste0(plotName, "_interactive"), height = "100%")
                    )
                  )
                } else {
                  list(
                    fluidRow(
                      align = "right",
                      downloadButton(outputId = paste0("download_", plotName), label = "")
                    ),
                    div(
                      style = paste0("height:", plot_h, "; width: 100%;"),
                      plotOutput(plotName, height = "100%")
                    )
                  )
                }
              }else{
                visNetwork::visNetworkOutput(paste0(plotName, "_interactive"))
                                             #,height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)
              }
            })
            
            # === PLOTS ===
            
            observe({
              
              canBe3D <- grepl("pca|plsda|tsne|umap|ica", plotName) & !grepl("scree|perm|cv", plotName)
              if(canBe3D){
                whichAnal <- stringr::str_match(plotName, "pca|plsda|tsne|umap|ica")[,1]
                is3D <- !input[[paste0(whichAnal, "_2d3d")]]
              }else{
                is3D <- plotName %in% c("network")
              }
              
              if(TRUE){
                if(!(plotName %in% c("network",
                                     "wordcloud",
                                     "ml_roc"))){
                  try({
                    dbg_plots("observe plotName=", plotName,
                              " class=", paste(class(myplot), collapse = "|"),
                              " length=", length(myplot))
                    if(is.null(myplot)) {
                      dbg_plots("observe skip plotName=", plotName, " reason=myplot NULL")
                    } else {
                      if(plotName == "heatmap_plot") {
                        if (is.list(myplot) && is.function(myplot$heatmap_static)) {
                          myplot$heatmap_static()
                        } else {
                          dbg_plots("observe skip heatmap static: invalid heatmap object")
                        }
                      } else {
                        if(length(myplot$layers[[1]]$data) > 0){
                          myplot$data = myplot$layers[[1]]$data
                        }
                        if(input$plot_mzlabels & !isTRUE(input$ggplotly) & (
                          any(grepl("mz|m/z", names(myplot$data)))
                        )){
                          if(length(myplot$layers[[1]]$mapping) > 0){
                            myplot$mapping = myplot$layers[[1]]$mapping
                          }
                          myX = rlang::quo_get_expr(myplot$mapping[['x']])
                          myY = rlang::quo_get_expr(myplot$mapping[['y']])
                          myText = rlang::quo_get_expr(myplot$mapping[['text']])
                          myCol = rlang::quo_get_expr(myplot$mapping[['colour']])
                          flip = grepl("tt|fc|aov|var|samp|corr|cliffd", plotName)
                          
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
                      
                      # output[[plotName]] <- shiny::renderPlot({
                      #   suppressWarnings(myplot)
                      # })  
                      observe({
                        # Dynamically update dimensions for each plot based on plotName
                        dims <- plotDims[[plotName]]
                        if (is.null(dims)) dims <- list(width = NULL, height = NULL)
                        dims$width <- session$clientData[[paste0("output_", plotName, "_width")]]
                        dims$height <- session$clientData[[paste0("output_", plotName, "_height")]]
                        plotDims[[plotName]] <- dims
                      })
                      
                      output[[plotName]] <- shiny::renderPlot({
                        dbg_plots("renderPlot start plotName=", plotName)

                        showtext::showtext_opts(dpi = screen_dpi)

                        # -- fix ticks? --
                        myplot <- myplot + theme(
                          axis.ticks = element_line(colour = "black", linewidth = .5),
                          axis.ticks.length = unit(0.075, "cm")
                        )

                        tryCatch({
                          if(plotName == "heatmap_plot"){
                            if (is.list(myplot) && is.function(myplot$heatmap_static)) {
                              myplot$heatmap_static()
                            } else {
                              stop("Invalid heatmap plot object: missing heatmap_static()")
                            }
                          }else{
                            suppressWarnings(print(myplot))
                          }
                        }, error = function(e) {
                          msg <- paste0("Plot rendering failed (", plotName, "): ", conditionMessage(e))
                          dbg_plots(msg)
                          try(metshiAlert(msg), silent = TRUE)
                          try(shiny::showNotification(msg, type = "error"), silent = TRUE)
                          stop(e)
                        })
                      }, res = screen_dpi)
                    }
                  }, silent = F)
                }
                
                plotFn <- paste0(c(gsub(":|,:", "_", mSet$settings$cls.name), 
                                   plotName), collapse="_") 
                if(grepl(x=plotFn, "ml")){
                  plotFn <- paste(plotFn, 
                                  mSet$analSet$ml$last$method, 
                                  mSet$analSet$ml$last$name, sep = "_")
                }
                
                output[[paste0("download_", plotName)]] <- downloadHandler(
                  filename = function() paste0(plotFn, if(input$plotsvg) ".svg" else ".png"),
                  content = function(file){
                    if(plotName == "heatmap_plot"){
                      saveFun(file=file)
                      if (is.list(myplot) && is.function(myplot$heatmap_static)) {
                        suppressWarnings(myplot$heatmap_static())
                        dev.off()  
                      } else {
                        stop("Invalid heatmap plot object for download")
                      }
                    }else{
                      dims <- get_plot_dims_px(session, plotName)
                      dbg_plots("download plotName=", plotName,
                                " width_px=", dims$width_px,
                                " height_px=", dims$height_px,
                                " plotsvg=", input$plotsvg)

                      if(input$plotsvg){
                        width_in <- dims$width_px / screen_dpi
                        height_in <- dims$height_px / screen_dpi
                        ggsave(
                          filename = file,
                          plot = myplot,
                          device = "svg",
                          width = width_in,
                          height = height_in,
                          units = "in"
                        )
                      }else{
                        width_in <- dims$width_px / screen_dpi
                        height_in <- dims$height_px / screen_dpi
                        showtext::showtext_opts(dpi = download_dpi)
                        ggsave(
                          filename = file,
                          plot = myplot,
                          device = "png",
                          width = width_in,
                          height = height_in,
                          units = "in",
                          dpi = download_dpi
                        )
                        showtext::showtext_opts(dpi = screen_dpi)
                      }
                    }
                  }
                )
                
                output[[paste0(plotName, "_interactive")]] <- 
                  
                  if(plotName == "network"){
                    print("!network change")
                    visNetwork::renderVisNetwork({
                      myplot %>% visNetwork::visExport(type="png", name=plotFn, label="download png")
                      })
                  }else if(plotName == "wordcloud"){
                    wordcloud2::renderWordcloud2(myplot)
                  }else{
                    plotly::renderPlotly({
                      dbg_plots("renderPlotly start plotName=", plotName, " is3D=", is3D, " canBe3D=", canBe3D)
                      if(!is3D & plotName != "heatmap_plot"){
                        myplot <- withCallingHandlers(
                          plotly::ggplotly(
                            myplot,
                            tooltip = "text"
                            #height = session$clientData[[empty]]/if(isSquare) 1.4 else 2
                          ),
                          warning = function(w) {
                            msg <- conditionMessage(w)
                            if (startsWith(msg, "Ignoring unknown aesthetics:")) {
                              aes_txt <- trimws(sub("^Ignoring unknown aesthetics:\\s*", "", msg))
                              aes_txt <- gsub("\\s+and\\s+", ",", aes_txt)
                              aes_txt <- gsub("\\s+", "", aes_txt)
                              aes_vec <- unlist(strsplit(aes_txt, ",", fixed = TRUE))
                              aes_vec <- aes_vec[nzchar(aes_vec)]
                              
                              if (length(aes_vec) > 0 && all(aes_vec %in% c("text", "key"))) {
                                invokeRestart("muffleWarning")
                              }
                            }
                          }
                        )
                        dbg_plots("renderPlotly ggplotly ok plotName=", plotName)
                      }
                      if(plotName != "heatmap_plot"){
                        myplot <- if(grepl("venn", plotName)) plotly::ggplotly(suppressWarnings(myplot)) %>% plotly::layout(autosize = TRUE,
                                                                                           xaxis = emptyax,
                                                                                           yaxis = emptyax,
                                                                                           showlegend=input$legend) else myplot %>% plotly::layout(showlegend=input$legend) 
                      }else{
                        if (is.list(myplot) && "heatmap_interactive" %in% names(myplot)) {
                          myplot <- suppressWarnings(myplot$heatmap_interactive)
                          if (is.function(myplot)) {
                            myplot <- suppressWarnings(myplot())
                          }
                        } else {
                          stop("Invalid heatmap plot object: missing heatmap_interactive")
                        }
                        # %>% plotly::layout(height = session$clientData[[empty]]/1.4,
                        #width = session$clientData[[empty]])
                      }
                      if(plotName != "heatmap_plot" && !grepl("venn", plotName)){
                        myplot <- plotly::layout(myplot, autosize = TRUE)
                      }
                      dbg_plots("renderPlotly layout ok plotName=", plotName)
                      
                      # Newer plotly versions require explicit event registration to
                      # receive `event_data("plotly_click")` without warnings.
                      myplot <- plotly::event_register(myplot, "plotly_click")
                      dbg_plots("renderPlotly event_register ok plotName=", plotName)

                      # Some browsers/shiny layouts can leave plotly containers hidden;
                      # force visibility at render time and log dimensions in console.
                      myplot <- htmlwidgets::onRender(
                        myplot,
                        "function(el, x) {
                           try {
                             console.log('[metshi][plotly] onRender', el.id, el.style.visibility, el.offsetWidth, el.offsetHeight);
                             el.style.visibility = 'visible';
                             if (el.parentElement) el.parentElement.style.visibility = 'visible';
                             var inner = el.querySelectorAll('.plot-container,.svg-container,.main-svg');
                             for (var i = 0; i < inner.length; i++) {
                               inner[i].style.visibility = 'visible';
                             }
                             setTimeout(function(){
                               try { if (window.Plotly) window.Plotly.Plots.resize(el); } catch(e) {}
                             }, 60);
                           } catch (e) {}
                         }"
                      )

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
                      
                      myplot <- suppressWarnings({
                        suppressWarnings(myplot %>%
                          plotly::config(
                            toImageButtonOptions = list(
                              format = if(input$plotsvg) "svg" else "png",
                              filename = paste0(plotFn, "_interactive")
                            ))
                          )
                      })
                      dbg_plots("renderPlotly done plotName=", plotName)
                      myplot
                    })
                  }
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
