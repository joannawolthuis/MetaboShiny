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
          toWrap = switch(do,
                          general = {
                            # make sidebar
                            # make pca, plsda, ml(make plotmanager do that)
                            # update select input bars with current variable and covariables defined in excel
                            if(is.null(mSet)){
                              interface$mode <- NULL
                            }else{
                              varNormPlots <- MetaboShiny::ggplotNormSummary(mSet = mSet,
                                                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              
                              sampNormPlots <- MetaboShiny::ggplotSampleNormSummary(mSet,
                                                                                    cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              
                              list(var1=varNormPlots$tl, var2=varNormPlots$bl,
                                   var3=varNormPlots$tr, var4=varNormPlots$br,
                                   samp1=sampNormPlots$tl, samp2=sampNormPlots$bl, 
                                   samp3=sampNormPlots$tr, samp4=sampNormPlots$br)
                              
                            }},
                          venn = {
                            # get user input for how many top values to use for venn
                            top = input$venn_tophits
                            if(nrow(venn_yes$now) > 4 | nrow(venn_yes$now) <= 1){
                              MetaboShiny::metshiAlert("Can only take more than 2 and less than five analyses!")
                              list()
                            }else{
                              p <- MetaboShiny::ggPlotVenn(mSet = mSet,
                                                           venn_yes = as.list(venn_yes),
                                                           top = input$venn_tophits,
                                                           cols = lcl$aes$mycols,
                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              lcl$vectors$venn_lists <<- p$info
                              shiny::updateSelectizeInput(session, "intersect_venn", choices = names(lcl$vectors$venn_lists))
                              list(venn_plot = p$plot)
                            }
                          },
                          enrich = {
                            NULL
                          },
                          summary = {
                            p = MetaboShiny::ggplotSummary(mSet, my_selection$mz, 
                                                           shape.fac = input$shape_var, 
                                                           cols = lcl$aes$mycols, 
                                                           cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                           styles = input$ggplot_sum_style,
                                                           add_stats = input$ggplot_sum_stats, 
                                                           color.fac = input$col_var, 
                                                           text.fac = input$txt_var)
                            
                            list(summary_plot = p)
                          },
                          pattern = {
                            p = MetaboShiny::ggPlotPattern(mSet,
                                                           n = input$pattern_topn,
                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                            
                            list(pattern_plot = p)
                          },
                          aov = { # render manhattan-like plot for UI
                            p = MetaboShiny::ggPlotAOV(mSet,
                                                       cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 20)
                            
                            list(aov_plot = p)
                          },
                          volc = {
                            # render volcano plot with user defined colours
                            p = MetaboShiny::ggPlotVolc(mSet,
                                                        cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                        20)
                            
                            list(volc_plot = p)
                          },
                          tsne = {
                            mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                              "ipca" # interactive PCA (old name, i like tpca more :P )
                            }else{
                              "normal" # normal pca
                            }
                            
                            if("tsne" %in% names(mSet$analSet)){
                              if(input$tsne_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                                # render 2d plot
                                p = MetaboShiny::plotPCA.2d(mSet, 
                                                            cols = lcl$aes$mycols,
                                                            pcx = 1,
                                                            pcy = 2, 
                                                            type = "tsne",
                                                            mode = mode,
                                                            shape.fac = input$shape_var,
                                                            col.fac = input$col_var
                                )
                                
                              }else{
                                # render 3d plot
                                p = MetaboShiny::plotPCA.3d(mSet, 
                                                            lcl$aes$mycols,
                                                            pcx = 1,
                                                            pcy = 2,
                                                            pcz = 3,
                                                            type = "tsne",
                                                            mode = mode,
                                                            shape.fac = input$shape_var,
                                                            font = lcl$aes$font,
                                                            col.fac = input$col_var,
                                                            cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }
                            }
                            list(tsne_plot = p)
                          },
                          pca = {
                            if("pca" %in% names(mSet$analSet)){
                              scree = MetaboShiny::ggPlotScree(mSet,
                                                               cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              # chekc which mode we're in
                              mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                                "ipca" # interactive PCA (old name, i like tpca more :P )
                              }else{
                                "normal" # normal pca
                              }
                              
                              if(input$pca_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                                # render 2d plot
                                pca = MetaboShiny::plotPCA.2d(mSet, 
                                                              cols = lcl$aes$mycols,
                                                              pcx = input$pca_x,
                                                              pcy = input$pca_y, 
                                                              mode = mode,
                                                              type = "pca",
                                                              col.fac = input$col_var,
                                                              shape.fac = input$shape_var)
                                loadings = MetaboShiny::plotPCAloadings.2d(mSet,pcx = input$pca_x,
                                                                           pcy = input$pca_y, 
                                                                           type = "pca",
                                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }else{
                                # render 3d plot
                                pca = MetaboShiny::plotPCA.3d(mSet, 
                                                              pcx = input$pca_x,
                                                              pcy = input$pca_y,
                                                              pcz = input$pca_z, 
                                                              type = "pca",
                                                              col.fac = input$col_var,
                                                              mode = mode,
                                                              cols = lcl$aes$mycols,
                                                              shape.fac = input$shape_var,
                                                              font = lcl$aes$font,
                                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                                loadings = MetaboShiny::plotPCAloadings.3d(mSet,
                                                                           pcx = input$pca_x,
                                                                           pcy = input$pca_y,
                                                                           pcz = input$pca_z, 
                                                                           type = "pca",
                                                                           font = lcl$aes$font,
                                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }
                              list(plot_pca = pca, plot_pca_loadings = loadings, pca_scree = scree)
                            }else{
                              NULL
                            }
                          },
                          plsda = {
                            
                            if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
                              
                              mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                                "ipca" # interactive PCA (old name, i like tpca more :P )
                              }else{
                                "normal" # normal pca
                              }
                              
                              # render cross validation plot
                              cv <- MetaboShiny::ggPlotClass(mSet, cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              # render permutation plot
                              perm <- MetaboShiny::ggPlotPerm(mSet,cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              # see PCA - render 2d or 3d plots, just with plsda as mode instead
                              if(input$plsda_2d3d | !input$ggplotly){
                                # 2d
                                plsda <- MetaboShiny::plotPCA.2d(mSet, lcl$aes$mycols,
                                                                 pcx = input$plsda_x,
                                                                 pcy = input$plsda_y, 
                                                                 type = "plsda",
                                                                 mode = mode,
                                                                 col.fac = input$col_var,
                                                                 shape.fac = input$shape_var,
                                                                 cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                                
                                loadings <- MetaboShiny::plotPCAloadings.2d(mSet,pcx = input$plsda_x,
                                                                            pcy = input$plsda_y, 
                                                                            type = "plsda",
                                                                            cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }else{
                                # 3d
                                plsda <- MetaboShiny::plotPCA.3d(mSet, lcl$aes$mycols,
                                                                 pcx = input$plsda_x,
                                                                 pcy = input$plsda_y,
                                                                 pcz = input$plsda_z, 
                                                                 type = "plsda",
                                                                 mode = mode,
                                                                 col.fac = input$col_var,
                                                                 shape.fac = input$shape_var,
                                                                 font = lcl$aes$font, 
                                                                 cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                                
                                loadings <- MetaboShiny::plotPCAloadings.3d(mSet,
                                                                            pcx = input$plsda_x,
                                                                            pcy = input$plsda_y,
                                                                            pcz = input$plsda_z, 
                                                                            type = "plsda",
                                                                            font = lcl$aes$font,
                                                                            cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }
                              list(plot_plsda = plsda, 
                                   plot_plsda_loadings = loadings,
                                   plsda_cv_plot = cv,
                                   plsda_perm_plot = perm)
                            }else{NULL}
                          },
                          ml = {
                            if("ml" %in% names(mSet$analSet)){
                              ml_roc = MetaboShiny::ggPlotROC(data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc,
                                                              attempts = input$ml_attempts,
                                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              
                              barplot_data <- MetaboShiny::ggPlotBar(data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$bar,
                                                                     attempts = input$ml_attempts,
                                                                     cf =gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                                     topn = input$ml_top_x,
                                                                     ml_name = mSet$analSet$ml$last$name,
                                                                     ml_type = mSet$analSet$ml$last$method)
                              
                              ml_barplot <- barplot_data$plot
                              lcl$tables$ml_bar <<- barplot_data$mzdata
                            }else{
                              shiny::hideTab(session = session, inputId = "ml2", target = "res")
                            }
                            list(ml_roc = ml_roc, 
                                 ml_bar = ml_barplot)
                          },
                          multigroup = {
                            p = MetaboShiny::ggplotSummary(mSet, my_selection$mz, 
                                                           shape.fac = input$shape_var, 
                                                           cols = lcl$aes$mycols, 
                                                           cf=gbl$functions$color.functions[[lcl$aes$spectrum]], 
                                                           mode = "multi",
                                                           styles = input$ggplot_sum_style,
                                                           add_stats = input$ggplot_sum_stats, color.fac = input$col_var, text.fac = input$txt_var)
                            list(summary_plot = p)
                          },
                          meba = {
                            p = MetaboShiny::ggplotMeba(mSet, my_selection$mz,
                                                        draw.average = T,
                                                        cols = lcl$aes$mycols,
                                                        cf=gbl$functions$color.functions[[lcl$aes$spectrum]])
                            list(summary_plot = p)
                          },
                          tt = {
                            # render manhattan-like plot for UI
                            p = MetaboShiny::ggPlotTT(mSet,
                                                      gbl$functions$color.functions[[lcl$aes$spectrum]], 
                                                      20)
                            
                            list(tt_plot = p)
                          },
                          fc = {
                            # render manhattan-like plot for UI
                            p <- MetaboShiny::ggPlotFC(mSet,
                                                       gbl$functions$color.functions[[lcl$aes$spectrum]], 20)
                            
                            list(fc_plot = p)
                          },
                          heatmap = {
                            
                            breaks = seq(min(mSet$dataSet$norm), 
                                         max(mSet$dataSet$norm), 
                                         length = 256/2)
                            
                            p = {
                              
                              if(!is.null(mSet$analSet$heatmap$matrix)){
                                # create heatmap object
                                hmap <- suppressWarnings({
                                  if(input$heatlimits){
                                    heatmaply::heatmaply(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                                         Colv = mSet$analSet$heatmap$my_order,
                                                         Rowv = T,
                                                         branches_lwd = 0.3,
                                                         margins = c(60, 0, NA, 50),
                                                         col = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                         col_side_colors = mSet$analSet$heatmap$translator[,!1],
                                                         col_side_palette = mSet$analSet$heatmap$colors,
                                                         subplot_widths = c(.9,.1),
                                                         subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                                         column_text_angle = 90,
                                                         xlab = "Sample",
                                                         ylab = "m/z",
                                                         showticklabels = c(T,F),
                                                         limits = c(min(mSet$dataSet$norm), max(mSet$dataSet$norm)),
                                                         #symm=F,symkey=F,
                                                         symbreaks=T
                                                         #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                                    )
                                  }else{
                                    heatmaply::heatmaply(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                                         Colv = mSet$analSet$heatmap$my_order,
                                                         Rowv = T,
                                                         branches_lwd = 0.3,
                                                         margins = c(60, 0, NA, 50),
                                                         colors = gbl$functions$color.functions[[lcl$aes$spectrum]](256),
                                                         col_side_colors = mSet$analSet$heatmap$translator[,!1],
                                                         col_side_palette = mSet$analSet$heatmap$colors,
                                                         subplot_widths = c(.9,.1),
                                                         subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                                         column_text_angle = 90,
                                                         xlab = "Sample",
                                                         ylab = "m/z",
                                                         #showticklabels = c(T,F),
                                                         #symm=F,symkey=F,
                                                         symbreaks=T
                                                         #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                                    )
                                    
                                    # stats::heatmap(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                    #                      Colv = mSet$analSet$heatmap$my_order,
                                    #                      Rowv = T,
                                    #                      col = gbl$functions$color.functions[[lcl$aes$spectrum]](256),
                                    #                      ColSideColors = as.character(mSet$analSet$heatmap$colors[match(mSet$analSet$heatmap$translator[,!1][[1]],
                                    #                                                                                     names(mSet$analSet$heatmap$colors)
                                    #                                                                                     )]),
                                    #                      xlab = "Sample",
                                    #                      ylab = "m/z"
                                    # )
                                    
                                  }
                                })
                                hmap$x$layout$annotations[[1]]$text <- ""
                                # save the order of mzs for later clicking functionality
                                lcl$vectors$heatmap <<- hmap$x$layout[[if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")) "yaxis2" else "yaxis3"]]$ticktext 
                                # return
                                hmap
                              }else{
                                data = data.frame(text = "No significant hits available!\nPlease try alternative source statistics below.")
                                ggplot(data) + geom_text(aes(label = text), x = 0.5, y = 0.5, size = 10) +
                                  theme(text = element_text(family = lcl$aes$font$family)) + theme_bw()
                              }
                            }
                            list(heatmap = p)
                          }, 
                          power = {
                            p = {
                              if("power" %in% names(mSet$analSet)){
                                MetaboShiny::ggPlotPower(mSet, 
                                                         max_samples = max(mSet$analSet$power[[1]]$Jpred),
                                                         comparisons = names(mSet$analSet$power),
                                                         cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }else{
                                NULL
                              }
                            }
                            list(power_plot = p)
                          },
                          wordcloud = {
                            if(nrow(lcl$tables$wordcloud_filt) > 0){
                              topWords = if(input$wordcloud_topWords > nrow(lcl$tables$wordcloud_filt)) nrow(lcl$tables$wordcloud_filt) else input$wordcloud_topWords
                              output$wordcloud <-  wordcloud2::renderWordcloud2({ 
                                wordcloud2::wordcloud2(data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,], color = "random-light", size=.7, shape = "circle")
                              })
                              p = {
                                wcdata = data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,]
                                colnames(wcdata)[2] <- "freq"
                                MetaboShiny::ggPlotWordBar(wcdata = wcdata,
                                                           cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                              }
                            }
                            list(wordbar = p)
                          }
          )
          mapply(function(myplot, plotName){
            
            if(grepl("aov|tt|fc|pattern|asca|volc", plotName)){
              whichAnal <- stringr::str_match(plotName, "aov|tt|fc|pattern|asca|volc")[,1]
              if(is.null(mSet$analSet[[whichAnal]]$sig.mat)){
                data = data.frame(text = "No significant hits!")
                myplot = ggplot(data) + geom_text(aes(label = text), x = 0.5, y = 0.5, size = 10)
              }
            }
            
            isSquare <- grepl("pca|plsda|tsne|roc|heatmap", plotName) & !grepl("scree|cv|perm|venn", plotName)
           
             # === WRAPPER ===
            empty <- if(grepl(plotName, pattern="var|samp")) "output_empty2_width" else "output_empty3_width"
         
             output[[paste0(plotName, "_wrap")]] <- shiny::renderUI({
                  list(conditionalPanel(
                    condition = 'input.ggplotly == true',
                    plotlyOutput(paste0(plotName, "_interactive"), height = "100%")),
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
                 is3D <- F
               }
               
               is3D <- plotName == "heatmap"
               
                if(!is3D){
                   myplot <- myplot + 
                     gbl$functions$plot.themes[[lcl$aes$theme]](base_size = 15) + 
                     ggplot2::theme(legend.position="none",
                                    axis.line = ggplot2::element_line(colour = 'black', size = .5),
                                    plot.title = ggplot2::element_text(hjust = 0.5,
                                                                       vjust = 0.1,
                                                                       size=lcl$aes$font$title.size*1.2,
                                                                       face="bold"),
                                    text = ggplot2::element_text(family = lcl$aes$font$family))
                   if(grepl("venn", plotName)){
                     myplot <- myplot + 
                       ggplot2::theme_void() +
                       ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                      legend.position="none",
                                      text = ggplot2::element_text(family = lcl$aes$font$family))
                        
                   }
                    
                }
                 
               if(!is.null(session$clientData[[empty]])){
                 
                 try({
                   output[[plotName]] <- shiny::renderPlot({
                     if(plotName == "heatmap"){
                       data = data.frame(text = "Currently only available in 'plotly' mode!\nPlease switch in the sidebar.")
                       ggplot(data) + geom_text(aes(label = text), x = 0.5, y = 0.5, size = 10) +
                         theme(text = element_text(family = lcl$aes$font$family)) + theme_bw()
                     }else{
                       myplot
                     }
                   })  
                 }, silent = F)
                 
                 plotFn <- paste0(c(gsub(":|,:", "_", mSet$dataSet$cls.name), 
                                    plotName), collapse="_")
                 
                 emptyax <- list(
                   title = "",
                   zeroline = FALSE,
                   showline = FALSE,
                   showticklabels = FALSE,
                   showgrid = FALSE
                 )
                 
                 output[[paste0(plotName, "_interactive")]] <- plotly::renderPlotly({
 
                   if(!is3D){
                     myplot <- plotly::ggplotly(myplot,
                                                tooltip = "text", 
                                                height = session$clientData[[empty]]/if(isSquare) 1.4 else 2)  
                   }
                   if(plotName != "heatmap"){
                     myplot <- if(grepl("venn", plotName)) myplot %>% layout(xaxis = emptyax,
                                                                             yaxis = emptyax,
                                                                             showlegend=F) else myplot %>% layout(showlegend=F) 
                   }else{
                     myplot <- myplot %>% layout(height = session$clientData[[empty]]/1.4,
                                                     width = session$clientData[[empty]])
                   }
                   myplot %>%
                     config(
                       toImageButtonOptions = list(
                         format = if(input$plotsvg) "svg" else "png",
                         filename = paste0(plotFn, "_interactive")
                       ))
                 })
                 
                 output[[paste0("download_", plotName)]] <- downloadHandler(
                   filename = function() paste0(plotFn, if(input$plotsvg) ".svg" else ".png"),
                   content = function(file){
                     ggplot2::ggsave(file, plot = myplot)
                     }
                   )
                 }
               })
            # ====
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
