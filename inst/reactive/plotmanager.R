observe({
  lcl$functions$plotRender <<- if(input$ggplotly) plotly::renderPlotly else shiny::renderPlot
  lcl$functions$plotOutput <<- if(input$ggplotly) plotly::plotlyOutput else shiny::plotOutput
})

# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
plotmanager <- shiny::reactiveValues()

# preload pca/plsda
shiny::observe({
  if(is.null(plotmanager$reload)){
    NULL # if not reloading anything, nevermind
  }else{
    if(!is.null(mSet)){
      success = F
      try({
        for(do in plotmanager$make){
          suppressWarnings({
            toWrap = switch(do,
                   general = {
                     # make sidebar
                     # make pca, plsda, ml(make plotmanager do that)
                     # update select input bars with current variable and covariables defined in excel
                     if(is.null(mSet)){
                       interface$mode <- NULL
                     }else{
                       varNormPlots <- MetaboShiny::ggplotNormSummary(mSet = mSet,
                                                                      plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                                      font = lcl$aes$font,
                                                                      cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                                      plotlyfy = input$ggplotly)
                       
                       output$var1 <- lcl$functions$plotRender(varNormPlots$tl)
                       output$var2 <- lcl$functions$plotRender(varNormPlots$bl)
                       output$var3 <- lcl$functions$plotRender(varNormPlots$tr)
                       output$var4 <- lcl$functions$plotRender(varNormPlots$br)
                       sampNormPlots <- MetaboShiny::ggplotSampleNormSummary(mSet,
                                                                             plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                                             font = lcl$aes$font,
                                                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                                             plotlyfy = input$ggplotly)
                       output$samp1 <- lcl$functions$plotRender(sampNormPlots$tl)
                       output$samp2 <- lcl$functions$plotRender(sampNormPlots$bl)
                       output$samp3 <- lcl$functions$plotRender(sampNormPlots$tr)
                       output$samp4 <- lcl$functions$plotRender(sampNormPlots$br)
                       
                       return(c(paste0("var", c(1:4)), 
                                paste0("samp", c(1:4)))
                       )
                       
                     }},
                   venn = {
                     NULL
                   },
                   enrich = {
                       NULL
                     },
                   summary = {
                     output$summary_plot <- lcl$functions$plotRender({
                       MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                  styles = input$ggplot_sum_style,
                                                  add_stats = input$ggplot_sum_stats, color.fac = input$col_var, text.fac = input$txt_var,
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font, plotlyfy = input$ggplotly)
                     })
                     return("summary_plot")
                   },
                   pattern = {
                     output$pattern_plot <- lcl$functions$plotRender({
                       # --- ggplot ---
                       MetaboShiny::ggPlotPattern(mSet,n = input$pattern_topn,
                                                  cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  plotlyfy = input$ggplotly,font = lcl$aes$font)
                     })
                     return("pattern_plot")
                   },
                   aov = { # render manhattan-like plot for UI
                       output$aov_overview_plot <- lcl$functions$plotRender({
                         # --- ggplot ---
                         MetaboShiny::ggPlotAOV(mSet,
                                                cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                                plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                plotlyfy = input$ggplotly, font = lcl$aes$font)
                       })
                       return("aov_overview_plot")
                   },
                   volc = {
                     # render volcano plot with user defined colours
                     output$volc_plot <- lcl$functions$plotRender({
                       # --- ggplot ---
                       MetaboShiny::ggPlotVolc(mSet,
                                               cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                               20,
                                               plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                               plotlyfy = input$ggplotly,font = lcl$aes$font)
                     })
                     return("volc_plot")
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
                         output$plot_tsne <- lcl$functions$plotRender({
                           MetaboShiny::plotPCA.2d(mSet, 
                                                   cols = lcl$aes$mycols,
                                                   pcx = 1,
                                                   pcy = 2, 
                                                   type = "tsne",
                                                   mode = mode,
                                                   shape.fac = input$shape_var,
                                                   plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                   plotlyfy=TRUE,
                                                   col.fac = input$col_var,
                                                   font = lcl$aes$font
                                                   #,cf = gbl$functions$color.functions[[lcl$aes$spectrum]]
                           )
                         })
                       }else{
                         # render 3d plot
                         output$tsne_plot <- lcl$functions$plotRender({
                           MetaboShiny::plotPCA.3d(mSet, 
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
                         })
                       }
                     }
                     return("tsne_plot")
                   },
                   pca = {
                     if("pca" %in% names(mSet$analSet)){
                       output$pca_scree <- lcl$functions$plotRender({
                         MetaboShiny::ggPlotScree(mSet,
                                                  cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font)
                       })
                       # chekc which mode we're in
                       mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                         "ipca" # interactive PCA (old name, i like tpca more :P )
                       }else{
                         "normal" # normal pca
                       }
                       
                       if(input$pca_2d3d | !input$ggplotly){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         output$plot_pca <- lcl$functions$plotRender({
                           MetaboShiny::plotPCA.2d(mSet, 
                                                   cols = lcl$aes$mycols,
                                                   pcx = input$pca_x,
                                                   pcy = input$pca_y, 
                                                   mode = mode,
                                                   type = "pca",
                                                   col.fac = input$col_var,
                                                   shape.fac = input$shape_var,
                                                   plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                   plotlyfy=TRUE,
                                                   font = lcl$aes$font
                           )
                         })
                         output$plot_pca_loadings <- lcl$functions$plotRender({
                           MetaboShiny::plotPCAloadings.2d(mSet,pcx = input$pca_x,
                                              pcy = input$pca_y, 
                                              type = "pca",
                                              plotlyfy=TRUE,
                                              font = lcl$aes$font,
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]]
                           )
                         })
                       }else{
                         # render 3d plot
                         output$plot_pca <- lcl$functions$plotRender({
                           MetaboShiny::plotPCA.3d(mSet, 
                                                   pcx = input$pca_x,
                                                   pcy = input$pca_y,
                                                   pcz = input$pca_z, 
                                                   type = "pca",
                                                   col.fac = input$col_var,
                                                   mode = mode,
                                                   shape.fac = input$shape_var,
                                                   font = lcl$aes$font,
                                                   cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                         })
                         output$plot_pca_loadings <- lcl$functions$plotRender({
                           MetaboShiny::plotPCAloadings.3d(mSet,
                                              pcx = input$pca_x,
                                              pcy = input$pca_y,
                                              pcz = input$pca_z, 
                                              type = "pca",
                                              font = lcl$aes$font,
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                         })
                       }
                     }else{
                       NULL
                     }
                     return(c("plot_pca", "plot_pca_loadings", "pca_scree"))
                   },
                   plsda = {
                     
                     if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
                       
                       mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                         "ipca" # interactive PCA (old name, i like tpca more :P )
                       }else{
                         "normal" # normal pca
                       }
                       
                       # render cross validation plot
                       output$plsda_cv_plot <- lcl$functions$plotRender({
                         MetaboShiny::ggPlotClass(mSet, cf = gbl$functions$color.functions[[lcl$aes$spectrum]], plotlyfy = F,
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font)
                       })
                       # render permutation plot
                       output$plsda_perm_plot <- lcl$functions$plotRender({
                         MetaboShiny::ggPlotPerm(mSet,cf = gbl$functions$color.functions[[lcl$aes$spectrum]], plotlyfy = F,
                                                 plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                 font = lcl$aes$font)
                       })
                       # see PCA - render 2d or 3d plots, just with plsda as mode instead
                       if(input$plsda_2d3d | !input$ggplotly){
                         # 2d
                         output$plot_plsda <- lcl$functions$plotRender({
                           MetaboShiny::plotPCA.2d(mSet, lcl$aes$mycols,
                                                   pcx = input$plsda_x,
                                                   pcy = input$plsda_y, 
                                                   type = "plsda",
                                                   mode = mode,
                                                   col.fac = input$col_var,
                                                   shape.fac = input$shape_var,
                                                   plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                   font = lcl$aes$font,
                                                   cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                           
                         })
                         output$plot_plsda_loadings <- lcl$functions$plotRender({
                           MetaboShiny::plotPCAloadings.2d(mSet,pcx = input$plsda_x,
                                              pcy = input$plsda_y, 
                                              type = "plsda",
                                              plotlyfy=TRUE,
                                              font = lcl$aes$font,
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]]
                           )
                           
                         })
                       }else{
                         # 3d
                         output$plot_plsda <- lcl$functions$plotRender({
                           MetaboShiny::plotPCA.3d(mSet, lcl$aes$mycols,
                                                   pcx = input$plsda_x,
                                                   pcy = input$plsda_y,
                                                   pcz = input$plsda_z, 
                                                   type = "plsda",
                                                   mode = mode,
                                                   col.fac = input$col_var,
                                                   shape.fac = input$shape_var,
                                                   font = lcl$aes$font, 
                                                   cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                         })
                         output$plot_plsda_loadings <- lcl$functions$plotRender({
                           MetaboShiny::plotPCAloadings.3d(mSet,
                                              pcx = input$plsda_x,
                                              pcy = input$plsda_y,
                                              pcz = input$plsda_z, 
                                              type = "plsda",
                                              font = lcl$aes$font,
                                              cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                         })
                       }
                     }else{NULL}
                     return(c("plot_plsda", "plot_plsda_loadings"))
                   },
                   ml = {
                     if("ml" %in% names(mSet$analSet)){
                       output$ml_roc <- lcl$functions$plotRender({
                         MetaboShiny::ggPlotROC(lcl$tables$ml_roc,
                                                input$ml_attempts,
                                                gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                plotlyfy = input$ggplotly,font = lcl$aes$font)
                       })
                       barplot_data <- MetaboShiny::ggPlotBar(mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$bar,
                                                              input$ml_attempts,
                                                              gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                              input$ml_top_x,
                                                              ml_name = mSet$analSet$ml$last$name,
                                                              ml_type = mSet$analSet$ml$last$method,
                                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                              plotlyfy = input$ggplotly,font = lcl$aes$font)
                       
                       
                       ml_barplot <- barplot_data$plot
                       lcl$tables$ml_bar <<- barplot_data$mzdata
                       output$ml_bar <- lcl$functions$plotRender({
                         if(input$ggplotly){
                           plotly::ggplotly(ml_barplot)
                         }else{
                           ml_barplot 
                         }
                       })
                     }else{
                       shiny::hideTab(session = session, inputId = "ml2", target = "res")
                     }
                     return(c("ml_roc", "ml_bar"))
                   },
                   multigroup = {
                     output$summary_plot <- lcl$functions$plotRender({
                       MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, 
                                                  cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]], mode = "multi",
                                                  styles = input$ggplot_sum_style,
                                                  add_stats = input$ggplot_sum_stats, color.fac = input$col_var, text.fac = input$txt_var,
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font, plotlyfy = input$ggplotly)  
                     })
                     return("summary_plot")
                   },
                   asca = {
                     NULL
                   },
                   meba = {
                     output$meba_plot <- lcl$functions$plotRender({
                       MetaboShiny::ggplotMeba(mSet, my_selection$mz,
                                               draw.average = T,
                                               cols = lcl$aes$mycols,
                                               cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                                               plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                               font = lcl$aes$font, plotlyfy = input$ggplotly)  
                     })
                     return("meba_plot")
                   },
                   tt = {
                     # render manhattan-like plot for UI
                     output$tt_overview_plot <- lcl$functions$plotRender({
                       # --- ggplot ---
                       MetaboShiny::ggPlotTT(mSet,
                                             gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                             plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                             plotlyfy=input$ggplotly,font = lcl$aes$font)
                     })
                     return("tt_overview_plot")
                   },
                   fc = {
                     # render manhattan-like plot for UI
                     output$fc_overview_plot <- lcl$functions$plotRender({
                       # --- ggplot ---
                       MetaboShiny::ggPlotFC(mSet,
                                             gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                             plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                             plotlyfy=input$ggplotly,font = lcl$aes$font)
                     })
                     return("fc_overview_plot")
                   },
                   heatmap = {
                     
                     breaks = seq(min(mSet$dataSet$norm), max(mSet$dataSet$norm), length = 256/2)
                     
                     output$heatmap <- lcl$functions$plotRender({
                       
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
                                                  showticklabels = c(T,F),
                                                  #symm=F,symkey=F,
                                                  symbreaks=T
                                                  #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                             )
                           }
                         })
                         # save the order of mzs for later clicking functionality
                         lcl$vectors$heatmap <<- hmap$x$layout[[if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")) "yaxis2" else "yaxis3"]]$ticktext 
                         # return
                         hmap
                       }else{
                         data = data.frame(text = "No significant hits available!\nPlease try alternative source statistics below.")
                         ggplot(data) + geom_text(aes(label = text), x = 0.5, y = 0.5, size = 10) +
                           theme(text = element_text(family = lcl$aes$font$family)) + theme_bw()
                       }
                     })
                     return("heatmap_overview_plot")
                   }, 
                   power = {
                     output$power_plot <- lcl$functions$plotRender({
                       if("power" %in% names(mSet$analSet)){
                         MetaboShiny::ggPlotPower(mSet, 
                                                  max_samples = input$power_nsamp,
                                                  comparisons = input$power_comps,
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font,
                                                  cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                  plotlyfy=input$ggplotly)
                       }else{
                         NULL
                       }
                     })
                     return("power_plot")
                   },
                   wordcloud = {
                     if(nrow(lcl$tables$wordcloud_filt) > 0){
                       topWords = if(input$wordcloud_topWords > nrow(lcl$tables$wordcloud_filt)) nrow(lcl$tables$wordcloud_filt) else input$wordcloud_topWords
                       output$wordcloud <-  wordcloud2::renderWordcloud2({ 
                         wordcloud2::wordcloud2(data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,], color = "random-light", size=.7, shape = "circle")
                       })
                       output$wordbar <- lcl$functions$plotRender({
                         wcdata = data.table::as.data.table(lcl$tables$wordcloud_filt)[order(n, decreasing = T)][1:topWords,]
                         colnames(wcdata)[2] <- "freq"
                         MetaboShiny::ggPlotWordBar(wcdata = wcdata,
                                                    cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                    plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                    plotlyfy=input$ggplotly,
                                                    font = lcl$aes$font)
                       })  
                     }
                     return("wordbar")
                   }
            )
            lapply(toWrap, function(plot){
              output[[paste0(plot, "_wrap")]] <- shiny::renderUI({ 
                lcl$functions$plotOutput("plot") 
                })
            })
          })
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
