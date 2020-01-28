# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
datamanager <- shiny::reactiveValues()


# preload pca/plsda
shiny::observe({
  if(is.null(datamanager$reload)){
    NULL # if not reloading anything, nevermind
  }else{
    if(!is.null(mSet)){
      mSet.old <- mSet
      success = F
      try({
        for(do in datamanager$reload){
          suppressWarnings({
            switch(do,
                   general = {
                     # reload sidebar
                     # reload pca, plsda, ml(make datamanager do that)
                     # update select input bars with current variable and covariables defined in excel
                     if(is.null(mSet)){
                       interface$mode <- NULL
                     }else{
                       varNormPlots <- MetaboShiny::ggplotNormSummary(mSet = mSet,
                                                                      plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                                      font = lcl$aes$font,
                                                                      cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       
                       output$var1 <- shiny::renderPlot(varNormPlots$tl)
                       output$var2 <- shiny::renderPlot(varNormPlots$bl)
                       output$var3 <- shiny::renderPlot(varNormPlots$tr)
                       output$var4 <- shiny::renderPlot(varNormPlots$br)
                       sampNormPlots <- MetaboShiny::ggplotSampleNormSummary(mSet,
                                                                             plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                                             font = lcl$aes$font,
                                                                             cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       output$samp1 <- shiny::renderPlot(sampNormPlots$tl)
                       output$samp2 <- shiny::renderPlot(sampNormPlots$bl)
                       output$samp3 <- shiny::renderPlot(sampNormPlots$tr)
                       output$samp4 <- shiny::renderPlot(sampNormPlots$br)
                       
                       if(is.null(mSet$dataSet$exp.type)){
                         mSet$dataSet$exp.type <- "1f" # one factor, binary class
                       }  
                       
                       shiny::showNotification("Updating interface...")
                       datamanager$reload <- "statspicker"
                       interface$mode <<- mSet$dataSet$exp.type
                       output$curr_name <- shiny::renderText({mSet$dataSet$cls.name})
                       shiny::updateNavbarPage(session, "statistics", selected = "inf")
                       origcount = nrow(mSet$storage$orig$data$norm)
                       output$samp_count <- shiny::renderText({
                         paste0(as.character(nrow(mSet$dataSet$norm)),
                                if(nrow(mSet$dataSet$norm) == origcount) "" else paste0("/",as.character(origcount)))
                         
                       })
                       shiny::updateSelectInput(session, "storage_choice", 
                                                choices = setdiff(names(mSet$storage)[sapply(mSet$storage, function(x) "settings" %in% names(x))], 'orig')
                       )
                     }
                     
                     shiny::updateSelectInput(session, "stats_type", 
                                              selected = {
                                                if(grepl(mSet$dataSet$exp.type, pattern = "^1f\\w")){
                                                  gsub(mSet$dataSet$exp.type, pattern="^1f\\w", "1f")
                                                }else{
                                                  mSet$dataSet$exp.type
                                                }
                                              })
                     
                     if(!any(duplicated(mSet$dataSet$covars$individual))){
                       shiny::updateSelectInput(session, "stats_type", 
                                                choices = list("one factor"="1f", 
                                                               "two factors"="2f")) 
                     }else{
                       shiny::updateSelectInput(session, "stats_type", 
                                                choices = list("one factor"="1f", 
                                                               "two factors"="2f",
                                                               "time series"="t",
                                                               "time series + one factor"="t1f"))
                       shiny::updateSelectInput(session, "time_var", selected = mSet$dataSet$time.var, 
                                                choices = c("label", 
                                                            colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]))
                     }
                     shiny::updateSelectInput(session, "stats_var", selected = mSet$dataSet$exp.var, 
                                              choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]))
                     shiny::updateSelectInput(session, "shape_var", 
                                              selected = "label",
                                              choices = c("label", 
                                                          colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]))
                     shiny::updateSelectInput(session, "col_var", 
                                              selected = "label", 
                                              choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]))
                     shiny::updateSelectInput(session, "txt_var", 
                                              selected = "sample", 
                                              choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]))
                     shiny::updateSelectInput(session, "subset_var", 
                                              choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]))
                     shiny::updateSelectInput(session, "ml_include_covars", 
                                              choices = c(colnames(mSet$dataSet$covars)[!(colnames(mSet$dataSet$covars) %in% c("label", "sample", "individual"))]))
                     
                     wordcloud_filters = file.path(lcl$paths$work_dir, "wordcloud")
                     if(dir.exists(wordcloud_filters)){
                       filter_files = list.files(wordcloud_filters, full.names = T)
                       items = lapply(filter_files, function(file){
                         data.table::fread(file)$word
                       })
                       itemNames = gsub(pattern = "\\.csv", replacement = "", x = basename(filter_files))
                       names(items) <- itemNames
                       gbl$vectors$wordcloud$filters <<- append(gbl$vectors$wordcloud$filters, items)
                       shiny::updateSelectInput(session, inputId = "wordcloud_filter", choices = c("stopwords","metabolomics","default", itemNames))
                     }
                     
                     shiny::setProgress(0.7)
                     
                     if(mSet$metshiParams$prematched){
                       search_button$on <- FALSE
                     }else{
                       search_button$on <- TRUE
                     }
                     if("tt" %in% names(mSet$analSet)){
                       if("V" %in% colnames(mSet$analSet$tt$sig.mat)){
                         shiny::updateCheckboxInput(session, "tt_nonpar", value = T)
                       }else{
                         shiny::updateCheckboxInput(session, "tt_nonpar", value = F)
                       }
                     }else{
                       shiny::updateCheckboxInput(session, "tt_nonpar", value = F)
                     }
                     
                     if(grepl(mSet$dataSet$exp.type, pattern = "1f.")){
                       pairs = MetaboShiny::expand.grid.unique(levels(mSet$dataSet$cls), 
                                                               levels(mSet$dataSet$cls))
                       pairs = paste(pairs[,1], "vs.", pairs[,2])
                       shiny::updateSelectInput(session, inputId = "power_comps", choices = pairs, selected = pairs[1])
                     }
                     
                     switch(mSet$dataSet$exp.type,
                            "1fb"=shinyWidgets::updateRadioGroupButtons(session, "heattable", choices = list("T-test"="tt", 
                                                                                                             "Fold-change analysis"="fc"), 
                                                                        selected = "tt"),
                            "1fm"=shinyWidgets::updateRadioGroupButtons(session, "heattable", choices = c(ANOVA="aov"), selected = "aov"),
                            "2f"=shinyWidgets::updateRadioGroupButtons(session, "heattable", choices = list(ANOVA="aov2", 
                                                                                                            ASCA="asca"), selected = "aov2"),
                            "t1f"=shinyWidgets::updateRadioGroupButtons(session, "heattable", choices = list(ANOVA="aov2", 
                                                                                                             ASCA="asca",
                                                                                                             MEBA="meba"), selected = "aov2"),
                            "t"=shinyWidgets::updateRadioGroupButtons(session, "heattable", choices = list(ANOVA="aov2",
                                                                                                           MEBA="meba"), selected = "aov2"))
                   },
                   statspicker = {
                     output$stats_picker <- shiny::renderUI({
                       if(input$stats_type == "2f"){
                         shiny::selectizeInput("stats_var", 
                                               label="Experimental variables:", 
                                               multiple = T,
                                               selected = mSet$dataSet$exp.var, 
                                               choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, 
                                                                                                              MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]),
                                               options = list(maxItems = 2),
                                               width = "80%") 
                       }else if(input$stats_type == "t"){
                         list()
                       }else{
                         shiny::selectizeInput("stats_var", 
                                               label="Experimental variable:", 
                                               selected = mSet$dataSet$exp.var, 
                                               choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, 
                                                                                                              MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]),
                                               multiple = F,
                                               width = "80%")
                       }
                     })  
                   },
                   venn = {
                     if("storage" %in% names(mSet)){
                       # save previous mset
                       mset_name = mSet$dataSet$cls.name
                       # TODO: use this in venn diagram creation
                       mSet$storage[[mset_name]] <- list(analysis = mSet$analSet)
                       # - - - - -
                       analyses = names(mSet$storage)
                       venn_no$start <- rbindlist(lapply(analyses, function(name){
                         analysis = mSet$storage[[name]]$analysis
                         analysis_names = names(analysis)
                         # - - -
                         exclude = c("tsne", "heatmap", "type")
                         analysis_names <- setdiff(analysis_names, exclude)
                         if(length(analysis_names) == 0){
                           return(data.table::data.table())
                         }
                         # - - -
                         with.subgroups <- intersect(analysis_names, c("ml", "plsr", "pca"))
                         if(length(with.subgroups) > 0){
                           extra_names <- lapply(with.subgroups, function(anal){
                             switch(anal,
                                    ml = {
                                      which.mls <- setdiff(names(analysis$ml),"last")
                                      ml.names = sapply(which.mls, function(meth){
                                        if(length(analysis$ml[[meth]]) > 0){
                                          paste0(meth, " - ", names(analysis$ml[[meth]]))
                                        }
                                      })
                                      unlist(ml.names)
                                    },
                                    plsr = {
                                      c ("plsda - PC1", "plsda - PC2", "plsda - PC3")
                                    },
                                    pca = {
                                      c ("pca - PC1", "pca - PC2", "pca - PC3")
                                    })
                           })
                           analysis_names <- c(setdiff(analysis_names, c("ml", "plsr", "plsda", "pca")), unlist(extra_names))
                         }
                         # - - -
                         data.frame(
                           paste0(analysis_names, " (", name, ")")
                         )
                       }))
                       venn_no$now <- venn_no$start
                     }else{
                       venn_no$start <- data.frame(names(mSet$analSet))
                       venn_no$now <- venn_no$start
                     }
                   },
                   pattern = {
                     output$pattern_plot <- plotly::renderPlotly({
                       # --- ggplot ---
                       MetaboShiny::ggPlotPattern(mSet,n = input$pattern_topn,
                                                  cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  plotlyfy=TRUE,font = lcl$aes$font)
                     })
                     # render results table
                     output$pattern_tab <-DT::renderDataTable({
                       # -------------
                       DT::datatable(mSet$analSet$corr$cor.mat,
                                     selection = 'single',
                                     autoHideNavigation = T,
                                     options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                     })
                   },
                   aov = {
                     which_aov = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
                     
                     if(which_aov %in% names(mSet$analSet)){
                       
                       keep <- switch(which_aov,
                                      aov = colnames(mSet$analSet$aov$sig.mat) %in% c("p.value", "FDR", "Fisher's LSD"),
                                      aov2 = grepl("adj\\.p|Adj", colnames(mSet$analSet$aov2$sig.mat)))
                       
                       output$aov_tab <- DT::renderDataTable({
                         DT::datatable(if(is.null(mSet$analSet[[which_aov]]$sig.mat)){
                           data.table::data.table("No significant hits found")
                         }else{
                           if(sum(keep) == 1){
                             tbl = data.table::data.table(rn=rownames(mSet$analSet[[which_aov]]$sig.mat),
                                                          "adj. p-value"=mSet$analSet[[which_aov]]$sig.mat[,keep])
                             
                           }else{
                             mSet$analSet[[which_aov]]$sig.mat[,keep]
                           }
                         },
                         selection = 'single',
                         autoHideNavigation = T,
                         options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       })
                       
                       # render manhattan-like plot for UI
                       output$aov_overview_plot <- plotly::renderPlotly({
                         # --- ggplot ---
                         MetaboShiny::ggPlotAOV(mSet,
                                                cf = gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                                plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                plotlyfy=TRUE, font = lcl$aes$font)
                       })
                     }
                   },
                   volc = {
                     # render volcano plot with user defined colours
                     output$volc_plot <- plotly::renderPlotly({
                       # --- ggplot ---
                       MetaboShiny::ggPlotVolc(mSet,
                                               cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                                               20,
                                               plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                               plotlyfy=TRUE,font = lcl$aes$font)
                     })
                     # render results table
                     output$volc_tab <-DT::renderDataTable({
                       # -------------
                       DT::datatable(mSet$analSet$volc$sig.mat,
                                     selection = 'single',
                                     autoHideNavigation = T,
                                     options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                     })
                   },
                   tsne = {
                     mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                       "ipca" # interactive PCA (old name, i like tpca more :P )
                     }else{
                       "normal" # normal pca
                     }
                     
                     if("tsne" %in% names(mSet$analSet)){
                       if(input$tsne_2d3d){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         output$plot_tsne <- plotly::renderPlotly({
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
                         output$plot_tsne <- plotly::renderPlotly({
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
                   },
                   pca = {
                     if("pca" %in% names(mSet$analSet)){
                       # render PCA variance per PC table for UI
                       output$pca_tab <-DT::renderDataTable({
                         pca.table <- data.table::as.data.table(round(mSet$analSet$pca$variance * 100.00,
                                                                      digits = 2),
                                                                keep.rownames = T)
                         colnames(pca.table) <- c("Principal Component", "% variance")
                         
                         DT::datatable(pca.table,
                                       selection = 'single',
                                       autoHideNavigation = T,
                                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       })
                       # render PCA loadings tab for UI
                       output$pca_load_tab <-DT::renderDataTable({
                         pca.loadings <- mSet$analSet$pca$rotation[,c(input$pca_x,
                                                                      input$pca_y,
                                                                      input$pca_z)]
                         #colnames(pca.loadings)[1] <- "m/z"
                         DT::datatable(pca.loadings,
                                       selection = 'single',
                                       autoHideNavigation = T,
                                       options = list(lengthMenu = c(5, 10, 15), pageLength = 10))
                       })
                       output$pca_scree <- renderPlot({
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
                       
                       if(input$pca_2d3d){ # check if switch button is in 2d or 3d mode
                         # render 2d plot
                         output$plot_pca <- plotly::renderPlotly({
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
                                                   #,cf = gbl$functions$color.functions[[lcl$aes$spectrum]]
                           )
                         })
                       }else{
                         # render 3d plot
                         output$plot_pca <- plotly::renderPlotly({
                           MetaboShiny::plotPCA.3d(mSet, 
                                                   lcl$aes$mycols,
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
                       }
                     }else{
                       NULL
                     } # do nothing
                   },
                   plsda = {
                     
                     if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
                       
                       mode <- if(mSet$dataSet$exp.type %in% c("2f", "t", "t1f")){ # if time series mode
                         "ipca" # interactive PCA (old name, i like tpca more :P )
                       }else{
                         "normal" # normal pca
                       }
                       
                       # render cross validation plot
                       output$plsda_cv_plot <- renderPlot({
                         MetaboShiny::ggPlotClass(mSet, cf = gbl$functions$color.functions[[lcl$aes$spectrum]], plotlyfy = F,
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font)
                       })
                       # render permutation plot
                       output$plsda_perm_plot <- renderPlot({
                         MetaboShiny::ggPlotPerm(mSet,cf = gbl$functions$color.functions[[lcl$aes$spectrum]], plotlyfy = F,
                                                 plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                 font = lcl$aes$font)
                       })
                       # render table with variance per PC
                       output$plsda_tab <- DT::renderDataTable({
                         # - - - -
                         plsda.table <- data.table::as.data.table(round(mSet$analSet$plsr$Xvar
                                                                        / mSet$analSet$plsr$Xtotvar
                                                                        * 100.0,
                                                                        digits = 2),
                                                                  keep.rownames = T)
                         colnames(plsda.table) <- c("Principal Component", "% variance")
                         plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
                         # -------------
                         DT::datatable(plsda.table,
                                       selection = 'single',
                                       autoHideNavigation = T,
                                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       })
                       # render table with PLS-DA loadings
                       output$plsda_load_tab <-DT::renderDataTable({
                         plsda.loadings <- mSet$analSet$plsda$vip.mat
                         colnames(plsda.loadings) <- paste0("PC", c(1:ncol(plsda.loadings)))
                         # -------------
                         DT::datatable(plsda.loadings[, c(input$plsda_x, input$plsda_y, input$plsda_z)],
                                       selection = 'single',
                                       autoHideNavigation = T,
                                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       })
                       # see PCA - render 2d or 3d plots, just with plsda as mode instead
                       if(input$plsda_2d3d){
                         # 2d
                         output$plot_plsda <- plotly::renderPlotly({
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
                       }else{
                         # 3d
                         output$plot_plsda <- plotly::renderPlotly({
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
                       }
                     }else{NULL}
                   },
                   ml = {
                     if("ml" %in% names(mSet$analSet)){
                       
                       shiny::showTab(session = session, 
                                      inputId = "ml2", 
                                      target = "res")
                       
                       roc_data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
                       
                       output$ml_roc <- plotly::renderPlotly({
                         plotly::ggplotly(MetaboShiny::ggPlotROC(roc_data,
                                                                 input$ml_attempts,
                                                                 gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                                 plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                                 plotlyfy=TRUE,font = lcl$aes$font))
                       })
                       
                       bar_data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$bar
                       
                       barplot_data <- MetaboShiny::ggPlotBar(bar_data,
                                                              input$ml_attempts,
                                                              gbl$functions$color.functions[[lcl$aes$spectrum]],
                                                              input$ml_top_x,
                                                              ml_name = mSet$analSet$ml$last$name,
                                                              ml_type = mSet$analSet$ml$last$method,
                                                              plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                              plotlyfy=TRUE,font = lcl$aes$font)
                       
                       
                       ml_barplot <- barplot_data$plot
                       lcl$tables$ml_bar <<- barplot_data$mzdata
                       
                       output$ml_bar <- plotly::renderPlotly({
                         
                         plotly::ggplotly(ml_barplot)
                         
                       })
                       choices = c()
                       methods <- setdiff(names(mSet$analSet$ml), "last")
                       for(method in methods){
                         model.names = names(mSet$analSet$ml[[method]])
                         choices <- c(choices, paste0(method, " - ", paste0(model.names)))
                       }
                       shiny::updateSelectInput(session, "show_which_ml", choices = choices, selected = paste0(mSet$analSet$ml$last$method, " - ", mSet$analSet$ml$last$name))
                       
                     }else{
                       shiny::hideTab(session = session, inputId = "ml2", target = "res")
                     }
                   },
                   asca = {
                     if("asca" %in% names(mSet$analSet)){
                       output$asca_tab <-DT::renderDataTable({ # render results table for UI
                         # -------------
                         DT::datatable(mSet$analSet$asca$sig.list$Model.ab,
                                       selection = 'single',
                                       colnames = c("Compound", "Leverage", "SPE"),
                                       autoHideNavigation = T,
                                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       })
                     }
                   },
                   meba = {
                     if("MB" %in% names(mSet$analSet)){
                       output$meba_tab <-DT::renderDataTable({
                         # -------------
                         DT::datatable(mSet$analSet$MB$stats,
                                       selection = 'single',
                                       colnames = c("Compound", "Hotelling/T2 score"),
                                       autoHideNavigation = T,
                                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       })
                     }
                   },
                   tt = {
                     # save results to table
                     res <<- mSet$analSet$tt$sig.mat
                     
                     if(is.null(res)){
                       res <<- data.table::data.table("No significant hits found")
                       mSet$analSet$tt <- NULL
                     }
                     # set buttons to proper thingy
                     # render results table for UI
                     output$tt_tab <-DT::renderDataTable({
                       # -------------
                       DT::datatable(res,
                                     selection = 'single',
                                     autoHideNavigation = T,
                                     options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       
                     })
                     # render manhattan-like plot for UI
                     output$tt_overview_plot <- plotly::renderPlotly({
                       # --- ggplot ---
                       MetaboShiny::ggPlotTT(mSet,
                                             gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                             plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                             plotlyfy=TRUE,font = lcl$aes$font)
                     })
                   },
                   fc = {
                     # save results table
                     res <<- mSet$analSet$fc$sig.mat
                     # if none found, give the below table...
                     if(is.null(res)) res <<- data.table::data.table("No significant hits found")
                     # render result table for UI
                     output$fc_tab <-DT::renderDataTable({
                       # -------------
                       DT::datatable(res,
                                     selection = 'single',
                                     autoHideNavigation = T,
                                     options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                       
                     })
                     # render manhattan-like plot for UI
                     output$fc_overview_plot <- plotly::renderPlotly({
                       # --- ggplot ---
                       MetaboShiny::ggPlotFC(mSet,
                                             gbl$functions$color.functions[[lcl$aes$spectrum]], 20,
                                             plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                             plotlyfy=TRUE,font = lcl$aes$font)
                     })
                   },
                   heatmap = {
                     
                     breaks = seq(min(mSet$dataSet$norm), max(mSet$dataSet$norm), length = 256/2)
                     
                     output$heatmap <- plotly::renderPlotly({
                       
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
                         p <- ggplot(data) + geom_text(aes(label = text), x = 0.5, y = 0.5, size = 10) +
                           theme(text = element_text(family = lcl$aes$font$family)) + theme_bw()
                         plotly::ggplotly(p)
                       }
                     })
                   }, 
                   power = {
                     output$power_plot <- plotly::renderPlotly({
                       if("power" %in% names(mSet$analSet)){
                         MetaboShiny::ggPlotPower(mSet, 
                                                  max_samples = input$power_nsamp,
                                                  comparisons = input$power_comps,
                                                  plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                                                  font = lcl$aes$font,
                                                  cf = gbl$functions$color.functions[[lcl$aes$spectrum]])
                       }else{
                         
                         NULL
                       }
                     })
                     
                   },
                   wordcloud = {
                     output$wordcloud <-  wordcloud2::renderWordcloud2({ 
                       if(nrow(lcl$tables$wordcloud_filt) > 0){
                         topWords = if(input$wordcloud_topWords > nrow(lcl$tables$wordcloud_filt)) nrow(lcl$tables$wordcloud_filt) else input$wordcloud_topWords
                         wordcloud2::wordcloud2(lcl$tables$wordcloud_filt[order(n, decreasing = T)][1:topWords,], color = "random-light", size=.7, shape = "circle")
                       }
                     })
                     # TODO: fix barchart
                     # output$wordbar_desc <- renderPlotly({MetaboShiny::ggPlotWordBar(wcdata = wcdata,
                     #                                                                 cf = gbl$functions$color.functions[[lcl$aes$spectrum]],
                     #                                                                 plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                     #                                                                 plotlyfy = TRUE, 
                     #                                                                 font = lcl$aes$font)})
                   }
            )
          })
        }
        success = T
      })
      if(success){
        mSet <<- mSet
      }else{
        MetaboShiny::metshiAlert("Data plotting failed!")
        mSet <<- mSet.old
      }
    }
    datamanager$reload <- NULL # set reloading to 'off'
  }
})
