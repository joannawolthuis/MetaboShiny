# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
datamanager <- reactiveValues()

# preload pca/plsda
observe({
  
  if(is.null(datamanager$reload)){
    
    NULL # if not reloading anything, nevermind
  }else{
    if(!exists("mSet")){
      
      NULL
    }
    switch(datamanager$reload,
           
           general = {
             # change interface
             if(mSet$dataSet$cls.num <= 1){
               interface$mode <- NULL } 
             else if(mSet$dataSet$cls.num == 2){
               interface$mode <- "bivar"}
             else{
               interface$mode <- "multivar"}
             # reload sidebar
             output$curr_name <- renderText({mSet$dataSet$cls.name}) 
             # reload pca, plsda, ml(make datamanager do that)
             # update select input bars with current variable and covariables defined in excel
             updateSelectInput(session, "first_var", selected = mSet$dataSet$cls.name, choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
             updateSelectInput(session, "second_var", choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
             updateSelectInput(session, "subset_var", choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, MARGIN = 2, function(col) length(unique(col)) < global$constants$max.cols))]))
             output$curr_name <- renderText({mSet$dataSet$cls.name}) 
             # if _T in sample names, data is time series. This makes the time series swap button visible. 
             if(all(grepl(pattern = "_T\\d", x = rownames(mSet$dataSet$norm)))){
               timebutton$status <- "on"
             }else{
               timebutton$status <- "off"
             }
             # show a button with t-test or fold-change analysis if data is bivariate. hide otherwise.
             # TODO: add button for anova/other type of sorting...
             if(mSet$dataSet$cls.num == 2 ){
               heatbutton$status <- "ttfc"
             }else{
               heatbutton$status <- NULL
             }
             # tab loading
             if(mSet$dataSet$cls.num <= 1){
               interface$mode <- NULL } 
             else if(mSet$dataSet$cls.num == 2){
               interface$mode <- "bivar"}
             else{
               interface$mode <- "multivar"}
           },
           venn = {
             if("storage" %in% names(mSet)){
               analyses = names(mSet$storage)
               venn_no$start <- rbindlist(lapply(analyses, function(name){
                 analysis = mSet$storage[[name]]$analysis
                 analysis_names = names(analysis)
                 # - - -
                 with.subgroups <- intersect(analysis_names, c("ml", "plsr"))
                 if(length(with.subgroups) > 0){
                   extra_names <- lapply(with.subgroups, function(anal){
                     switch(anal,
                            ml = {
                              which.mls <- intersect(c("rf", "ls"), names(analysis$ml))
                              ml.names = sapply(which.mls, function(meth){
                                if(length(analysis$ml[[meth]]) > 0){
                                  paste0(meth, " - ", names(analysis$ml[[meth]]))
                                }
                              })
                              unlist(ml.names)
                            },
                            plsr = {
                              c ("plsda - PC1", "plsda - PC2", "plsda - PC3")
                            })
                   })
                   analysis_names <- c(setdiff(analysis_names, c("ml", "plsr", "plsda")), unlist(extra_names))
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
           aov = {
             if(!is.null(input$timecourse_trigger)){
               present = switch(input$timecourse_trigger,
                                {"aov2" %in% names(mSet$analSet)},
                                {"aov" %in% names(mSet$analSet)})
               if(present){
                 if(input$timecourse_trigger){ # send time series anova to normal anova storage
                   which.anova <- "aov2"
                   keep <- grepl("adj\\.p", colnames(mSet$analSet$aov2$sig.mat))
                 }else{
                   which.anova = "aov"
                   keep <- c("p.value", "FDR", "Fisher's LSD")
                 }
               }
             }else{
               present = "aov" %in% names(mSet$analSet)
               if(present){
                 which.anova = "aov"
                 keep <- c("p.value", "FDR", "Fisher's LSD")
               }
             }
             
             if(present){
               # render results table for UI
               output$aov_tab <- DT::renderDataTable({
                 DT::datatable(if(is.null(mSet$analSet[[which.anova]]$sig.mat)){
                   data.table("No significant hits found")
                 }else{mSet$analSet[[which.anova]]$sig.mat[,keep]
                 }, 
                 selection = 'single',
                 autoHideNavigation = T,
                 options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                })
             }
           },
           volc = {
             # render volcano plot with user defined colours
             output$volc_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotVolc(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
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
           pca = {
             if("pca" %in% req(names(mSet$analSet))){
               # create PCA legend plot
               # TODO: re-enable this plot, it was clickable so you could filter out certain groups
               output$pca_legend <- plotly::renderPlotly({
                 frame <- data.table(x = c(1), 
                                     y = mSet$dataSet$cls.num)
                 p <- ggplot(data=frame,
                             aes(x, 
                                 y, 
                                 color=factor(y),
                                 fill=factor(y)
                             )
                 ) + 
                   geom_point(shape = 21, size = 5, stroke = 5) +
                   scale_colour_manual(values=global$vectors$mycols) +
                   theme_void() + 
                   theme(legend.position="none")
                 # --- return ---
                 ggplotly(p, tooltip = NULL) %>% config(displayModeBar = F)
               })
               # render PCA variance per PC table for UI
               output$pca_tab <-DT::renderDataTable({
                 pca.table <- as.data.table(round(mSet$analSet$pca$variance * 100.00,
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
               # chekc which mode we're in
               mode <- if("timecourse_trigger" %in% names(req(input))){
                 if(input$timecourse_trigger){ # if time series mode
                   "ipca" # interactive PCA (old name, i like tpca more :P )
                 }else{
                   "pca" # normal pca
                 } 
               }else{
                 "pca"
               }
               
               # - - - - -
               if(input$pca_2d3d){ # check if switch button is in 2d or 3d mode
                 # render 2d plot
                 output$plot_pca <- plotly::renderPlotly({
                   plotPCA.2d(mSet, global$vectors$mycols,
                              pcx = input$pca_x,
                              pcy = input$pca_y, mode = mode,
                              shape.fac = input$second_var)
                 })
               }else{
                 # render 3d plot
                 output$plot_pca <- plotly::renderPlotly({
                   plotPCA.3d(mSet, global$vectors$mycols,
                              pcx = input$pca_x,
                              pcy = input$pca_y,
                              pcz = input$pca_z, mode = mode,
                              shape.fac = input$second_var)
                 })
               }
             }else{
               
             } # do nothing
           },
           plsda = {
             
             if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
               
               # render cross validation plot
               output$plsda_cv_plot <- renderPlot({
                 ggPlotClass(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], plotlyfy = F)
               })
               # render permutation plot
               output$plsda_perm_plot <- renderPlot({
                 ggPlotPerm(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]], plotlyfy = F)
               })
               # render table with variance per PC
               output$plsda_tab <- DT::renderDataTable({
                 # - - - -
                 plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
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
                   plotPCA.2d(mSet, global$vectors$mycols,
                              pcx = input$plsda_x,
                              pcy = input$plsda_y, mode = "plsda",
                              shape.fac = input$second_var)
                 })
               }else{
                 # 3d
                 output$plot_plsda <- plotly::renderPlotly({
                   plotPCA.3d(mSet, global$vectors$mycols,
                              pcx = input$plsda_x,
                              pcy = input$plsda_y,
                              pcz = input$plsda_z, mode = "plsda",
                              shape.fac = input$second_var)
                 })
               }
             }else{NULL}
           },
           ml = {
             if("ml" %in% names(mSet$analSet)){
               roc_data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
               
               output$ml_roc <- plotly::renderPlotly({
                 plotly::ggplotly(ggPlotROC(roc_data, 
                                            input$ml_attempts, 
                                            global$functions$color.functions[[getOptions("user_options.txt")$gspec]]))
               })
               
               bar_data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$bar
               
               output$ml_bar <- plotly::renderPlotly({
                 
                 plotly::ggplotly(ggPlotBar(bar_data, 
                                            input$ml_attempts, 
                                            global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 
                                            input$ml_top_x, 
                                            ml_name = mSet$analSet$ml$last$name,
                                            ml_type = mSet$analSet$ml$last$method))
               })
             }else{NULL}
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
               res <<- data.table("No significant hits found")
               mSet$analSet$tt <<- NULL
             }
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
               ggPlotTT(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           },
           fc = {
             # save results table
             res <<- mSet$analSet$fc$sig.mat 
             # if none found, give the below table...
             if(is.null(res)) res <<- data.table("No significant hits found")
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
               ggPlotFC(global$functions$color.functions[[getOptions("user_options.txt")$gspec]], 20)
             })
           },
           heatmap = {
             
             breaks = seq(min(mSet$dataSet$norm), max(mSet$dataSet$norm), length = 256/2)
             
              output$heatmap <- plotly::renderPlotly({
                 
                 if(!is.null(mSet$analSet$heatmap$matrix)){
                   # create heatmap object
                   hmap <- suppressWarnings({
                     if(input$heatmap_scaleall){
                       heatmaply::heatmaply(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                            Colv = mSet$analSet$heatmap$my_order, 
                                            Rowv = T,
                                            branches_lwd = 0.3,
                                            margins = c(60, 0, NA, 50),
                                            col = global$functions$color.functions[[getOptions("user_options.txt")$gspec]],
                                            col_side_colors = mSet$analSet$heatmap$translator[,!1],
                                            col_side_palette = mSet$analSet$heatmap$colors,
                                            subplot_widths = c(.9,.1),
                                            subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                            column_text_angle = 90,
                                            xlab = "Sample",
                                            ylab = "m/z",
                                            showticklabels = c(T,F),
                                            limits = c(min(mSet$dataSet$norm), max(mSet$dataSet$norm))
                                            #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                       )  
                     }else{
                       heatmaply::heatmaply(mSet$analSet$heatmap$matrix[1:if(input$heatmap_topn < nrow(mSet$analSet$heatmap$matrix)) input$heatmap_topn else nrow(mSet$analSet$heatmap$matrix),],
                                            Colv = mSet$analSet$heatmap$my_order, 
                                            Rowv = T,
                                            branches_lwd = 0.3,
                                            margins = c(60, 0, NA, 50),
                                            colors = global$functions$color.functions[[getOptions("user_options.txt")$gspec]](256),
                                            col_side_colors = mSet$analSet$heatmap$translator[,!1],
                                            col_side_palette = mSet$analSet$heatmap$colors,
                                            subplot_widths = c(.9,.1),
                                            subplot_heights = if(mSet$analSet$heatmap$my_order) c(.1, .05, .85) else c(.05,.95),
                                            column_text_angle = 90,
                                            xlab = "Sample",
                                            ylab = "m/z",
                                            showticklabels = c(T,F)
                                            #label_names = c("m/z", "sample", "intensity") #breaks side colours...
                       )
                     }
                     
                  })
                   
                   # save the order of mzs for later clicking functionality
                   global$vectors$heatmap <<- hmap$x$layout$yaxis3$ticktext
                   
                   # return
                   hmap
                 }else{
                   data = data.frame(text = "No significant hits available!\nPlease try alternative source statistics below.")
                   p <- ggplot(data) + geom_text(aes(label = text), x = 0.5, y = 0.5, size = 10) +
                     theme(text = element_text(family = global$constants$font.aes$font)) + theme_bw()
                   plotly::ggplotly(p)
                 }
               })
           })
    # - - - - 
    datamanager$reload <- NULL # set reloading to 'off'
  }
})