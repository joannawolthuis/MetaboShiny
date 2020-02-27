# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
tablemanager <- shiny::reactiveValues()

# pmake pca/plsda
shiny::observe({
  if(is.null(tablemanager$make)){
    NULL # if not makeing anything, nevermind
  }else{
    if(!is.null(mSet)){
      success = F
      try({
        for(do in tablemanager$make){
          suppressWarnings({
            switch(do,
                   vennrich = {
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
                         exclude = c("tsne", "heatmap", "type", "enrich")
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
                       lcl$vectors$analyses <<- analysis_names
                     }else{
                       venn_no$start <- data.frame(names(mSet$analSet))
                       venn_no$now <- venn_no$start
                     }
                   },
                   pattern = {
                     # render results table
                     output$pattern_tab <-DT::renderDataTable({
                       MetaboShiny::metshiTable(content = mSet$analSet$corr$cor.mat)
                     })
                   },
                   aov = {
                     which_aov = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
                     
                     if(which_aov %in% names(mSet$analSet)){
                       
                       keep <- switch(which_aov,
                                      aov = colnames(mSet$analSet$aov$sig.mat) %in% c("p.value", "FDR", "Fisher's LSD"),
                                      aov2 = grepl("adj\\.p|Adj", colnames(mSet$analSet$aov2$sig.mat)))
                       
                       con =if(is.null(mSet$analSet[[which_aov]]$sig.mat)){
                         data.table::data.table("No significant hits found")
                       }else{
                         if(sum(keep) == 1){
                           tbl = data.table::data.table(rn=rownames(mSet$analSet[[which_aov]]$sig.mat),
                                                        "adj. p-value"=mSet$analSet[[which_aov]]$sig.mat[,keep])
                           
                         }else{
                           mSet$analSet[[which_aov]]$sig.mat[,keep]
                         }
                       }
                       
                       output$aov_tab <- DT::renderDataTable({
                         MetaboShiny::metshiTable(content = con)
                       })
                     }
                   },
                   volc = {
                     # render results table
                     output$volc_tab <-DT::renderDataTable({
                       res <- if(is.null(mSet$analSet$volc$sig.mat)) data.table::data.table("No significant hits found") else{
                         rownames(mSet$analSet$volc$sig.mat) <<- gsub(rownames(mSet$analSet$volc$sig.mat), pattern = "^X", replacement = "")
                         rownames(mSet$analSet$volc$sig.mat) <<- gsub(rownames(mSet$analSet$volc$sig.mat), pattern = "(\\d+\\.\\d+)(\\.+)", replacement = "\\1/")
                         res = mSet$analSet$volc$sig.mat
                       }
                       MetaboShiny::metshiTable(content = res)
                     })
                   },
                   tsne = {
                     NULL
                     },
                   pca = {
                     if("pca" %in% names(mSet$analSet)){
                       # render PCA variance per PC table for UI
                       output$pca_tab <-DT::renderDataTable({
                         pca.table <- data.table::as.data.table(round(mSet$analSet$pca$variance * 100.00,
                                                                      digits = 2),
                                                                keep.rownames = T)
                         colnames(pca.table) <- c("Principal Component", "% variance")
                         MetaboShiny::metshiTable(content = pca.table)
                       })
                       # render PCA loadings tab for UI
                       output$pca_load_tab <-DT::renderDataTable({
                         pca.loadings <- mSet$analSet$pca$rotation[,c(input$pca_x,
                                                                      input$pca_y,
                                                                      input$pca_z)]
                         MetaboShiny::metshiTable(content = pca.loadings)
                       })
                     }else{
                       NULL
                     } # do nothing
                   },
                   plsda = {
                     
                     if("plsda" %in% names(mSet$analSet)){ # if plsda has been performed...
                       
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
                         MetaboShiny::metshiTable(content = plsda.table)
                       })
                       # render table with PLS-DA loadings
                       output$plsda_load_tab <-DT::renderDataTable({
                         plsda.loadings <- mSet$analSet$plsda$vip.mat
                         colnames(plsda.loadings) <- paste0("PC", c(1:ncol(plsda.loadings)))
                         MetaboShiny::metshiTable(content = plsda.loadings[, c(input$plsda_x, input$plsda_y, input$plsda_z)])
                       })
                     }else{NULL}
                   },
                   ml = {
                     if("ml" %in% names(mSet$analSet)){
                       roc_data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
                       roc_data$perf <- data.table::as.data.table(roc_data$perf)
                       
                       tbl = unique(roc_data$perf[,c("AUC_PAIR", "comparison", "attempt")])
                       
                       output$ml_overview_tab <- DT::renderDataTable({
                         MetaboShiny::metshiTable(content = tbl)
                       })                  
                     }
                   },
                   asca = {
                     if("asca" %in% names(mSet$analSet)){
                       output$asca_tab <-DT::renderDataTable({ 
                         con = mSet$analSet$asca$sig.list$Model.ab
                         colnames(con) <- c("Leverage", "SPE")
                         MetaboShiny::metshiTable(content = con)
                       })
                     }
                   },
                   meba = {
                     if("MB" %in% names(mSet$analSet)){
                       output$meba_tab <-DT::renderDataTable({
                         con = mSet$analSet$MB$stats
                         colnames(con) <- c("Hotelling/T2 score")
                         MetaboShiny::metshiTable(content = con)
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
                       MetaboShiny::metshiTable(content = res)
                     })
                   },
                   fc = {
                     # if none found, give the below table...
                     res = mSet$analSet$fc$sig.mat
                     res <- if(is.null(mSet$analSet$fc$sig.mat)) data.table::data.table("No significant hits found")
                     
                     # render result table for UI
                     output$fc_tab <-DT::renderDataTable({
                       MetaboShiny::metshiTable(content = res)
                     })
                   },
                   heatmap = {
                     NULL
                   }, 
                   power = {
                     NULL
                   },
                   wordcloud = {
                     NULL
                   }
            )
          })
        }
        success = T
      })
      if(!success){
        MetaboShiny::metshiAlert("Table rendering failed!")
      }
    }
    tablemanager$make <- NULL # set makeing to 'off'
  }
})
