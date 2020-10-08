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
            toWrap <- switch(do,
                   vennrich = {
                     # - - - - -
                     analyses = names(mSet$storage)
                     venn_no$start <- report_no$start <- data.table::rbindlist(lapply(analyses, function(name){
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
                     report_no$now <- report_no$start
                     lcl$vectors$analyses <<- unlist(venn_no$start[,1])
                     # ---
                     lapply(c("mummi_anal", "heattable", "network_table"), function(inputID){
                       shiny::updateSelectInput(session,
                                                inputID, 
                                                choices = {
                                                 allChoices = as.character(lcl$vectors$analyses)
                                                 if(inputID == "heattable" | inputID == "network_table"){
                                                   allChoices[grepl(mSet$settings$cls.name, allChoices, fixed=TRUE)]  
                                                 }else{
                                                   allChoices
                                                 }
                                                })  
                     })
                     # --- 
                     list()
                   },
                   enrich = {
                     enrich$overview <- if("mummi.resmat" %in% names(mSet$analSet$enrich)){
                       mSet$analSet$enrich$mummi.resmat 
                     }else{
                       mSet$analSet$enrich$mummi.gsea.resmat 
                     } 
                     list()
                   },
                   corr = {
                     res =if(is.null(mSet$analSet$corr$cor.mat)){
                       data.table::data.table("No significant hits found")
                     }else{
                       res = mSet$analSet$corr$cor.mat
                     }
                     list(corr_tab = res)
                   },
                   aov = {
                     which_aov = if(mSet$settings$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
                     
                     if(which_aov %in% names(mSet$analSet)){
                       
                       keep <- switch(which_aov,
                                      aov = colnames(mSet$analSet$aov$sig.mat) %in% c("p.value", "FDR", "Fisher's LSD"),
                                      aov2 = grepl("adj\\.p|Adj", colnames(mSet$analSet$aov2$sig.mat)))
                       
                       res =if(is.null(mSet$analSet[[which_aov]]$sig.mat)){
                         data.table::data.table("No significant hits found")
                       }else{
                         if(sum(keep) == 1){
                           tbl = data.table::data.table(rn=rownames(mSet$analSet[[which_aov]]$sig.mat),
                                                        "adj. p-value"=mSet$analSet[[which_aov]]$sig.mat[,keep])
                           
                         }else{
                           mSet$analSet[[which_aov]]$sig.mat[,keep]
                         }
                       }
                       list(aov_tab = res)
                     }else{
                       list()
                     }
                   },
                   volcano = {
                     # render results table
                     res <- if(is.null(mSet$analSet$volcano$sig.mat)) data.table::data.table("No significant hits found") else{
                       rownames(mSet$analSet$volcano$sig.mat) <- gsub(rownames(mSet$analSet$volcano$sig.mat), pattern = "^X", replacement = "")
                       rownames(mSet$analSet$volcano$sig.mat) <- gsub(rownames(mSet$analSet$volcano$sig.mat), pattern = "(\\d+\\.\\d+)(\\.+)", replacement = "\\1/")
                       mSet$analSet$volcano$sig.mat}
                     list(volcano_tab = res)
                   },
                   tsne = {
                     NULL
                   },
                   pca = {
                     if("pca" %in% names(mSet$analSet)){
                       # render PCA variance per PC table for UI
                       pca.table <- data.table::as.data.table(round(mSet$analSet$pca$variance * 100.00,
                                                                    digits = 2),
                                                              keep.rownames = T)
                       colnames(pca.table) <- c("Principal Component", "% variance")
                       
                       # render PCA loadings tab for UI
                       pca.loadings <- mSet$analSet$pca$rotation[,c(input$pca_x,
                                                                    input$pca_y,
                                                                    input$pca_z)]
                       list(pca_load_tab = pca.loadings,
                            pca_tab = pca.table)
                     }else{
                       list()
                     } # do nothing
                   },
                   plsda = {
                     if("plsda" %in% names(mSet$analSet)){
                       # render table with variance per PC
                       plsda.table <- data.table::as.data.table(round(mSet$analSet$plsr$Xvar
                                                                      / mSet$analSet$plsr$Xtotvar
                                                                      * 100.0,
                                                                      digits = 2),
                                                                keep.rownames = T)
                       colnames(plsda.table) <- c("Principal Component", "% variance")
                       plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
                       # render table with PLS-DA loadings
                       plsda.loadings <- mSet$analSet$plsda$vip.mat
                       colnames(plsda.loadings) <- paste0("PC", c(1:ncol(plsda.loadings)))
                       plsda.loadings = plsda.loadings[, c(input$plsda_x, input$plsda_y, input$plsda_z)]
                       list(plsda_tab = plsda.table, 
                            plsda_load_tab = plsda.loadings)
                     }else{
                       list()
                     }
                   },
                   ml = {
                     if("ml" %in% names(mSet$analSet)){
                       roc_data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
                       roc_data$perf <- data.table::as.data.table(roc_data$perf)
                       res = unique(roc_data$perf[,c("AUC_PAIR", "comparison", "attempt")])
                       lcl$tables$ml_roc_all <<- res
                       list(ml_overview_tab = res)
                     }else{
                       list()
                     }
                   },
                   asca = {
                     if("asca" %in% names(mSet$analSet)){
                       res = mSet$analSet$asca$sig.list$Model.ab
                       colnames(res) <- c("Leverage", "SPE")
                       list(asca_tab = res)
                     }else{
                       list()
                     }
                   },
                   meba = {
                     if("MB" %in% names(mSet$analSet)){
                       res = mSet$analSet$MB$stats
                       colnames(res) <- c("Hotelling/T2 score")
                       list(meba_tab = res)
                     }else{
                       list()
                     }
                   },
                   tt = {
                     # save results to table
                     res <- mSet$analSet$tt$sig.mat
                     if(is.null(res)){
                       res <- data.table::data.table("No significant hits found")
                       mSet$analSet$tt <- NULL
                     }
                     # set buttons to proper thingy
                     list(tt_tab = res)
                   },
                   fc = {
                     # if none found, give the below table...
                     # save results to table
                     res <- mSet$analSet$fc$sig.mat
                     if(is.null(res)){
                       res <- data.table::data.table("No significant hits found")
                       mSet$analSet$fc <- NULL
                     }
                     list(fc_tab = res)
                   },
                   heatmap = {
                     NULL
                   }, 
                   power = {
                     NULL
                   }
            )
          })
        }
        success = T
      })
      
      if(!success){
        metshiAlert("Table rendering failed!")
      }else{
        mapply(function(mytable, tableName){
          output[[tableName]] <- DT::renderDataTable({
            subbed = gsub("\\+", "", rownames(mytable))
            rns = rownames(mytable)
            if(subbed[1] %in% colnames(mSet$dataSet$norm)){ # check if a mz table
              stars = sapply(mSet$report$mzStarred[subbed]$star, 
                             function(hasStar) if(hasStar) "★" else "")
              starCol = data.table::data.table("★" = stars)
              mytable = cbind(starCol, mytable)
            }
            metshiTable(content = mytable, rownames = rns)
          }, server = FALSE)
        }, toWrap, names(toWrap)) 
      } 
    }
    tablemanager$make <- NULL # set makeing to 'off'
  }
})
