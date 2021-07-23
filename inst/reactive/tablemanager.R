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
          #suppressWarnings({
            toWrap <- switch(do,
                   vennrich = {
                     # - - - - -
                     analyses = names(mSet$storage) 
                     analyses_table = data.table::rbindlist(lapply(analyses, function(name){
                       analysis = mSet$storage[[name]]$analSet
                       analysis_names = names(analysis)
                       # - - -
                       exclude = c("tsne", "heatmap", "type", "enrich", "power", "network")
                       analysis_names <- setdiff(analysis_names, exclude)
                       if(length(analysis_names) == 0){
                         return(data.table::data.table())
                       }
                       # - - -
                       with.subgroups <- intersect(analysis_names, c("ml", "plsr", "pca"))
                       extra_names <- if(length(with.subgroups) > 0){
                         lapply(with.subgroups, function(anal){
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
                                    c ("plsda - Component 1", "plsda - Component 2", "plsda - Component 3")
                                  },
                                  pca = {
                                    c ("pca - PC1", "pca - PC2", "pca - PC3")
                                  })
                         })
                       }else{ list() }
                       analysis_names <- c(setdiff(analysis_names, c("ml", "plsr", "plsda", "pca")), unlist(extra_names), "all m/z")
                       
                       # - - -
                       data.frame(
                         name = paste0(analysis_names, " (", name, ")"),
                         threshold = c("any")
                       )
                     }))
                     #if(ncol(venn_no$start) > 0){
                     
                     lcl$vectors$analyses <<- unlist(analyses_table[,1])
                     #}else{
                     #  lcl$vectors$analyses <<- c()
                     #}
                     # ---
                     lapply(c("mummi_anal", "heattable", "network_table", "ml_specific_mzs"), function(inputID){
                       shiny::updateSelectizeInput(session,
                                                inputID, 
                                                choices = {
                                                 ch = allChoices = as.character(lcl$vectors$analyses)
                                                 if(inputID %in% c("heattable", "network_table")){
                                                   ch = allChoices[grepl(mSet$settings$cls.name, allChoices, fixed=TRUE)]  
                                                   }else{
                                                   ch = allChoices
                                                   if(inputID == "ml_specific_mzs"){
                                                     ch = c("no", "manual", ch)
                                                   }
                                                   }
                                                 ch
                                                })#, server = T) 
                     })
                     # --- 
                     venn_no$start <- report_no$start <- analyses_table
                     venn_no$now <- venn_no$start
                     report_no$now <- report_no$start
                     
                     list()
                   },
                   enrich = {
                     enrich$overview <<- if(!is.null(mSet$analSet$enrich$"mummi.resmat")){
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
                   featsel = {
                     decision = mSet$analSet$featsel[[1]]$finalDecision
                     res = data.frame(decision = decision[decision != "Rejected"], 
                                      row.names = names(decision[decision != "Rejected"]))
                     list(featsel_tab = res) 
                   },
                   pca = {
                     if("pca" %in% names(mSet$analSet)){
                       # render PCA variance per PC table for UI
                       pca.table <- data.table::as.data.table(round(mSet$analSet$pca$variance * 100.00,
                                                                    digits = 2),
                                                              keep.rownames = T)
                       colnames(pca.table) <- c("Principal Component", "% variance")
                       
                       # render PCA loadings tab for UI
                       pca.loadings <- mSet$analSet$pca$rotation[,as.numeric(c(input$pca_x,
                                                                               input$pca_y,
                                                                               input$pca_z))]
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
                       colnames(plsda.table) <- c("Component", "% variance")
                       plsda.table[, "Component"] <- paste0("Component ", 1:nrow(plsda.table))
                       # render table with PLS-DA loadings
                       plsda.loadings <- mSet$analSet$plsda$vip.mat
                       colnames(plsda.loadings) <- paste0("Component ", c(1:ncol(plsda.loadings)))
                       plsda.loadings = plsda.loadings[, as.numeric(c(input$plsda_x, input$plsda_y, input$plsda_z))]
                       list(plsda_tab = plsda.table, 
                            plsda_load_tab = plsda.loadings)
                     }else{
                       list()
                     }
                   },
                   ml = {
                     ###
                     
                     if("ml" %in% names(mSet$analSet)){
                       data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]
                       if(!is.null(data$res$prediction)){
                         data$res$shuffled = FALSE
                         data$res = list(data$res)
                       }
                       
                       ml_performance_rows = lapply(1:length(data$res), function(i){
                         res = data$res[[i]]
                         ml_performance = getMLperformance(res, 
                                                           pos.class = input$ml_plot_posclass,
                                                           x.metric=input$ml_plot_x,
                                                           y.metric=input$ml_plot_y)
                         ml_performance$coords$shuffled = c(res$shuffled)
                         ml_performance$coords$run = i
                         ml_performance
                       })
                       coords = data.table::rbindlist(lapply(ml_performance_rows, function(x) x$coords))
                       ml_performance = list(coords = coords,
                                             names = ml_performance_rows[[1]]$names)
                       
                       split_coords = split(coords, coords$run)
                       
                       ml_tbl_rows = lapply(split_coords, function(tbl){
                         shuffled = unique(tbl$shuffled)
                         training.rows = data.table::data.table()
                         if(!shuffled){
                           training = tbl[grep("Fold", `Test set`)]
                           spl.folds = split(training, training$`Test set`)
                           fold.rows = lapply(spl.folds, function(test_set){
                             data.table::data.table(AUC = pracma::trapz(test_set$x, 
                                                                        test_set$y),
                                                    "Test set" = unique(test_set$`Test set`))
                           })
                           training.rows = data.table::rbindlist(fold.rows)
                         }
                         # testing
                         test_name = if(shuffled) paste0("Shuffled test #", unique(tbl$run) - 1) else "Test"
                         test_set = tbl[`Test set`=="Test"]
                         testing.row = data.table::data.table(AUC = pracma::trapz(test_set$x, 
                                                                                  test_set$y),
                                                              "Test set" = test_name)
                         rbind(training.rows, testing.row)
                       })
                       res = data.table::rbindlist(ml_tbl_rows)
                       lcl$tables$ml_roc_all <<- res
                       # params 
                       if("params" %in% names(data)){
                         params = data.table::data.table(param = gsub("ml_", "", names(data$params)), value = data$params)
                       }else{
                         params = data.table::data.table(unavailable = "No parameters saved for this model...")
                       }
                       
                       no_shuffle = data$res[[which(unlist(sapply(data$res, function(x) !x$shuffle)))]]
                       res2 = no_shuffle$importance
                       rownames(res2) = gsub("^X", "", rownames(res2))
                       rownames(res2) = gsub("\\.$", "-", rownames(res2))
                       res2 = data.frame(importance=res2[,1], row.names=rownames(res2))
                       lcl$tables$ml_imp <<- res2
                       
                       list(ml_overview_tab = res,
                            ml_param_tab = params,
                            ml_importance_tab = res2)
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
                   combi = {
                     # save results to table
                     res <- mSet$analSet$combi$sig.mat
                     if(is.null(res)){
                       res <- data.table::data.table("No significant hits found")
                       mSet$analSet$combi <- NULL
                     }
                     res = as.data.frame(res)
                     rownames(res) <- res$rn
                     res$rn <- NULL
                     # set buttons to proper thingy
                     list(combi_tab = res)
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
          #})
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
              # starico = "★"
              starico = '<i class=\"fa fa-star\" role=\"presentation\" aria-label=\"star icon\"></i>'
              stars = sapply(mSet$report$mzStarred[subbed]$star, 
                             function(hasStar) if(hasStar) starico else "")
              starCol = data.table::data.table(starico = stars)
              colnames(starCol) = starico
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
