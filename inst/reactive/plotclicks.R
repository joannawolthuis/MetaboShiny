# triggers when a plotly plot is clicked by user
shiny::observeEvent(plotly::event_data("plotly_click", priority = "event"), {
  
  d <<- plotly::event_data("plotly_click", priority = "event") # get click details (which point, additional included info, etc...
  print(d)
  
  for(pietype in c("add", "iso", "db")){
    try({
      if(input$tab_search == "match_filters_tab" & input$match_filters == paste0("pie_",pietype)){
        i = d$pointNumber + 1
        showsubset = as.character(pieinfo[[pietype]]$Var.1[i])
        mzMode = MetaboShiny::getIonMode(my_selection$mz, lcl$paths$patdb)
        if(pietype == "add"){
          if(!(showsubset %in% result_filters$add[[mzMode]])){
            result_filters$add[[mzMode]] <- c(result_filters$add[[mzMode]], showsubset)
          }else{
            curr_filt = result_filters$add[[mzMode]]
            result_filters$add[[mzMode]] <- curr_filt[curr_filt != showsubset]
            }
        }else{
          if(!(showsubset %in% result_filters[[pietype]])){
            result_filters[[pietype]] <- c(result_filters[[pietype]], showsubset)
          }else{
            curr_filt = result_filters[[pietype]]
            result_filters[[pietype]] <- curr_filt[curr_filt != showsubset]  
          }
        }
        search$go <- T
      }
    }, silent = F)
  }
  
  curr_tab <- switch(input$statistics,
                     dimred = {
                       input$dimred
                     }, permz = {
                       input$permz
                     }, overview = {
                       input$overview
                     }, ml = "ml")
  
  if(req(curr_tab) %in% c("tt", "pca","plsda", "fc", "rf", "aov", "volc")){ # these cases need the same processing and use similar scoring systems
    if('key' %not in% colnames(d)) return(NULL)
    mzs <- switch(curr_tab,
                  tt = names(mSet$analSet$tt$p.value),
                  fc = names(mSet$analSet$fc$fc.log),
                  plsda = rownames(mSet$analSet$plsr$loadings),
                  pca = rownames(mSet$analSet$pca$rotation),
                  aov = if(mSet$timeseries)rownames(mSet$analSet$aov2$sig.mat) else names(mSet$analSet$aov$p.value),
                  volc = rownames(mSet$analSet$volcano$sig.mat)
    )
    if(d$key %not in% mzs) return(NULL)
    my_selection$mz <- d$key
    
  }else if(req(curr_tab) == "ml"){ # makes ROC curves and boxplots clickable
    switch(input$ml_results, roc = { # if roc, check the curve numbers of the roc plot
      attempt = d$curveNumber + 1
      xvals <- mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
      if(attempt > 0){
        output$ml_tab <- DT::renderDataTable({
          imp <- data.table::as.data.table(xvals$imp[[attempt]], keep.rownames = T)
          colnames(imp) <- c("mz", "importance")
          imp <- imp[importance > 0,]
          lcl$tables$ml_roc <<- data.frame(importance = imp$importance,
                                           row.names = gsub(imp$mz,
                                                            pattern = "`|`",
                                                            replacement=""))
          DT::datatable(lcl$tables$ml_roc,
                        selection = 'single',
                        autoHideNavigation = T,
                        options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
        })
      }
    }, bar = { # for bar plot just grab the # bar clicked
      try({
        supposed_mz <- as.character(lcl$tables$ml_bar[d$x,"mz"][[1]])
        if(supposed_mz %in% colnames(mSet$dataSet$preproc)){
          my_selection$mz <- supposed_mz
        }
      })
    })}else if(req(curr_tab) == "heatmap"){#grepl(pattern = "heatmap", x = curr_tab)){ # heatmap requires the table used to make it saved to global (hmap_mzs)
      if(!is.null(d$y)){
        if(d$y > length(lcl$vectors$heatmap)) return(NULL)
        my_selection$mz <- lcl$vectors$heatmap[d$y]  
      }
    }else if(curr_tab == "pattern"){
      my_selection$mz <- rownames(mSet$analSet$corr$cor.mat)[d$curveNumber + 1]
    }else if(curr_tab == "enrich"){
      
      # TODO: make non-redundant..
      curr_pw <- rownames(mSet$analSet$enrich$mummi.resmat)[d$pointNumber + 1]
      pw_i <- which(mSet$analSet$enrich$path.nms == curr_pw)
      cpds = mSet$analSet$enrich$path.hits[[pw_i]]
      hit_tbl = data.table::as.data.table(mSet$analSet$enrich$matches.res)
      myHits <- hit_tbl[Matched.Compound %in% cpds]
      myHits$Mass.Diff <- as.numeric(myHits$Mass.Diff)/(as.numeric(myHits$Query.Mass)*1e-6)
      colnames(myHits) <- c("rn", "identifier", "adduct", "dppm")
      
      enrich$current <- myHits
      
    }
})
