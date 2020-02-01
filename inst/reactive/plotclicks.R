# triggers when a plotly plot is clicked by user
shiny::observeEvent(plotly::event_data("plotly_click"), {
  
  d <- plotly::event_data("plotly_click") # get click details (which point, additional included info, etc...

  try({
    for(pietype in c("add", "iso", "db")){
      if(input$tab_iden_4 == paste0("pie_",pietype) | input$tab_iden_5 == paste0("pie_",pietype) ){
        i = d$pointNumber + 1
        showsubset = as.character(pieinfo[[pietype]]$Var.1[i])
        result_filters[[pietype]] <- showsubset
        search$go <- T
      }
    }
  }, silent=T)
  
  curr_tab <- switch(input$statistics,
                     dimred = {
                       input$dimred
                     }, permz = {
                       input$permz
                     }, overview = {
                       input$overview
                     }, ml = "ml")
  
  if(req(curr_tab) %in% c("tt", "fc", "rf", "aov", "volc")){ # these cases need the same processing and use similar scoring systems
    if('key' %not in% colnames(d)) return(NULL)
    mzs <- switch(curr_tab,
                  tt = names(mSet$analSet$tt$p.value),
                  fc = names(mSet$analSet$fc$fc.log),
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
    }
})
