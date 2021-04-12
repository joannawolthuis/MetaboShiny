# triggers when a plotly plot is clicked by user
shiny::observeEvent(plotly::event_data("plotly_click", priority = "event"), {
  
  d <<- plotly::event_data("plotly_click", priority = "event") # get click details (which point, additional included info, etc...
  
  for(pietype in c("add", "iso", "db")){
    try({
      if(input$tab_search == "match_filters_tab" & input$match_filters == paste0("pie_",pietype)){
        i = d$pointNumber + 1
        showsubset = as.character(pieinfo[[pietype]]$Var.1[i])
        mzMode =if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
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
  
  curr_tab <- input$statistics
  
  if(curr_tab %in% c("tt", 
                     "pca",
                     "heatmap",
                     "ml",
                     "plsda", 
                     "fc", 
                     "rf", 
                     "aov", 
                     "corr",
                     "volcano",
                     "network",
                     "enrich",
                     "venn",
                     "meba",
                     "combi",
                     "asca")){ 
    
    if(curr_tab == "ml" & input$ml_results == "roc"){
      attempt = as.numeric(d$key[[1]])
      if(length(attempt) > 0 & !is.na(attempt)){
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
            MetaboShiny:: metshiTable(lcl$tables$ml_roc)
          })
        }
      }
    }else{
      if(curr_tab %in% c("heatmap", "enrich","venn","network")){
        switch(curr_tab,
               heatmap = {
                 if(!is.null(d$y)){
                   if(d$y > length(lcl$vectors$heatmap)) return(NULL)
                   my_selection$mz <<- lcl$vectors$heatmap[d$y]  
                   plotmanager$make <- "summary"
                 }
               },
               network = {
                 if(!is.null(d$y)){
                   if(d$y > length(lcl$vectors$network_heatmap)) return(NULL)
                   my_selection$mz <<- lcl$vectors$network_heatmap[d$y]  
                   plotmanager$make <- "summary"
                 }
               },
               enrich = {
                 curr_pw <- rownames(enrich$overview)[d$pointNumber + 1]
                 pw_i <- which(mSet$analSet$enrich$path.nms == curr_pw)
                 cpds = unlist(mSet$analSet$enrich$path.hits[[pw_i]])
                 hit_tbl = data.table::as.data.table(mSet$analSet$enrich$mumResTable)
                 myHits <- hit_tbl[Matched.Compound %in% unlist(cpds)]
                 myHits$Mass.Diff <- as.numeric(myHits$Mass.Diff)/(as.numeric(myHits$Query.Mass)*1e-6)
                 colnames(myHits) <- c("rn", "identifier", "adduct", "dppm")
                 enrich$current <- myHits
               },
               venn = {
                 if("key" %in% colnames(d)){
                   picked_intersection = d$key[[1]]
                   groups = stringr::str_split(picked_intersection, "<br />")[[1]]
                   shiny::updateSelectInput(session = session, 
                                            inputId = "intersect_venn", 
                                            selected = groups)   
                 }
               })
      }else{
        if('key' %not in% colnames(d)) return(NULL)
        #if(gsub(d$key[[1]],pattern="`",replacement="") %not in% colnames(mSet$dataSet$proc)) return(NULL)
        my_selection$mz <<- gsub(d$key[[1]],pattern="`",replacement="")
        plotmanager$make <- "summary"
      }
    }
  }
})

shiny::observeEvent(input$network_interactive_selected, {
  my_selection$mz <<- input$network_interactive_selected
  plotmanager$make <- "summary"
})
