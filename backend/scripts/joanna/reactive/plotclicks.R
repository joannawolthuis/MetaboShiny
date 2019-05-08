# triggers when a plotly plot is clicked by user
observeEvent(plotly::event_data("plotly_click"),{

  d <- plotly::event_data("plotly_click") # get click details (which point, additional included info, etc..)

  if(input$tab_iden_4 == "pie_add"){

    i = d$pointNumber + 1
    show.adduct = local$vectors$pie_add[i]

    if(nrow(local$tables$last_matches)>0){
      keep.rows <- which(local$tables$last_matches$adduct == show.adduct)
      shown_matches$table <- local$tables$last_matches[keep.rows,]
    }
    return(NULL)
  }

  if(input$tab_iden_4 == "pie_db"){
    i = d$pointNumber + 1
    show.db = local$vectors$pie_db[i]
    if(nrow(local$tables$last_matches)>0){
      keep.rows <- which(local$tables$last_matches$source == show.db)
      shown_matches$table <- local$tables$last_matches[keep.rows,]
    }
    return(NULL)
  }

  curr_tab <- switch(input$statistics,
                     dimred = {
                       input$dimred
                     }, permz = {
                       input$permz
                     }, overview = {
                       input$overview
                     })
  
  print(curr_tab)
  
  if(req(curr_tab ) %in% c("tt", "fc", "rf", "aov", "volc")){ # these cases need the same processing and use similar scoring systems
    if('key' %not in% colnames(d)) return(NULL)
    mzs <- switch(curr_tab,
                  tt = names(mSet$analSet$tt$p.value),
                  fc = names(mSet$analSet$fc$fc.log),
                  aov = switch(input$timecourse_trigger,
                               rownames(mSet$analSet$aov2$sig.mat),
                               names(mSet$analSet$aov$p.value)),
                  volc = rownames(mSet$analSet$volcano$sig.mat)
    )
    if(d$key %not in% mzs) return(NULL)
    local$curr_mz <<- d$key
    # - return -
    output[[paste0(curr_tab, "_specific_plot")]] <- plotly::renderPlotly({
      # --- ggplot ---
      ggplotSummary(mSet, local$curr_mz, shape.fac = input$shape_var, cols = local$aes$mycols,cf=global$functions$color.functions[[local$aes$spectrum]],
                    styles = input$ggplot_sum_style, add_stats = input$ggplot_sum_stats,
                    col.fac = input$col_var, txt.fac = input$txt_var,
                    plot.theme = global$functions$plot.themes[[local$aes$theme]],
                    font = local$aes$font)
    })
  }else if(req(curr_tab) == "pca"){ # deprecated - used to hide and show certain groups
    if(!"z" %in% names(d)){
      which_group = 1
      which_group <- d$curveNumber + 1
      traceLoc <- length(unique(mSet$dataSet$cls)) + 1
      scatter <- pca_plot$x$attrs[[traceLoc]]
      idx <- scatter$color == which_group
      if(pca_plot$x$data[[which_group]]$visible == "legendonly"){
        pca_plot$x$data[[which_group]]$visible = TRUE
        scatter$visible[idx] <- T
      }else{ # hide
        pca_plot$x$data[[which_group]]$visible = "legendonly"
        scatter$visible[idx] <- F
      }
      pca_plot$x$attrs[[traceLoc]] <- scatter
      pca_plot <<- pca_plot
      output$plot_pca <- plotly::renderPlotly({pca_plot})
    }}else if(req(curr_tab) == "ml"){ # makes ROC curves and boxplots clickable
      switch(input$ml_results, roc = { # if roc, check the curve numbers of the roc plot
        attempt = d$curveNumber - 1
        xvals <- mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
        if(attempt > 1){
          ml_type <- xvals$type[[1]]
          model <- xvals$models[[attempt]]

          output$ml_tab <- DT::renderDataTable({
            imp <- as.data.table(caret::varImp(model)$importance, keep.rownames = T)
            colnames(imp) <- c("mz", "importance")
            imp <- imp[importance > 0,]
            local$tables$ml_roc <<- data.frame(importance = imp$importance, 
                                                row.names = gsub(imp$mz, 
                                                                 pattern = "`|`", 
                                                                 replacement=""))
            DT::datatable(local$tables$ml_roc,
                          selection = 'single',
                          autoHideNavigation = T,
                          options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
          })
        }
      }, bar = { # for bar plot just grab the # bar clicked
        local$curr_mz <<- as.character(local$tables$ml_bar[d$x,"mz"][[1]])
      })}else if(grepl(pattern = "heatmap", x = curr_tab)){ # heatmap requires the table used to make it saved to global (hmap_mzs)
        req(local$vectors$heatmap)
        if(d$y > length(local$vectors$heatmap)) return(NULL)
        local$curr_mz <<- local$vectors$heatmap[d$y]
      }

  # render curent miniplot based on current compound
  output$curr_plot <- plotly::renderPlotly({
    # --- ggplot ---
    ggplotSummary(mSet, local$curr_mz, shape.fac = input$shape_var, cols = local$aes$mycols,
                  cf=global$functions$color.functions[[local$aes$spectrum]],
                  styles = input$ggplot_sum_style,
                  add_stats = input$ggplot_sum_stats,
                  col.fac = input$col_var, txt.fac = input$txt_var,plot.theme = global$functions$plot.themes[[local$aes$theme]],
                  font = local$aes$font)
  })

  # change current compound in text
  output$curr_mz <- renderText(local$curr_mz)
})
