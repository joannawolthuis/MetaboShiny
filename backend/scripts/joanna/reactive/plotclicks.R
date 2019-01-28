# triggers when a plotly plot is clicked by user
observeEvent(plotly::event_data("plotly_click"),{ 
  
  d <- plotly::event_data("plotly_click") # get click details (which point, additional included info, etc..)
  
  if(req(input$statistics ) %in% c("tt", "fc", "rf", "aov", "volc")){ # these cases need the same processing and use similar scoring systems
    if('key' %not in% colnames(d)) return(NULL)
    mzs <- switch(input$statistics, 
                  tt = names(mSet$analSet$tt$p.value),
                  fc = names(mSet$analSet$fc$fc.log),
                  aov = switch(input$timecourse_trigger,
                               rownames(mSet$analSet$aov2$sig.mat),
                               names(mSet$analSet$aov$p.value)),
                  volc = rownames(mSet$analSet$volcano$sig.mat)
    )
    if(d$key %not in% mzs) return(NULL)
    curr_cpd <<- d$key
    # - return -
    output[[paste0(input$statistics, "_specific_plot")]] <- plotly::renderPlotly({
      # --- ggplot ---
      ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]],
                    style = input$ggplot_sum_style, scatter = as.logical(input$ggplot_sum_scatter), add_stats = input$ggplot_sum_stats,
                    col.fac = input$col_var, txt.fac = input$txt_var)
    })
  }else if(req(input$statistics ) == "pca"){ # deprecated - used to hide and show certain groups
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
    }}else if(req(input$statistics ) == "ml"){ # makes ROC curves and boxplots clickable
      switch(input$ml_results, roc = { # if roc, check the curve numbers of the roc plot
        attempt = d$curveNumber - 1
        xvals <- mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
        if(attempt > 1){
          ml_type <- xvals$type[[1]]
          model <- xvals$models[[attempt]]
          output$ml_tab <- switch(ml_type,
                                  rf = { # random forest specific data fetching
                                    importance = as.data.table(model$importance, keep.rownames = T)
                                    rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
                                    rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)] # order descending
                                    # - - - return - - -
                                    ml_tab <<- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
                                    DT::renderDataTable({ # render importance table for selected model
                                      DT::datatable(rf_tab,
                                                    selection = 'single',
                                                    autoHideNavigation = T,
                                                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                                    })
                                  }, 
                                  ls = { # lasso specific data fetching
                                    tab = model$beta
                                    keep = which(tab[,1] != 0)
                                    tab_new = data.frame("beta" = tab[keep,1],
                                                         "absbeta" = abs(tab[keep,1]), # use the absolute beta as additional measure (min or plus importance is similar for statistical validity i think)
                                                         row.names = rownames(tab)[keep])
                                    colnames(tab_new) <- c("beta", "abs_beta")
                                    ml_tab <<- tab_new[order(tab_new[,1],decreasing = T),] # order descending
                                    DT::renderDataTable({ #  render importance table for selected model
                                      DT::datatable(ml_tab,
                                                    selection = 'single',
                                                    autoHideNavigation = T,
                                                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                                    })
                                  })
          
        }
      }, bar = { # for bar plot just grab the # bar clicked
        curr_cpd <<- as.character(mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$bar[d$x,"mz"][[1]])
      })}else if(grepl(pattern = "heatmap", x = input$statistics)){ # heatmap requires the table used to make it saved to global (hmap_mzs)
        req(global$vectors$heatmap)
        if(d$y > length(global$vectors$heatmap)) return(NULL)
        curr_cpd <<- global$vectors$heatmap[d$y]
      }
  
  # render curent miniplot based on current compound
  output$curr_plot <- plotly::renderPlotly({
    # --- ggplot ---
    ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols, 
                  cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]],
                  style = input$ggplot_sum_style, scatter = as.logical(input$ggplot_sum_scatter),
                  add_stats = input$ggplot_sum_stats,
                  col.fac = input$col_var, txt.fac = input$txt_var)
  })
  
  # change current compound in text
  output$curr_cpd <- renderText(curr_cpd)
})