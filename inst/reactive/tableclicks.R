# check which adducts are currently selected by user
shiny::observe({
  wanted.adducts <- lcl$vectors$calc_adducts[input$magicball_add_tab_rows_selected]
  # ---------
  lcl$vectors$add_list <<- wanted.adducts
})

# creates observers for click events in the tables defined above
lapply(c("tt",
         "fc",
         "aov",
         "rf",
         "asca",
         "meba",
         "pca_load",
         "plsda_load",
         "enrich_pw",
         "ml",
         "mummi_detail",
         "venn"), FUN=function(table){
  shiny::observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
    curr_row = input[[paste0(table, "_tab_rows_selected")]]
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column

    if (is.null(curr_row)) return()

    res_tbl <- data.table::as.data.table(switch(table,
                                                tt = mSet$analSet$tt$sig.mat,
                                                fc = mSet$analSet$fc$sig.mat,
                                                pca_load = mSet$analSet$pca$rotation,
                                                plsda_load = mSet$analSet$plsda$vip.mat,
                                                ml = lcl$tables$ml_roc, #TODO: fix this, now in global
                                                asca = mSet$analSet$asca$sig.list$Model.ab,
                                                aov = {
                                                  if(!is.null(input$timecourse_trigger)){
                                                    switch(input$timecourse_trigger,
                                                           mSet$analSet$aov2$sig.mat,
                                                           mSet$analSet$aov$sig.mat)
                                                  }else{
                                                    mSet$analSet$aov$sig.mat
                                                  }
                                                },
                                                rf = vip.score,
                                                enrich_pw = enrich_overview_tab,
                                                meba = mSet$analSet$MB$stats,
                                                plsda_vip = plsda_tab,
                                                mummi_detail = lcl$tables$mummi_detail,
                                                venn = lcl$tables$venn_overlap), keep.rownames = T)
    if(nrow(res_tbl) > 0){

      # get current selected compound from the original table (needs to be available in global env)
      my_selection$mz <- res_tbl[curr_row, rn]
      
      outplot_name <- paste0(table, "_specific_plot")

      # send plot to relevant spot in UI
      output[[outplot_name]] <- plotly::renderPlotly({
        # --- ggplot ---
        if(table == 'meba'){ # meba needs a split by time
          MetaboShiny::ggplotMeba(mSet, my_selection$mz,
                     draw.average = T,
                     cols = lcl$aes$mycols,
                     cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                     plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                     font = lcl$aes$font
                     )
        }else if(table == 'asca'){ # asca needs a split by time
          MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, 
                        cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]], mode = "ts",
                        styles = input$ggplot_sum_style,
                        add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                        plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                        font = lcl$aes$font)
        }else{ # regular boxplot
          if(!is.null(input$timecourse_trigger)){
            if(input$timecourse_trigger){
              MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]], mode = "ts",
                            styles = input$ggplot_sum_style,
                            add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                            plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                            font = lcl$aes$font)
            }else{
              MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                            styles = input$ggplot_sum_style,
                            add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                            plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                            font = lcl$aes$font)
            }
          }else{
            MetaboShiny::ggplotSummary(mSet, my_selection$mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                          styles = input$ggplot_sum_style,
                          add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                          plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                          font = lcl$aes$font)
          }
        }
      })
    }
  })
})

# mummichog enrichment
lapply(c("pos", "neg"), function(mode){
  shiny::observeEvent(input[[paste0("mummi_", mode, "_tab_rows_selected")]],{
    curr_row <- input[[paste0("mummi_", mode, "_tab_rows_selected")]] # get current row
    if (is.null(curr_row)) return()
    curr_pw <- rownames(lcl$vectors[[paste0("mummi_", mode)]]$sig)[curr_row]
    cpds <- lcl$vectors[[paste0("mummi_", mode)]]$pw2cpd[[curr_pw]]
    mzs <- lcl$vectors[[paste0("mummi_", mode)]]$cpd2mz[cpds]
    keep <- sapply(mzs, function(x) !is.null(x))
    mzs <- mzs[keep]
    mzs <- unique(unlist(mzs))
    # - - - - - - - -
    tbl = data.frame("p-value" = if(mSet$dataSet$cls.num == 2){
      mSet$analSet$tt$sig.mat[match(mzs, rownames(mSet$analSet$tt$sig.mat)),"p.value"]
    }else{
      mSet$analSet$aov$sig.mat[match(mzs, rownames(mSet$analSet$tt$sig.mat)),"p.value"]
    })
    rownames(tbl) <- mzs
    # - - - - - - - -
    lcl$tables$mummi_detail <<- tbl
    # - - - - - - - -
    output$mummi_detail_tab <- DT::renderDataTable({
      DT::datatable(tbl, selection = 'single')
    })
  })
})

# triggers on clicking a row in the match results table
shiny::observeEvent(input$match_tab_rows_selected,{
  curr_row <- input$match_tab_rows_selected # get current row
  if (is.null(curr_row)) return()
  # - - - - - - - - - - - - - - - - - - - - - -
  my_selection$name <- shown_matches$forward_unique[curr_row,'name'][[1]] # get current structure
  updateTextInput(session,
                  "pm_query",
                   value = my_selection$name)
  my_selection$form <- unlist(shown_matches$forward_unique[curr_row,'baseformula']) # get current formula
  my_selection$struct <- unlist(shown_matches$forward_unique[curr_row,'structure']) # get current formula
})


shiny::observeEvent(input$browse_tab_rows_selected,{
  curr_row <- input$browse_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_def <- browse_content$table[curr_row, description]
  output$browse_definition <- shiny::renderText(curr_def)
  my_selection$revstruct <- browse_content$table[curr_row,c('structure')][[1]]
})

# triggers on clicking a row in the reverse hit results table
shiny::observeEvent(input$hits_tab_rows_selected,{
  curr_row <- input$hits_tab_rows_selected # get current row
  if (is.null(curr_row)) return()
  my_selection$mz <- shown_matches$reverse[curr_row,'query_mz'][[1]]
})