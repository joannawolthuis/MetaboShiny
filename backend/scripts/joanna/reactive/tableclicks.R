# check which adducts are currently selected by user
observe({
  wanted.adducts <- lcl$vectors$calc_adducts[input$magicball_add_tab_rows_selected]
  # ---------
  lcl$vectors$add_list <<- wanted.adducts
})


# which table names to check for user click events
res.update.tables <<- c("tt",
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
                        "venn")

# creates observers for click events in the tables defined above
lapply(unique(res.update.tables), FUN=function(table){
  observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
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
      lcl$curr_mz <<- res_tbl[curr_row, rn]
      
      outplot_name <- paste0(table, "_specific_plot")

      # send plot to relevant spot in UI
      output[[outplot_name]] <- plotly::renderPlotly({
        # --- ggplot ---
        if(table == 'meba'){ # meba needs a split by time
          ggplotMeba(mSet, lcl$curr_mz,
                     draw.average = T,
                     cols = lcl$aes$mycols,
                     cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                     plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                     font = lcl$aes$font
                     )
        }else if(table == 'asca'){ # asca needs a split by time
          ggplotSummary(mSet, lcl$curr_mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]], mode = "ts",
                        styles = input$ggplot_sum_style,
                        add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                        plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                        font = lcl$aes$font)
        }else{ # regular boxplot
          if(!is.null(input$timecourse_trigger)){
            if(input$timecourse_trigger){
              ggplotSummary(mSet, lcl$curr_mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]], mode = "ts",
                            styles = input$ggplot_sum_style,
                            add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                            plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                            font = lcl$aes$font)
            }else{
              ggplotSummary(mSet, lcl$curr_mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                            styles = input$ggplot_sum_style,
                            add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                            plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                            font = lcl$aes$font)
            }
          }else{
            ggplotSummary(mSet, lcl$curr_mz, shape.fac = input$shape_var, cols = lcl$aes$mycols, cf=gbl$functions$color.functions[[lcl$aes$spectrum]],
                          styles = input$ggplot_sum_style,
                          add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var,
                          plot.theme = gbl$functions$plot.themes[[lcl$aes$theme]],
                          font = lcl$aes$font)
          }

        }
      })
      
      datamanager$reload <- "mz"
      
    }

  })
})

# mummichog enrichment
lapply(c("pos", "neg"), function(mode){
  observeEvent(input[[paste0("mummi_", mode, "_tab_rows_selected")]],{
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
observeEvent(input$match_tab_rows_selected,{
  curr_row <<- input$match_tab_rows_selected # get current row
  if (is.null(curr_row)) return()
  # - - - - - - - - - - - - - - - - - - - - - -
  try({
    curr_name <<- shown_matches$forward[curr_row,'name'][[1]]
    updateTextInput(session, 
                    "pm_query",
                    value = curr_name)
    # - - - - - - - - - - - - - - - - - - - - - -
    curr_def <<- shown_matches$forward[curr_row,'description'][[1]] # get current definition (hidden in table display but not deleted)
    output$curr_definition <- renderText(curr_def) # render definition
    curr_struct <<- shown_matches$forward[curr_row,'structure'][[1]] # get current structure
    output$curr_struct <- renderPlot({plot.mol(curr_struct,style = "cow")}) # plot molecular structure
    curr_formula <<- shown_matches$forward[curr_row,'baseformula'][[1]] # get current formula
    output$curr_formula <- renderText({curr_formula}) # render text of current formula
  })
})


observeEvent(input$browse_tab_rows_selected,{
  curr_row <- input$browse_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_def <- lcl$tables$browse_table[curr_row, description]
  output$browse_definition <- renderText(curr_def)
  
  search_cmd <- lcl$tables$browse_table[curr_row,c('structure')][[1]]
  
  if(!mSet$metshiParams$prematched){
    print("Please perform pre-matching first to enable this feature!")
    return(NULL)
  }else{
    lcl$tables$hits_table <<- unique(get_prematches(who = search_cmd,
                                                    what = "map.structure", #map.mz as alternative
                                                    patdb = lcl$paths$patdb)[,c("query_mz", "adduct", "%iso", "dppm")])
    shown_matches$reverse <- if(nrow(lcl$tables$hits_table) > 0){
      lcl$tables$hits_table
    }else{
      data.table('name' = "Didn't find anything ( •́ .̫ •̀ )")
    }
  }
})

# triggers on clicking a row in the reverse hit results table
observeEvent(input$hits_tab_rows_selected,{
  curr_row <<- input$hits_tab_rows_selected # get current row
  if (is.null(curr_row)) return()
  # - - - - - - - - - - - - - - - - - - - - - -
  try({
    lcl$curr_mz <<- shown_matches$reverse[curr_row,'query_mz'][[1]]
    datamanager$reload <- "mz"
    
  })
})