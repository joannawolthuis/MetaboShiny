# check which adducts are currently selected by user
observe({
  # --------------
  wanted.adducts <- global$vectors$calc_adducts[input$magicball_add_tab_rows_selected]
  # ---------
  global$vectors$add_list <<- wanted.adducts
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
    # get current selected compound from the original table (needs to be available in global env)
    curr_cpd <<- data.table::as.data.table(switch(table,
                                                  tt = mSet$analSet$tt$sig.mat,
                                                  fc = mSet$analSet$fc$sig.mat,
                                                  pca_load = mSet$analSet$pca$rotation,
                                                  plsda_load = mSet$analSet$plsda$vip.mat,
                                                  ml = global$tables$ml_roc, #TODO: fix this, now in global
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
                                                  mummi_detail = global$tables$mummi_detail,
                                                  venn = global$tables$venn_overlap), keep.rownames = T)[curr_row, rn]
    # print current compound in sidebar
    output$curr_cpd <- renderText(curr_cpd)
    
    # make miniplot for sidebar with current compound
    output$curr_plot <- plotly::renderPlotly({
      # --- ggplot ---
      ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]],
                    styles = input$ggplot_sum_style,
                    add_stats = input$ggplot_sum_stats, col.fac = input$col_var,txt.fac = input$txt_var)
    })
    
    outplot_name <- paste0(table, "_specific_plot")
    
    # send plot to relevant spot in UI
    output[[outplot_name]] <- plotly::renderPlotly({
      # --- ggplot ---
      if(table == 'meba'){ # meba needs a split by time
        ggplotMeba(curr_cpd, draw.average = T, cols = global$vectors$mycols,cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]])
      }else if(table == 'asca'){ # asca needs a split by time
        ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]], mode = "ts",
                      styles = input$ggplot_sum_style,
                      add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var)
      }else{ # regular boxplot
        if(!is.null(input$timecourse_trigger)){
          if(input$timecourse_trigger){
            ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]], mode = "ts",
                          styles = input$ggplot_sum_style,
                          add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var)
          }else{
            ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]],
                          styles = input$ggplot_sum_style,
                          add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var)
          }
        }else{
          ggplotSummary(curr_cpd, shape.fac = input$shape_var, cols = global$vectors$mycols, cf=global$functions$color.functions[[getOptions("user_options.txt")$gspec]],
                        styles = input$ggplot_sum_style,
                        add_stats = input$ggplot_sum_stats, col.fac = input$col_var, txt.fac = input$txt_var)
        }
        
      }
    })
  })
})

# mummichog enrichment
lapply(c("pos", "neg"), function(mode){
  observeEvent(input[[paste0("mummi_", mode, "_tab_rows_selected")]],{
    curr_row <- input[[paste0("mummi_", mode, "_tab_rows_selected")]] # get current row
    if (is.null(curr_row)) return()
    curr_pw <- rownames(global$vectors[[paste0("mummi_", mode)]]$sig)[curr_row]
    cpds <- global$vectors[[paste0("mummi_", mode)]]$pw2cpd[[curr_pw]]
    mzs <- global$vectors[[paste0("mummi_", mode)]]$cpd2mz[cpds]
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
    global$tables$mummi_detail <<- tbl
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
  try({
    curr_name <<- global$tables$last_matches[curr_row,'name'][[1]]
    updateTextInput(session, "pm_query", value = curr_name)
    # write to clipboard
    if(input$auto_copy){
      clipr::write_clip(curr_name)
      print('copied to clipboard ( ˘ ³˘)♥')
    }
    # -----------------------------
    curr_def <<- global$tables$last_matches[curr_row,'description'] # get current definition (hidden in table display but not deleted)
    output$curr_definition <- renderText(curr_def$description) # render definition
    curr_struct <<- global$tables$last_matches[curr_row,'structure'][[1]] # get current structure
    output$curr_struct <- renderPlot({plot.mol(curr_struct,style = "cow")}) # plot molecular structure
    curr_formula <<- global$tables$last_matches[curr_row,'baseformula'][[1]] # get current formula
    output$curr_formula <- renderText({curr_formula}) # render text of current formula
  })
})