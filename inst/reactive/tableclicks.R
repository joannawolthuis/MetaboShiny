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
         "volc",
         "asca",
         "meba",
         "pca_load",
         "plsda_load",
         "enrich_pw",
         "ml",
         "pattern",
         "mummi_detail",
         "venn"), FUN=function(table){

  shiny::observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
    curr_row = input[[paste0(table, "_tab_rows_selected")]]
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    
    lcl$curr_table <<- table
    
    if (is.null(curr_row)) return()

    which_aov = if(mSet$dataSet$exp.type %in% c("t", "2f", "t1f")) "aov2" else "aov"
    
    res_tbl <- data.table::as.data.table(switch(table,
                                                pattern = mSet$analSet$corr$cor.mat,
                                                tt = mSet$analSet$tt$sig.mat,
                                                fc = mSet$analSet$fc$sig.mat,
                                                volc = mSet$analSet$volcano$sig.mat,
                                                pca_load = mSet$analSet$pca$rotation,
                                                plsda_load = mSet$analSet$plsda$vip.mat,
                                                ml = lcl$tables$ml_roc, #TODO: fix this, now in global
                                                asca = mSet$analSet$asca$sig.list$Model.ab,
                                                aov = mSet$analSet[[which_aov]]$sig.mat,
                                                enrich_pw = enrich$current,
                                                rf = vip.score,
                                                meba = mSet$analSet$MB$stats,
                                                plsda_vip = plsda_tab,
                                                mummi_detail = lcl$tables$mummi_detail,
                                                venn = lcl$tables$venn_overlap),
                                         keep.rownames = T)
    if(nrow(res_tbl) > 0){
      # get current selected compound from the original table (needs to be available in global env)
      my_selection$mz <- res_tbl[curr_row, rn]
      # --- ggplot ---
      if(table == 'meba'){ # meba needs a split by time
        plotmanager$make <- "meba"
      }else if(table %in% c('aov2', 'asca')){ # asca needs a split by time
        plotmanager$make <- "multigroup"
      }else{
        plotmanager$make <- "summary"
      }
    }
  })
})


shiny::observeEvent(input$enrich_tab_rows_selected,{
  curr_row <- input$enrich_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_pw <- rownames(mSet$analSet$enrich$mummi.resmat)[curr_row]
  pw_i <- which(mSet$analSet$enrich$path.nms == curr_pw)
  cpds = mSet$analSet$enrich$path.hits[[pw_i]]
  hit_tbl = data.table::as.data.table(mSet$analSet$enrich$matches.res)
  myHits <- hit_tbl[Matched.Compound %in% cpds]
  myHits$Mass.Diff <- as.numeric(myHits$Mass.Diff)/(as.numeric(myHits$Query.Mass)*1e-6)
  colnames(myHits) <- c("rn", "identifier", "adduct", "dppm")
  
  enrich$current <- myHits
  
})

# triggers on clicking a row in the match results table
shiny::observeEvent(input$match_tab_rows_selected,{
  curr_row <- input$match_tab_rows_selected # get current row
  if (is.null(curr_row)) return()
  # - - - - - - - - - - - - - - - - - - - - - -
  my_selection$name <- shown_matches$forward_unique[curr_row,'name'][[1]] # get current structure
 
  my_selection$form <- unlist(shown_matches$forward_unique[curr_row,'baseformula']) # get current formula
  my_selection$struct <- unlist(shown_matches$forward_unique[curr_row,'structure']) # get current formula

  toClipboard = switch(input$matchAutocopy,
                       SMILES = my_selection$struct,
                       formula = my_selection$form,
                       name = my_selection$name)
  if(input$autocopy){
    #clipr::write_clip(toClipboard, allow_non_interactive = TRUE)
    shinyjs::runjs(stringr::str_wrap(gsubfn::fn$paste("
  var textArea = document.createElement('textarea');
  textArea.style.position = 'fixed';
  textArea.style.top = 0;
  textArea.style.left = 0;

  textArea.style.width = '2em';
  textArea.style.height = '2em';

  textArea.style.padding = 0;

  textArea.style.border = 'none';
  textArea.style.outline = 'none';
  textArea.style.boxShadow = 'none';

  textArea.style.background = 'transparent';

  textArea.value = '$toClipboard';

  document.body.appendChild(textArea);
  textArea.focus();
  textArea.select();

  try {
    var successful = document.execCommand('copy');
    var msg = successful ? 'successful' : 'unsuccessful';
    console.log('Copying text command was ' + msg);
  } catch (err) {
    console.log('Oops, unable to copy');
  }

  document.body.removeChild(textArea);"), width=10000))
    shiny::updateTextInput(session,
                           "wordcloud_searchTerm",
                           value = as.character(toClipboard))
  }
  
  })


shiny::observeEvent(input$browse_tab_rows_selected,{
  curr_row <- input$browse_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_def <- browse_content$table[curr_row, description]
  output$desc_ui <- shiny::renderText(curr_def)
  my_selection$struct <- browse_content$table[curr_row,c('structure')][[1]]
  my_selection$form <- unlist(browse_content$table[curr_row,'formula']) # get current formula
})

# triggers on clicking a row in the reverse hit results table
shiny::observeEvent(input$hits_tab_rows_selected,{
  curr_row <- input$hits_tab_rows_selected # get current row
  if (is.null(curr_row)) return()
  my_selection$mz <- shown_matches$reverse[curr_row,'query_mz'][[1]]
})