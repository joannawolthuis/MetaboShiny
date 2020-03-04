output$match_tab <- DT::renderDataTable({
  # don't show some columns but keep them in the original table, so they can be used
  # for showing molecule descriptions, structure
  empty = data.table::data.table(" " = rep("", nrow(shown_matches$forward_unique)))
  
  content = cbind(shown_matches$forward_unique, 
                  empty)
  
  targets = unique(c(which(colnames(content) %in% gbl$vectors$hide_match_cols),
                     which(colnames(content) == "isocat")))
  
  MetaboShiny::metshiTable(
              content = content,
              options = list(
                fixedColumns = list(
                rightColumns = 9,
                heightMatch = 'none'),
                columnDefs = list(list(visible = FALSE,
                                       targets = targets - 1))),
              rownames = F)
})

output$enrich_tab <-DT::renderDataTable({
  MetaboShiny::metshiTable(content = enrich$overview)
})

output$enrich_pw_tab <-DT::renderDataTable({
  MetaboShiny::metshiTable(content = enrich$current)
})

output$hits_tab <- DT::renderDataTable({
  MetaboShiny::metshiTable(content = shown_matches$reverse)
})

output$browse_tab <-DT::renderDataTable({
  MetaboShiny::metshiTable(content = browse_content$table,
              options = list(columnDefs = list(list(visible=FALSE, 
                                                    targets=which(colnames(browse_content$table) %in% c("description", "structure", "formula", "charge"))))
              ))
}, server=T)

# generate positive and negative adduct picker tabs (for csv creation)
# defaults are in the huge global object :-)
shiny::observe({
  modes = c("pos", "neg")
  lapply(modes, function(mode){
    output[[paste0(mode, "_add_tab")]] <- DT::renderDataTable({
      DT::datatable(gbl$vectors$pos_adducts,
                    selection = list(mode = 'multiple',
                                     selected = lcl$vectors[[paste0(mode, "_selected_add")]], target="row"),
                    options = list(pageLength = 5, dom = 'tp'),
                    rownames = F)
    })
  })
})

# adduct table editing from settings tab

output$ml_tab <- DT::renderDataTable({
  MetaboShiny::metshiTable(content = data.table::data.table("nothing selected"="Please select a model from ROC plot or left-hand table!"))
})

observeEvent(input$ml_overview_tab_rows_selected, {
  attempt = unique(mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc$perf)[input$ml_overview_tab_rows_selected,]$attempt
  xvals <- mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]$roc
  output$ml_tab <- DT::renderDataTable({
    imp <- data.table::as.data.table(xvals$imp[[attempt]], keep.rownames = T)
    colnames(imp) <- c("mz", "importance")
    imp <- imp[importance > 0,]
    lcl$tables$ml_roc <<- data.frame(importance = imp$importance,
                                     row.names = gsub(imp$mz,
                                                      pattern = "`|`",
                                                      replacement=""))
    MetaboShiny::metshiTable(content = lcl$tables$ml_roc)
  })
})

values = shiny::reactiveValues()

lapply(c("adducts", "adduct_rules"), function(prefix){
  switch(prefix,
         adducts = {
           show.tab = adducts
           tab.pfx = "add_tab"
         },
         adduct_rules = {
           show.tab = adduct_rules
           tab.pfx = "add_rule_tab"
         })
  
  shiny::observeEvent(input[[paste0("import_", prefix)]], {
    file.copy(input[[tab.pfx]]$datapath, 
              file.path(lcl$paths$work_dir, 
                        paste0(prefix, ".csv")), 
              overwrite = T)
    switch(prefix,
           adducts = {
             adducts <<- data.table::fread(input$add_tab$datapath)
             show.tab = adducts
           },
           adduct_rules = {
             adduct_rules <<- data.table::fread(input$add_rule_tab$datapath)
             show.tab = adduct_rules
           })
    output[[tab.pfx]] <- rhandsontable::renderRHandsontable({
      if (!exists(prefix))
        rhandsontable::rhandsontable(show.tab, stretchH = "all", useTypes = TRUE)
    })
  }) 
  
  shiny::observeEvent(input[[paste0("import_", prefix)]], {
    file.copy(input[[tab.pfx]]$datapath, 
              file.path(lcl$paths$work_dir, 
                        paste0(prefix,".csv")), 
              overwrite = T)
    
    switch(prefix,
           adducts = {
             adducts <<- data.table::fread(input$add_tab$datapath)
             show.tab = adducts
           },
           adduct_rules = {
             adduct_rules <<- data.table::fread(input$add_rule_tab$datapath)
             show.tab = adduct_rules
           })
    
    output[[paste0(prefix, "_tab")]] <- rhandsontable::renderRHandsontable({
      if (!is.null(show.tbl))
        rhandsontable::rhandsontable(show.tbl, stretchH = "all", useTypes = TRUE)
    })
  })
  
})

# uses rhandsontable for live table editing...
# TODO: fix
adduct_tab_data <- shiny::reactive({
  if (!is.null(input$adduct_tab)){
    DF = rhandsontable::hot_to_r(input$adduct_tab)
  } else {
    if (is.null(values[["adducts"]]))
      DF = adducts
    else
      DF = values[["adducts"]]
  }
  values[["adducts"]] = DF
  # ---------------
  DF
})

adduct_rules_tab_data <- shiny::reactive({
  if (!is.null(input$adduct_rules_tab)){
    DF = rhandsontable::hot_to_r(input$adduct_rules_tab)
  } else {
    if (is.null(values[["adduct_rules"]]))
      DF = adduct_rules
    else
      DF = values[["adduct_rules"]]
  }
  values[["adduct_rules"]] = DF
  # ---------------
  DF
})

output$adduct_tab <- rhandsontable::renderRHandsontable({
  DF = adduct_tab_data()
  if (!is.null(DF))
    rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
})

output$adduct_rules_tab <- rhandsontable::renderRHandsontable({
  # TODO: implement https://smartsview.zbh.uni-hamburg.de/rest
  DF = adduct_rules_tab_data()
  if (!is.null(DF))
    rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
})

observeEvent(input$save_adducts, {
  shiny::showNotification("Saving changes to adducts ...")
  fwrite(values$adducts, file = file.path(lcl$paths$work_dir, "adducts.csv"))
  fwrite(values$adduct_rules, file = file.path(lcl$paths$work_dir, "adduct_rules.csv"))
  adducts <<- fread(file.path(lcl$paths$work_dir, "adducts.csv"))
  adduct_rules <<- fread(file.path(lcl$paths$work_dir, "adduct_rules.csv"))
})

# output$magicball_add_tab <- DT::renderDataTable({
#   if(any(unlist(scanmode))){
#     DT::datatable(data.table::data.table(Adduct = if(all(unlist(scanmode))){
#       adducts$Name
#     }else{adducts[scanmode %in% Ion_mode]$Name}),
#     selection = list(mode = 'multiple',
#                      selected = lcl$vectors[[paste0(scanmode, "_selected_add")]], target="row"),
#     extensions = 'Scroller',
#     options = list(deferRender = TRUE,
#                    scrollY = 200,
#                    scroller = TRUE, dom = 'tp',
#                    columnDefs = list(list(className = 'dt-center', targets = "_all"))),
#     rownames = F)  
#   }
# })
