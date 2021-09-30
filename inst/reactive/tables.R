output$match_tab <- DT::renderDataTable({
  # don't show some columns but keep them in the original table, so they can be used
  # for showing molecule descriptions, structure
  empty = data.table::data.table(" " = rep("", nrow(shown_matches$forward_unique)))
  
  content = cbind(shown_matches$forward_unique, 
                  empty)
  
  targets = unique(c(which(colnames(content) %in% gbl$vectors$hide_match_cols),
                     which(colnames(content) == "isocat")))
  
  colnames(content)[colnames(content) == "compoundname"] <- "name"
  MetaboShiny::metshiTable(
              content = content,
              options = list(
                fixedColumns = list(
                rightColumns = 9,
                heightMatch = 'none'),
                columnDefs = list(list(visible = FALSE,
                                       targets = targets - 1))),
              rownames = F)
}, server = F)

output$enrich_tab <- DT::renderDataTable({
  MetaboShiny::metshiTable(content = enrich$overview)
}, server = F)

output$enrich_pw_tab <-DT::renderDataTable({
  MetaboShiny::metshiTable(content = enrich$current)
}, server = F)

output$hits_tab <- DT::renderDataTable({
  MetaboShiny::metshiTable(content = shown_matches$reverse)
}, server = F)

output$browse_tab <-DT::renderDataTable({
  MetaboShiny::metshiTable(content = browse_content$table,
              options = list(columnDefs = list(list(visible=FALSE, 
                                                    targets=which(colnames(browse_content$table) %in% c("description", "structure", "formula", "charge", "source"))))
              ))
}, server=T)

output$ml_queue_all <- DT::renderDataTable({
  MetaboShiny::metshiTable(content = data.table::data.table(name = names(ml_queue$jobs)))
}, server = T)

# adduct table editing from settings tab

output$ml_tab <- DT::renderDataTable({
  MetaboShiny::metshiTable(content = data.table::data.table("nothing selected" = "Please select a model from ROC plot or left-hand table!"))
}, server = T)

shiny::observeEvent(input$ml_overview_tab_rows_selected, {
  attempt <- lcl$tables$ml_roc_all[input$ml_overview_tab_rows_selected,]$attempt
  data <- mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]
  xvals <- data$roc
  output$ml_tab <- DT::renderDataTable({
    imp <- xvals$imp[[attempt]]
    imp <- data.table::data.table(mz = rownames(imp),
                                  importance = imp[[1]])
    imp <- imp[importance > 0,]
    fixed.mzs = gsub("^X", "", imp$mz)
    fixed.mzs = gsub("\\.$", "-", fixed.mzs)
    lcl$tables$ml_roc <<- data.frame(importance = imp$importance,
                                    row.names = fixed.mzs)
    MetaboShiny::metshiTable(content = lcl$tables$ml_roc)
  }, server = F)
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
    
    # file.copy(input[[tab.pfx]]$datapath, 
    #           file.path(lcl$paths$work_dir, 
    #                     paste0(prefix,".csv")), 
    #           overwrite = T)
    
    switch(prefix,
           adducts = {
             show.tab = adducts
           },
           adduct_rules = {
             show.tab = adduct_rules
           })
    
    output[[paste0(prefix, "_tab")]] <- rhandsontable::renderRHandsontable({
      if (!is.null(show.tab))
        rhandsontable::rhandsontable(show.tab, stretchH = "all", useTypes = TRUE)
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

shiny::observeEvent(input$save_adducts, {
  
  shiny::showNotification("Saving & applying changes to adducts ...")
  
  for(tab.pfx in c("adducts", "adduct_rules")){
    file.copy(input[[tab.pfx]]$datapath,
              file.path(lcl$paths$work_dir,
                        paste0(prefix,".csv")),
              overwrite = T)  
  }
  
  adducts <<- data.table::fread(file.path(lcl$paths$work_dir, "adducts.csv"))
  adduct_rules <<- data.table::fread(file.path(lcl$paths$work_dir, "adduct_rules.csv"))
  selAdd = intersect(input$score_add, 
                     adducts$Name)
  shinyWidgets::updatePickerInput("score_add", 
                                  choices=adducts$Name,
                                  selected = selAdd)
})
