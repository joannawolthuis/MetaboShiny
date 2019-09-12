output$match_tab <- DT::renderDataTable({
  # don't show some columns but keep them in the original table, so they can be used
  # for showing molecule descriptions, structure
  DT::datatable(shown_matches$forward$unique,
                selection = 'single',
                autoHideNavigation = T,
                #filter = "top",
                options = list(lengthMenu = c(5, 10, 15),
                               pageLength = 5,
                               #searchCols = lcl$default_search_columns,
                               #search = list(regex = FALSE, caseInsensitive = FALSE, search = lcl$default_search),
                               columnDefs = list(list(visible=FALSE, 
                                                      targets=c(which(colnames(shown_matches$forward$unique) %in% gbl$vectors$hide_match_cols),
                                                                which(colnames(shown_matches$forward$unique)=="isocat"))))
                )
  )
})

output$hits_tab <- DT::renderDataTable({
  DT::datatable(shown_matches$reverse,
                selection = 'single',
                autoHideNavigation = T,
                options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
})

output$browse_tab <-DT::renderDataTable({
  DT::datatable(browse_content$table,
                selection = 'single',
                autoHideNavigation = T,
                options = list(lengthMenu = c(5, 10, 15),
                               pageLength = 15,
                               columnDefs = list(list(visible=FALSE, 
                                                      targets=which(colnames(browse_content$table) %in% c("description", "structure", "formula", "charge"))))))
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

output$mzgroup_add_tab <- DT::renderDataTable({
  DT::datatable(data.table::data.table(adduct=adducts$Name),
                selection = list(mode = 'multiple',
                                 target="row"),
                options = list(pageLength = 5, dom = 'tp'),
                rownames = F)
})

# toggles when 'select all adducts' is pressed (filled circle)
shiny::observeEvent(input$sel_all_adducts, {
  lcl$vectors$neg_selected_adducts <<- c(1:nrow(lcl$vectors$pos_adducts))
  lcl$vectors$pos_selected_adducts <<- c(1:nrow(lcl$vectors$neg_adducts))
})

# triggers when 'select no adducts' is selected
shiny::observeEvent(input$sel_no_adducts, {
  lcl$vectors$neg_selected_adducts <<- c(0)
  lcl$vectors$pos_selected_adducts <<- c(0)
})

# triggers when common adducts are to be selected
shiny::observeEvent(input$sel_comm_adducts, {
  lcl$vectors$neg_selected_adducts <<- c(1:3, nrow(lcl$vectors$pos_adducts))
  lcl$vectors$pos_selected_adducts <<- c(1, 2, 14:15, nrow(lcl$vectors$neg_adducts))
})

# adduct table editing from settings tab

values = shiny::reactiveValues()

# TODO: fix and re-docment this
shiny::observeEvent(input$import_adducts, {
  DF = data.table::fread(input$add_tab$datapath)
  output$adduct_tab <- rhandsontable::renderRHandsontable({
    if (!is.null(DF))
      rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
  })
  output$adduct_upload_check <- shiny::renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath('www/yes.png')
    # Return a list containing the filename and alt text
    list(src = filename, width = 20,
         height = 20)
  }, deleteFile = FALSE)
})

# uses rhandsontable for live table editing...
# TODO: fix
adduct_tab_data <- shiny::reactive({
  if (!is.null(input$adduct_tab)){
    DF = rhandsontable::hot_to_r(input$adduct_tab)
  } else {
    if (is.null(values[["DF"]]))
      DF = adducts
    else
      DF = values[["DF"]]
  }
  values[["DF"]] = DF
  # ---------------
  DF
})

output$adduct_tab <- rhandsontable::renderRHandsontable({
  DF = adduct_tab_data()
  if (!is.null(DF))
    rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
})

shiny::observe({
  DF = adduct_tab_data()
  shinyFiles::shinyFileSave(input, "save_adducts", roots = c(home = '~'), session=session)
  fileinfo <- shinyFiles::parseFilePaths(roots = c(home = '~'), input$save_adducts)
  if (nrow(fileinfo) > 0) {
    switch(fileinfo$type,
           csv = fwrite(file = fileinfo$datapath, x = DF)
    )}
})

output$magicball_add_tab <- DT::renderDataTable({
  if(any(unlist(scanmode))){
    DT::datatable(data.table::data.table(Adduct = if(all(unlist(scanmode))){
      adducts$Name
    }else{adducts[scanmode %in% Ion_mode]$Name}),
    selection = list(mode = 'multiple',
                     selected = lcl$vectors[[paste0(scanmode, "_selected_add")]], target="row"),
    options = list(pageLength = 5, dom = 'tp',
                   columnDefs = list(list(className = 'dt-center', targets = "_all"))),
    rownames = F)  
  }
})
