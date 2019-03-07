# generate positive and negative adduct picker tabs (for csv creation)
# defaults are in the huge global object :-)
observe({
  modes = c("pos", "neg")
  lapply(modes, function(mode){
    output[[paste0(mode, "_add_tab")]] <- DT::renderDataTable({
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', 
                                     selected = global$vectors[[paste0(mode, "_selected_add")]], target="row"),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
})

# toggles when 'select all adducts' is pressed (filled circle)
observeEvent(input$sel_all_adducts, {
  global$vectors$neg_selected_adducts <<- c(1:nrow(global$vectors$pos_adducts))
  global$vectors$pos_selected_adducts <<- c(1:nrow(global$vectors$neg_adducts))
})

# triggers when 'select no adducts' is selected
observeEvent(input$sel_no_adducts, {
  global$vectors$neg_selected_adducts <<- c(0)
  global$vectors$pos_selected_adducts <<- c(0)
})

# triggers when common adducts are to be selected
observeEvent(input$sel_comm_adducts, {
  global$vectors$neg_selected_adducts <<- c(1:3, nrow(global$vectors$pos_adducts))
  global$vectors$pos_selected_adducts <<- c(1, 2, 14:15, nrow(global$vectors$neg_adducts))
})

# adduct table editing from settings tab

values = reactiveValues()

# TODO: fix and re-docment this
observeEvent(input$import_adducts, {
  DF = fread(input$add_tab$datapath)
  output$adduct_tab <- rhandsontable::renderRHandsontable({
    if (!is.null(DF))
      rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
  })
  output$adduct_upload_check <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath('www/yes.png')
    # Return a list containing the filename and alt text
    list(src = filename, width = 20,
         height = 20)
  }, deleteFile = FALSE)
})

# uses rhandsontable for live table editing...
# TODO: fix
adduct_tab_data <- reactive({
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

observe({
  DF = adduct_tab_data()
  shinyFileSave(input, "save_adducts", roots = c(home = '~'), session=session)
  fileinfo <- parseSavePath(roots = c(home = '~'), input$save_adducts)
  if (nrow(fileinfo) > 0) {
    switch(fileinfo$type,
           csv = fwrite(file = fileinfo$datapath, x = DF)
    )}
})