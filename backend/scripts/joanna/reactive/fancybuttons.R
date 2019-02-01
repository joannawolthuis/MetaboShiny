# reload plots (pca/plsda) if the 2d/3d button is triggered
observeEvent(input$pca_2d3d, {
  datamanager$reload <- "pca"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$plsda_2d3d, {
  datamanager$reload <- "plsda"
},ignoreInit = TRUE, ignoreNULL = T)

# tt

output$tt_parbutton <- shiny::renderUI({
  if("tt" %in% names(mSet$analSet)){
    if("V" %in% colnames(mSet$analSet$tt$sig.mat)){
      switchButton("tt_nonpar", "Non-parametric?", col="BW", type="YN", value = T)
    }else{
      switchButton("tt_nonpar", "Non-parametric?", col="BW", type="YN", value = F)
    }
  }else{
    output$tt_parbutton <- shiny::renderUI({
      switchButton("tt_nonpar", "Non-parametric?", col="BW", type="YN", value = F)
    })
  }
})

observeEvent(input$tt_nonpar, {
  statsmanager$calculate <- "tt"
  datamanager$reload <- "tt"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$tt_eqvar, {
  statsmanager$calculate <- "tt"
  datamanager$reload <- "tt"
},ignoreInit = TRUE, ignoreNULL = T)

# set default mode for heatmap top hits pick button (tt/fc or asca/meba)
heatbutton <- reactiveValues(status = "ttfc")

# render heatmap button
output$heatbutton <- renderUI({
  if(is.null(heatbutton$status)){
    NULL
  }else{
    switch(heatbutton$status,
           asmb = switchButton(inputId = "heatmode",
                               label = "Use data from:", 
                               value = TRUE, col = "BW", type = "ASMB"),
           ttfc = switchButton(inputId = "heatmode",
                               label = "Use data from:", 
                               value = TRUE, col = "BW", type = "TTFC")
    )
  }
})

observeEvent(input$heatmode, {
  statsmanager$calculate <- "heatmap"
  datamanager$reload <- "heatmap"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$heatsign, {
  statsmanager$calculate <- "heatmap"
  datamanager$reload <- "heatmap"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$heatmap_topn, {
  if(!is.null(mSet$analSet$heatmap)){
    datamanager$reload <- "heatmap" # just reload
  }
},ignoreInit = TRUE, ignoreNULL = T)

