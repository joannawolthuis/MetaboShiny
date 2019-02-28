# reload plots (pca/plsda) if the 2d/3d button is triggered
observeEvent(input$pca_2d3d, {
  datamanager$reload <- "pca"
}, ignoreNULL = T)

observeEvent(input$plsda_2d3d, {
  datamanager$reload <- "plsda"
}, ignoreNULL = T)

observeEvent(input$tt_nonpar, {
  statsmanager$calculate <- "tt"
  datamanager$reload <- "tt"
},ignoreInit = TRUE)

observeEvent(input$tt_eqvar, {
  statsmanager$calculate <- "tt"
  datamanager$reload <- "tt"
},ignoreInit = TRUE)

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
  print("lol 1")
  statsmanager$calculate <- "heatmap"
  datamanager$reload <- "heatmap"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$heatsign, {
  print("lol 2")
  statsmanager$calculate <- "heatmap"
  datamanager$reload <- "heatmap"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$heatmap_topn, {
  print("lol 3")
  if(!is.null(mSet$analSet$heatmap)){
    datamanager$reload <- "heatmap" # just reload
  }
},ignoreInit = TRUE, ignoreNULL = T)

