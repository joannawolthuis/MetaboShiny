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

observeEvent(input$wc_topn, {
  datamanager$reload <- "matchplots"
},ignoreInit = TRUE, ignoreNULL = T)

observeEvent(input$db_only, {
  if(input$db_only){
    print("a")
    setOption(lcl$paths$opt.loc, "mode", "dbonly")
    logged$status <- "dbonly"
  }else{
    print("b")
    setOption(lcl$paths$opt.loc, "mode", "complete")
    logged$status <- "logged"
  }
},ignoreInit = F, ignoreNULL = T)

