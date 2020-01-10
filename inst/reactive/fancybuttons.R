# reload plots (pca/plsda) if the 2d/3d button is triggered
shiny::observeEvent(input$pca_2d3d, {
  datamanager$reload <- "pca"
}, ignoreNULL = T)

shiny::observeEvent(input$plsda_2d3d, {
  datamanager$reload <- "plsda"
}, ignoreNULL = T)

shiny::observeEvent(input$tsne_2d3d, {
  datamanager$reload <- "tsne"
}, ignoreNULL = T)

shiny::observeEvent(input$tt_nonpar, {
  if(!is.null(input$permz)){
    if(input$permz == "tt"){
      statsmanager$calculate <- "tt"
      datamanager$reload <- "tt" 
    }  
  }
},ignoreInit = TRUE)

shiny::observeEvent(input$tt_eqvar, {
  statsmanager$calculate <- "tt"
  datamanager$reload <- "tt"
},ignoreInit = TRUE)

shiny::observeEvent(input$heatsign, {
  statsmanager$calculate <- "heatmap"
  datamanager$reload <- "heatmap"
})#,ignoreInit = F, ignoreNULL = T)

shiny::observeEvent(input$heatmap_topn, {
  if(!is.null(mSet$analSet$heatmap)){
    datamanager$reload <- "heatmap" # just reload
  }
})#,ignoreInit = F, ignoreNULL = T)

shiny::observeEvent(input$wc_topn_pm, {
  datamanager$reload <- "match_wordcloud_pm"
})#,ignoreInit = TRUE, ignoreNULL = T)


shiny::observeEvent(input$wc_topn, {
  datamanager$reload <- "match_wordcloud"
})#,ignoreInit = TRUE, ignoreNULL = T)

shiny::observeEvent(input$db_only, {
  if(input$db_only){
    MetaboShiny::setOption(lcl$paths$opt.loc, "mode", "dbonly")
    logged$status <- "dbonly"
  }else{
    MetaboShiny::setOption(lcl$paths$opt.loc, "mode", "complete")
    logged$status <- "logged"
  }
},ignoreInit = F, ignoreNULL = T)



