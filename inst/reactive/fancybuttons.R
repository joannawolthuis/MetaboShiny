# reload plots (pca/plsda) if the 2d/3d button is triggered
shiny::observeEvent(input$pca_2d3d, {
  plotmanager$make <- "pca"
}, ignoreNULL = T)

shiny::observeEvent(input$plsda_2d3d, {
  plotmanager$make <- "plsda"
}, ignoreNULL = T)

shiny::observeEvent(input$tsne_2d3d, {
  plotmanager$make <- "tsne"
}, ignoreNULL = T)

# shiny::observeEvent(input$tt_nonpar, {
#   if(!is.null(input$permz)){
#     if(input$permz == "tt"){
#       statsmanager$calculate <- "tt"
#       plotmanager$make <- "tt" 
#     }  
#   }
# },ignoreInit = TRUE)
# 
# shiny::observeEvent(input$tt_eqvar, {
#   statsmanager$calculate <- "tt"
#   plotmanager$make <- "tt"
# },ignoreInit = TRUE)
# 
# shiny::observeEvent(input$heatsign, {
#   statsmanager$calculate <- "heatmap"
#   plotmanager$make <- "heatmap"
# })#,ignoreInit = F, ignoreNULL = T)
# 
shiny::observeEvent(input$heatmap_topn, {
  if(!is.null(mSet$analSet$heatmap)){
    plotmanager$make <- "heatmap" # just reload
  }
})#,ignoreInit = F, ignoreNULL = T)

shiny::observeEvent(input$ml_top_x, {
  if(!is.null(mSet$analSet$ml)){
    plotmanager$make <- "ml" # just reload
  }
})#,ignoreInit = F, ignoreNULL = T)

 
# shiny::observeEvent(input$db_only, {
#   if(input$db_only){
#     MetaboShiny::setOption(lcl$paths$opt.loc, "mode", "dbonly")
#     logged$status <- "dbonly"
#   }else{
#     MetaboShiny::setOption(lcl$paths$opt.loc, "mode", "complete")
#     logged$status <- "logged"
#   }
# },ignoreInit = F, ignoreNULL = T)




