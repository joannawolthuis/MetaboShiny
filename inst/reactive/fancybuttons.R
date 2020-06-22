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

update_heatmap = shiny::reactive({
  input$heatmap_topn
}) %>% shiny::debounce(1000)

shiny::observeEvent(update_heatmap(), {
  if(!is.null(mSet$analSet$heatmap)){
    plotmanager$make <- "heatmap" # just reload
  }
})

update_ml = shiny::reactive({
  input$ml_top_x
}) %>% shiny::debounce(1000)

shiny::observeEvent(update_ml(), {
  if(!is.null(mSet$analSet$ml)){
    plotmanager$make <- "ml" # just reload
  }
}) #%>% shiny::debounce(1000)

 
# shiny::observeEvent(input$db_only, {
#   if(input$db_only){
#     MetaboShiny::setOption(lcl$paths$opt.loc, "mode", "dbonly")
#     logged$status <- "dbonly"
#   }else{
#     MetaboShiny::setOption(lcl$paths$opt.loc, "mode", "complete")
#     logged$status <- "logged"
#   }
# },ignoreInit = F, ignoreNULL = T)




