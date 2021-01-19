# reload plots (pca/plsda) if the 2d/3d button is triggered
scatters = c("pca", "plsda", "tsne", "umap", "ica")

lapply(scatters, function(anal){
  refreshers = c("2d3d", "ellipse", "x", "y", "z")
  lapply(refreshers, function(r){
    input_id = paste0(anal,"_",r)
    shiny::observeEvent(input[[input_id]], {
      plotmanager$make <- anal
    }, 
    ignoreNULL = T, 
    ignoreInit = T)    
  })
})

shiny::observeEvent(input$corr_topn, {
  plotmanager$make <- "corr"
}, ignoreNULL = T)

update_heatmap = shiny::reactive({
  input$heatmap_topn
}) %>% shiny::debounce(1000)

shiny::observeEvent(update_heatmap(), {
  if(!is.null(mSet$analSet$heatmap)){
    plotmanager$make <- "heatmap" # just reload
  }
})

shiny::observeEvent(input$plot_ml_mistake, {
  if(!is.null(mSet$analSet$ml)){
    plotmanager$make <- "ml_mistake" # just reload  
  }
},ignoreInit = T, ignoreNULL = T)

update_ml = shiny::reactive({
  input$ml_top_x
  input$ml_curve_type
}) %>% shiny::debounce(1000)

shiny::observeEvent(update_ml(), {
  if(!is.null(mSet$analSet$ml)){
    plotmanager$make <- tablemanager$make <- "ml" # just reload
  }
})




