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

topn_sliders <- c("corr", "fc", "tt", "aov", "heatmap")

lapply(topn_sliders, function(anal){
  r <- shiny::reactive({
    input[[paste0(anal,"_topn")]]
  }) %>% shiny::debounce(2000)
  shiny::observeEvent(r(), { # so it doesn't constantly update
    if(anal == "heatmap"){
      statsmanager$calculate <- anal
    }
    plotmanager$make <- anal
  }, ignoreNULL = T)
})

lapply(c("x", "y", "posclass"), function(axis){
  shiny::observeEvent(input[[paste0("ml_plot_", axis)]], {
    if(!is.null(mSet$analSet$ml)){
      plotmanager$make <- "ml" # just reload  
    }
  })  
})

shiny::observeEvent(input$reload_ml_stats, {
  if(!is.null(mSet$analSet$ml)){
    plotmanager$make <- tablemanager$make <- "ml" # just reload
  }
},ignoreInit = T, ignoreNULL = T)

update_ml = shiny::reactive({
  input$ml_topn
  input$ml_curve_type
}) %>% shiny::debounce(1000)

shiny::observeEvent(update_ml(), {
  if(!is.null(mSet$analSet$ml)){
    plotmanager$make <- tablemanager$make <- "ml" # just reload
  }
})




