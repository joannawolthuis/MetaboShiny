shiny::observe({
  if(!(input$ml_specific_mzs %in% c("no", " ", "", "none", "manual")) & !is.null(input$ml_specific_mzs)){
    mzs = getAllHits(mSet,
                     input$ml_specific_mzs)
    sigthresh = sum(mzs$significant)
    shiny::updateSliderInput(session, "ml_mzs_topn", max = nrow(mzs))
    output$ml_specific_mzs_sigcount = shiny::renderText(paste0("Significant hits: ", sigthresh))
  }
})