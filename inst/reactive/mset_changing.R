observeEvent(input$change_cls, {
  if(input$stats_type %in% c("t", "t1f")){
    shiny::updateCheckboxInput(session, "paired", value = T) # auto set for time series 
  }
  mSetter$do <- "change"
  datamanager$reload <- "general"
})

shiny::observeEvent(input$change_subset, {
  mSetter$do <- "subset"
  datamanager$reload <- "general"
  })

shiny::observeEvent(input$reset_subset, {
  mSetter$do <- "unsubset"
  datamanager$reload <- "general"
})

shiny::observeEvent(input$load_storage, {
  mSetter$do <- "load"
  datamanager$reload <- "general"
})
