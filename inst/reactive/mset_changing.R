observeEvent(input$change_cls, {
  if(input$stats_type %in% c("t", "t1f")){
    shiny::updateCheckboxInput(session, "paired", value = T) # auto set for time series 
  }
  mSetter$do <- "change"
  uimanager$refresh <- "general"
})

shiny::observeEvent(input$change_subset, {
  mSetter$do <- "subset"
  uimanager$refresh <- "general"
  })

shiny::observeEvent(input$reset_subset, {
  mSetter$do <- "unsubset"
  uimanager$refresh <- "general"
})

shiny::observeEvent(input$load_storage, {
  mSetter$do <- "load"
  uimanager$refresh <- "general"
})


observeEvent(input$show_prematched_mz_only, {
  if(mSet$metshiParams$prematched){
    mSetter$do <- "refresh"
    uimanager$refresh <- "general"  
  }
},ignoreInit = T)