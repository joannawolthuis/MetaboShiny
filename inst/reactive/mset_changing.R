observeEvent(input$change_cls, {
  mSetter$do <- "change"
  datamanager$reload <- "general"
})

shiny::observeEvent(input$paired, {
  if(!is.null(mSet)){
    if(input$paired){
      mSetter$do <- "pair"
    }else{
      mSetter$do <- "unsubset"
      mSet$dataSet$paired <<- FALSE
    }
  }
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
