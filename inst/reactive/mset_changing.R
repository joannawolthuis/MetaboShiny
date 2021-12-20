observeEvent(input$change_cls, {
  if(input$stats_type %in% c("t", "t1f")){
    shiny::updateCheckboxInput(session, "paired", value = T) # auto set for time series 
  }
  mSetter$do <- "change"
  #uimanager$refresh <- "general"
})

shiny::observeEvent(input$change_subset, {
  mSetter$do <- "subset"
  #uimanager$refresh <- "general"
  })

shiny::observeEvent(input$reset_subset, {
  mSetter$do <- "unsubset"
  #uimanager$refresh <- "general"
})

shiny::observeEvent(input$load_storage, {
  mSetter$do <- "load"
  #uimanager$refresh <- "general"
})

observeEvent(input$change_subset_mz, {
  if(mSet$metshiParams$prematched){
    mSetter$do <- "subset_mz"
   #uimanager$refresh <- "general"  
  }else{
    shiny::showNotification("No can do, please prematch first!")
  }
},ignoreInit = T)

observeEvent(input$remove_storage, {
  if(!is.null(input$storage_choice)){
    exp = input$storage_choice
    shinyWidgets::confirmSweetAlert(
      session = session,
      inputId = "remove_confirm",
      text = tags$div(
        tags$b("Click upper right ", icon("times"), " button to cancel."),br(),
        shiny::img(#class = "rotategem", 
          src = "gemmy_rainbow.png", 
          width = "70px", height = "70px"),
        br()
      ),
      btn_labels = c("No", "Yes"),
      title = gsubfn::fn$paste("Are you sure you want to remove '$exp'?"),
      #showCloseButton = T,
      html = TRUE
    )
  }
})

observeEvent(input$remove_confirm,{
  if(isTRUE(input$remove_confirm)){
    mSet$storage[[input$storage_choice]] <- NULL
    mSet <<- mSet
    uimanager$refresh <- "general"
  }
},ignoreNULL = T)
