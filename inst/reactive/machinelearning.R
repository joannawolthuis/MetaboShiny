shiny::observeEvent(input$ml_top_x, {
  if(!is.null(mSet)){
    plotmanager$make <- "ml"
  }
},ignoreInit = T, ignoreNULL = T)

shiny::observeEvent(input$show_which_ml,{
  if(!is.null(mSet) & input$show_which_ml != ""){
    split.name = strsplit(input$show_which_ml, split = " - ")[[1]]
    mSet$analSet$ml$last$method <<- split.name[[1]]
    mSet$analSet$ml$last$name <<- split.name[[2]]
    tablemanager$make <- c("vennrich", "ml")
    plotmanager$make <- "ml"
  }
},ignoreNULL = T, ignoreInit = T)

shiny::observeEvent(input$ml_plot_facet, {
  if(input$ml_plot_facet != "don't facet"){
    plotmanager$make <- "ml"
  }
},ignoreNULL = T, ignoreInit = T)


shiny::observeEvent(input$ml_batch_size_sampling, {
  if(input$ml_batch_size_sampling){
    shinyWidgets::updateRadioGroupButtons(session, 
                                          "ml_sampling", 
                                          choices = c(`don't` = "none", 
                                                      `<i class='fa fa-arrow-down'></i> downsample` = "down",
                                                      `ROSE` = "rose", 
                                                      #`SMOTE` = "smote",
                                                      `upsample <i class='fa fa-arrow-up'></i>` = "up"))
  }else{
    shinyWidgets::updateRadioGroupButtons(session, 
                                          "ml_sampling", 
                                          choices = c(`don't` = "none",
                                                      `<i class='fa fa-arrow-down'></i> downsample` = "down",
                                                      `ROSE` = "rose",
                                                      `SMOTE` = "smote",
                                                      `ADASYN` = "adasyn",
                                                      `upsample <i class='fa fa-arrow-up'></i>` = "up"))
    
  }
})

shiny::observe({
  if(!is.null(input$ml_method)){
    if(!input$ml_method %in% c("", " ")){
      sel_mdl = input$ml_method
      if(sel_mdl == "glm (logistic)") sel_mdl <- "glm"
      
      mdl = caret::getModelInfo()[[sel_mdl]]
      params <- mdl$parameters
      output$ml_params <- renderUI({
        list(
          shiny::helpText(mdl$label),
          shiny::hr(),
          h2("Tuning settings"),
          lapply(1:nrow(params), function(i){
            row = params[i,]
            list(
              shiny::textInput(inputId = paste0("ml_", row$parameter),
                               label = row$parameter,
                               value=if(input$ml_method=="glmnet"){
                                 switch(row$parameter,
                                        alpha = 1,
                                        lambda = "0:1:0.01")
                               }),
              shiny::helpText(paste0(row$label, " (", row$class, ")."))
            )
          })
        )
      })
    }
  }
})

shiny::observeEvent(input$ml_train_ss, {
  keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
  if(length(lcl$lists$ml_train) > 0){
    lcl$lists$ml_train[[length(lcl$lists$ml_train) + 1]] <<- list(input$subset_var,
                                                                  input$subset_group)
    subset.name <- paste0(sapply(lcl$lists$ml_train, function(ss_pair){
      paste(ss_pair[1], ss_pair[2], sep = "-")
    }),collapse = "\n")

  }else{
    subset.name <- paste(input$subset_var, 
                         input$subset_group, sep = "-")
    lcl$lists$ml_train <<- list(list(input$subset_var,
                                     input$subset_group))
  }
  output$ml_train_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$reset_ml_train, {
  subset.name <- "all"
  lcl$lists$ml_train <<- list()
  output$ml_train_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$ml_test_ss, {
  keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
  if(length(lcl$lists$ml_test) > 0){
    lcl$lists$ml_test[[length(lcl$lists$ml_test) + 1]] <<- list(input$subset_var,
                                                               input$subset_group)
    subset.name <- paste0(sapply(lcl$lists$ml_test, function(ss_pair){
      paste(ss_pair[1], ss_pair[2], sep = "-")
      }),collapse = "\n")
    
  }else{
    subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
    lcl$lists$ml_test <<- list(list(input$subset_var,
                                     input$subset_group))
  }
  output$ml_test_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$reset_ml_test, {
  subset.name <- "all"
  lcl$vectors$ml_test <<- NULL
  output$ml_test_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$queue_ml, {
  imp = shiny::isolate(shiny::reactiveValuesToList(input))
  ml_args = imp[grep(names(imp),pattern = "^ml_")]
  ml_args = ml_args[grep("clicked|current|rows|ss|tab_|mistake|curve|results", names(ml_args),invert = T)]
  ml_args$ml_train_subset <- lcl$lists$ml_train
  ml_args$ml_test_subset <- lcl$lists$ml_test
  # save to queue
  ml_queue$jobs[[ml_args$ml_name]] <- ml_args
  # new random name for next job? TODO: only do if user not specifying their own
  # shiny::updateTextInput(session, "ml_name", value = stringi::stri_rand_strings(1, 
  #                                                                               10, 
  #                                                                               pattern = "[A-Za-z0-9]"))
})

shiny::observeEvent(input$queue_ml_del, {
  rows = input$ml_queue_all_rows_selected
  ml_queue$jobs = ml_queue$jobs[-rows]
})

shiny::observeEvent(input$clear_ml_runs, {
  shinyWidgets::confirmSweetAlert(
    session = session,
    inputId = "clear_ml_sure",
    text = tags$div(
      tags$b("Click upper right ", icon("times"), " button to cancel."),br(),
      br(),
      shiny::img(class = "imagetop", 
                 src = "www/metshi_heart_bezel.png", 
                 height = "70px"),
      br()
    ),
    btn_labels = c("No", "Yes"),
    title = "Erase all machine learning results?",
    #showCloseButton = T,
    html = TRUE
  )
})

observeEvent(input$clear_ml_sure,{
  if(isTRUE(input$clear_ml_sure)){
    mSet$analSet$ml <<- NULL
    uimanager$refresh <- "ml"
  }
},ignoreNULL = T)

shiny::observeEvent(input$ml_queue_all_rows_selected, {
  params = ml_queue$jobs[[input$ml_queue_all_rows_selected]]
  param_dt = data.table::data.table(parameter = names(params),
                         value = params)
  output$ml_queue_sel <- DT::renderDataTable({
    MetaboShiny::metshiTable(content = param_dt)
  }, server = T)
  })