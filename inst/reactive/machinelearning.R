shiny::observeEvent(input$ml_top_x, {
  if(!is.null(mSet)){
    plotmanager$make <- "ml"
  }
},ignoreInit = T, ignoreNULL = T)

shiny::observeEvent(input$show_which_ml,{
  if(!is.null(mSet)){
    split.name = strsplit(input$show_which_ml, split = " - ")[[1]]
    mSet$analSet$ml$last$method <<- split.name[[1]]
    mSet$analSet$ml$last$name <<- split.name[[2]]
    tablemanager$make <- "ml"
    plotmanager$make <- "ml"
  }
},ignoreNULL = T, ignoreInit = T)

shiny::observe({
  if(!is.null(input$ml_method)){
    mdl = caret::getModelInfo()[[input$ml_method]]
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
})

output$ml_name <- shiny::renderUI({
  ml_name = trimws(paste0(if(input$ml_run_on_norm) "norm" else "orig",
                   " ",
                   input$ml_train_perc, "%train",
                   " ",
                   paste0("crossVal-",input$ml_perf_metr,input$ml_folds),
                   " ",
                   if(length(input$ml_batch_covars) > 0) paste0("batchSplit-", paste0(input$ml_batch_covars, collapse="AND")),
                   if(length(input$ml_batch_covars) > 0 & input$ml_batch_sampling != "none"){
                     paste0("balancedMethod-",
                            input$ml_batch_sampling, 
                            "-by-", paste0(input$ml_batch_covars,
                                           collapse = "+")," ")
                   }else{""},
                   if(length(input$ml_include_covars) > 0){
                     paste0("metadataInclude", paste0(input$ml_include_covars, 
                                                      collapse="+"), " ")
                   }else{ "" },
                   #if(input$ml_sampling != "none") paste0("balancedClasses-", input$ml_sampling, ""," ") else "",
                   if(!(input$ml_samp_distr %in% c(" ", "no"))){
                     paste0("usePrevTrainTest-", input$ml_samp_distr, " ")
                   }else{""},
                   if(!is.null(lcl$vectors$ml_train)){paste0("trainOn-", paste0(lcl$vectors$ml_train, collapse="="))} else "",
                   if(!is.null(lcl$vectors$ml_test)){paste0("testOn-", paste0(lcl$vectors$ml_test, collapse="="))} else ""
  ))
  
  caret.mdls <- caret::getModelInfo()
  caret.methods <- names(caret.mdls)
  tune.opts <- lapply(caret.methods, function(mdl) caret.mdls[[mdl]]$parameters)
  names(tune.opts) <- caret.methods
  
  meth.info <- caret.mdls[[input$ml_method]]
  params = meth.info$parameters
  
  tuneGrid = expand.grid(
    {
      lst = lapply(1:nrow(params), function(i){
        info = params[i,]
        inp.val = input[[paste0("ml_", info$parameter)]]
        # - - check for ranges - -
        if(grepl(inp.val, pattern=":")){
          split = strsplit(inp.val,split = ":")[[1]]
          inp.val <- seq(as.numeric(split[1]),
                         as.numeric(split[2]),
                         as.numeric(split[3]))
        }else if(grepl(inp.val, pattern = ",")){
          split = strsplit(inp.val,split = ",")[[1]]
          inp.val <- split
        }
        # - - - - - - - - - - - - -
        switch(as.character(info$class),
               numeric = as.numeric(inp.val),
               character = as.character(inp.val))
      })
      names(lst) = params$parameter
      #lst <- lst[sapply(lst,function(x)all(!is.na(x)))]
      lst
    })
  ml_settings = paste0(unlist(sapply(colnames(tuneGrid), function(x) if(!is.na(tuneGrid[[x]])) paste0(x,"=",tuneGrid[[x]]) else NULL)), collapse="&")
  if(ml_settings != ""){
    ml_name = paste0(ml_name, " ", ml_settings)
  }
  
  ml_name = trimws(ml_name)
  shiny::textInput("ml_name", 
                   label=shiny::h3("Name:"), 
                   value = ml_name)
})

shiny::observeEvent(input$ml_train_ss, {
  keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
  subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
  lcl$vectors$ml_train <<- c(input$subset_var,
                             input$subset_group)
  output$ml_train_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$reset_ml_train, {
  subset.name <- "all"
  lcl$vectors$ml_train <<- NULL
  output$ml_train_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$ml_test_ss, {
  keep.samples <- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
  subset.name <- paste(input$subset_var, input$subset_group, sep = "-")
  lcl$vectors$ml_test <<- c(input$subset_var, input$subset_group)
  output$ml_test_ss <- shiny::renderText(subset.name)
})

shiny::observeEvent(input$reset_ml_test, {
  subset.name <- "all"
  lcl$vectors$ml_test <<- NULL
  output$ml_test_ss <- shiny::renderText(subset.name)
})