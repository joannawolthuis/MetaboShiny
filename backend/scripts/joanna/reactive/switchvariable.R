# triggers when the 'change variable' dropdown menu is filled and button is clicked
observeEvent(input$change_cls, {

  # check if previous analysis storage already exists, if not, make it
  if(!("storage" %in% names(mSet))){
    mSet$storage <<- list()
  }

  mset_name = mSet$dataSet$cls.name

  # save previous analyses (should be usable in venn diagram later)
  mSet$storage[[mset_name]] <<- list(analysis = mSet$analSet)

  local$constants$last_mset <<- mset_name

  # change current variable of interest to user pick from covars table
  mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,input$stats_var, with=F][[1]])

  # adjust bivariate/multivariate (2, >2)...
  mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))

  # adjust name of experimental variable
  if(grepl(mSet$dataSet$cls.name, pattern = ":")){
    subset_name <- paste0(":", gsub(mSet$dataSet$cls.name, pattern = ".*:", replacement = ""))
  }else{
    subset_name <- ""
  }

  new_name <- paste0(input$stats_var, subset_name)
  mSet$dataSet$cls.name <<- new_name

  if(new_name %in% names(mSet$storage)){
    mSet$analSet <<- mSet$storage[[input$stats_var]]$analysis
  }else{
    # remove old analSet
    mSet$analSet <<- NULL
  }

  datamanager$reload <- "general"

  output$curr_name <- renderText({mSet$dataSet$cls.name})

  updateNavbarPage(session, "statistics", selected = "inf")


})
