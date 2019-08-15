# set default timemode switcher button mode to hidden
timebutton <- reactiveValues(status = "off")

observe({
  if(!is.null(mSet)){
    print(mSet$timeseries)
    if(is.null(mSet$timeseries)) mSet$timeseries <<- FALSE
    datamanager$reload <- "general"
  }else{
    # hide time series button
    timebutton$status <- "off"
    heatbutton$status <- "asmb"
  }
})

observe({
  if(timebutton$status == "off") updateCheckboxInput(session, "timecourse_trigger", label = NULL, value = FALSE)
})

# render time series swap button
output$timebutton <- renderUI({
  if (is.null(timebutton$status)) {
    input$timecourse_trigger <- FALSE
  } else{
    switch(timebutton$status,
           off = NULL,
           on = switchButton(inputId = "timecourse_trigger", # this is used to trigger time series dataswap
                             label = "Toggle time course mode?",
                             value = mSet$timeseries,
                             col = "BW",
                             type = "YN"))
  }
})


observeEvent(input$timecourse_trigger, {

  if(!("storage" %in% names(mSet))){
    mSet$storage <<- list()
  }

  if(input$timecourse_trigger){
    # change to timecourse mode

    mSet$timeseries <<- TRUE

    # subset
    subset <- if(grepl(mSet$dataSet$cls.name, pattern = ":")){
      strsplit(x = gsub(mSet$dataSet$cls.name, pattern = ".*:", replacement = ""), split = "\\+")[[1]]
    }else{
      c()
    }

    # save previous analyses (should be usable in venn diagram later)
    mset_name = get_mset_name(mainvar = gsub(mSet$dataSet$cls.name, pattern = ":.*|\\(.*", replacement = ""),
                              subsetvar = NULL,
                              subsetgroups = subset,
                              timeseries = T)

    # TODO: use this in venn diagram creation
    mSet$storage[[mset_name]] <<- list(analysis = mSet$analSet)
    mSet$dataSet$cls.name <<- mset_name

    # adjust mset design type (necessary for metaboanalystr)
    SetDesignType(mSet, "time")
    mSet$analSet$type <<- "time"
    # rename some factors of interest (your experimental variable, and 'time') as A and B
    facA <- as.factor(mSet$dataSet$covars[,gsub(mSet$dataSet$cls.name, pattern = ":.*|\\(.*", replacement = ""), with=F][[1]])
    facB <- mSet$dataSet$covars[,"time"][[1]]


    # change mSet experimental factors (these are used in ASCA/MEBA etc.)
    mSet$dataSet$exp.fac <<- as.factor(facA)
    mSet$dataSet$time.fac <<- as.factor(facB)
    mSet$dataSet$facA <<- mSet$dataSet$exp.fac;
    mSet$dataSet$facB <<- mSet$dataSet$time.fac;
    mSet$dataSet$facA.lbl <<- mSet$dataSet$cls.name
    mSet$dataSet$facB.lbl <<- "time"

    # change interface to timeseries mode (make 'interface' manager do it)
    interface$mode <<- "time"

    # change heatmap chooser to asca/meba because those are timeseries-specific
    heatbutton$status <<- "asmb"

    # REMOVE PREVIOUS ANALYSIS TO TRIGGER RELOAD (or the PCA won't reload)
    mSet$analSet$pca <<- NULL
    output$curr_name <- renderText({mSet$dataSet$cls.name})

  }else{

    mSet$timeseries <<- FALSE

    # change back to normal mode
    subset <- if(grepl(mSet$dataSet$cls.name, pattern = ":")){
      strsplit(x = gsub(mSet$dataSet$cls.name, pattern = ".*:", replacement = ""), split = "\\+")[[1]]
    }else{
      c()
    }

    # save previous analyses (should be usable in venn diagram later)
    mset_name = get_mset_name(mainvar = gsub(mSet$dataSet$cls.name, pattern = ":.*|\\(.*", replacement = ""),
                              subsetvar = NULL,
                              subsetgroups = subset,
                              timeseries = F)

    mSet$storage[[mset_name]] <<- list(analysis = mSet$analSet)
    mSet$analSet$type <<- "stat"

    # - - - - - - - - - - - -
    # rename experimental factors
    mSet$dataSet$cls.name <<- mset_name
    mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,gsub(mSet$dataSet$cls.name, pattern = ":.*|\\(.*", replacement = ""), with=F][[1]])
    mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
    # remove old analSet
    mSet$analSet$pca <<- NULL

    # change heatmap button back to t-test/fold-change bivariate mode
    heatbutton$status <- "ttfc"

    # reset interface
    if(mSet$dataSet$cls.num <= 1){
      interface$mode <- NULL }
    # check bivariate/multivariate and change interface back
    else if(mSet$dataSet$cls.num == 2){
      interface$mode <- "bivar"}
    else{
      interface$mode <- "multivar"}

    output$curr_name <- renderText({mSet$dataSet$cls.name})

  }

  mSet$last_mset <<- mset_name

  statsmanager$calculate <- input$statistics
  datamanager$reload <- input$statistics

}, ignoreInit = TRUE)


