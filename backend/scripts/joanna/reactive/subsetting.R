observeEvent(input$change_subset, {

  # save previous
  mset_name <- get_mset_name(mainvar = mSet$dataSet$cls.name,
                             subsetvar = NULL,
                             subsetgroups = NULL)
  mSet$storage[[mset_name]] <<- list(data = mSet$dataSet,
                                     analysis = mSet$analSet)

  # make new subset
  mset_name <- get_mset_name(mainvar = gsub(mSet$dataSet$cls.name,
                                            pattern = ":.*", replacement = ""),
                             subsetvar = input$subset_var,
                             subsetgroups = input$subset_group)

  if(mset_name %in% names(mSet$storage)){
    mSet$dataSet <<- mSet$storage[[mset_name]]$data
    mSet$analSet <<- mSet$storage[[mset_name]]$analysis
  }else{
    mSet$dataSet$cls.name <<- mset_name
    keep.samples <<- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    mSet$dataSet$covars <<- mSet$dataSet$covars[sample %in% keep.samples]
    keep.log.proc <<- rownames(mSet$dataSet$proc) %in% keep.samples
    mSet$dataSet$proc <<- mSet$dataSet$proc[keep.log.proc,]
    keep.log.norm <<- rownames(mSet$dataSet$norm) %in% keep.samples
    mSet$dataSet$norm <<- mSet$dataSet$norm[keep.log.norm,]
    keep.log.preproc <<- rownames(mSet$dataSet$preproc) %in% keep.samples
    mSet$dataSet$preproc <<- mSet$dataSet$preproc[keep.log.preproc,]
    mSet$dataSet$cls <<- droplevels(mSet$dataSet$cls[keep.log.norm])
    mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))

    if("facA" %in% names(mSet$dataSet)){
      mSet$dataSet$facA <<- mSet$dataSet$facA[keep.log.norm]
      mSet$dataSet$facB <<- mSet$dataSet$facB[keep.log.norm]
      mSet$dataSet$time.fac <<- mSet$dataSet$time.fac[keep.log.norm]
      mSet$dataSet$exp.fac <<- mSet$dataSet$exp.fac[keep.log.norm]
    }
  }
  local$last_mset <<- mset_name
  mSet$analSet <<- NULL

  covars <- colnames(mSet$dataSet$covars)
  subsettable.covars <- covars[which(sapply(covars, function(x){
    lvls <- unlist(mSet$dataSet$covars[,..x])
    level.count <- length(levels(as.factor(lvls)))
    keep <- level.count > 1 & level.count < length(lvls)
    # - -
    keep
  }))]
  updateSelectInput(session, "subset_var", choices = subsettable.covars)
  output$curr_name <- renderText({mSet$dataSet$cls.name})
  if(mSet$dataSet$cls.num <= 1){
    interface$mode <- NULL }
  else if(mSet$dataSet$cls.num == 2){
    interface$mode <- "bivar"}
  else{
    interface$mode <- "multivar"}

  if(!("pca" %in% names(mSet$analSet))){
    statsmanager$calculate <- "pca"
  }

  datamanager$reload <- "pca"
  })

observeEvent(input$reset_subset, {

  # save previous
  mSet$storage[[mSet$dataSet$cls.name]] <<- list(data = mSet$dataSet,
                                                 analysis = mSet$analSet)

  # load previous
  mset_name <- get_mset_name(mainvar = gsub(mSet$dataSet$cls.name, pattern = ":.*", replacement = ""),
                             subsetvar = NULL,
                             subsetgroups = NULL)

  if(mset_name %in% names(mSet$storage)){
    mSet$dataSet <<- mSet$storage[[mset_name]]$data
    mSet$analSet <<- mSet$storage[[mset_name]]$analysis
  }else{
    main.var <- gsub(mSet$dataSet$cls.name, pattern = ":.*", replacement = "")

    # get any non subsetted mset
    full.msets <- grep(names(mSet$storage), pattern = ":", invert = T)
    keep <- which(sapply(full.msets, function(mset){
      "data" %in% names(mSet$storage[[mset]])
    }))[[1]]
    mSet$dataSet <<- mSet$storage[[keep]]$data
    mSet$analSet <<- mSet$storage[[keep]]$analysis
    # change current variable of interest to user pick from covars table
    mSet$dataSet$cls <<- as.factor(mSet$dataSet$covars[,main.var, with=F][[1]])
    levels(mSet$dataSet$cls) <<- droplevels(mSet$dataSet$cls)

    # adjust bivariate/multivariate (2, >2)...
    mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
    mSet$dataSet$cls.name <<- mset_name
    mSet$analSet <- NULL
  }

  local$last_mset <<- mset_name

  covars <- colnames(mSet$dataSet$covars)
  subsettable.covars <- covars[which(sapply(covars, function(x){
    lvls <- unlist(mSet$dataSet$covars[,..x])
    level.count <- length(levels(as.factor(lvls)))
    keep <- level.count > 1 & level.count < length(lvls)
    # - -
    keep
  }))]

  updateSelectInput(session, "subset_var", choices = subsettable.covars)

  output$curr_name <- renderText({mSet$dataSet$cls.name})

  # reload current plot
  # for(tabgroup in c("dimred", "permz", "overview")){
  #   if(tabgroup %in% names(input)){
  #     statsmanager$calculate <- input[[tabgroup]]
  #     datamanager$reload <- input[[tabgroup]]
  #   }
  # }
  datamanager$reload <<- "general"

  updateNavbarPage(session, "statistics", selected = "inf")

})

observeEvent(input$subset_var, {
  lvls = levels(as.factor(mSet$dataSet$covars[[input$subset_var]]))
  updateSelectizeInput(session, "subset_group", choices = lvls)
},ignoreInit = T)
