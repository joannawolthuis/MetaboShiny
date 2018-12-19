observeEvent(input$change_subset, {
  
  # save previous
  mset_name <- get_mset_name(mainvar = mSet$dataSet$cls.name,
                             subsetvar = NULL,
                             subsetgroups = NULL)
  mSet$storage[[mset_name]] <<- list(data = mSet$dataSet,
                                     analysis = mSet$analSet)
  
  # make new subset    
  mset_name <- get_mset_name(mainvar = gsub(mSet$dataSet$cls.name, pattern = ":.*", replacement = ""),
                             subsetvar = input$subset_var,
                             subsetgroups = input$subset_group)
  
  if(mset_name %in% names(mSet$storage)){
    mSet$dataSet <<- mSet$storage[[mset_name]]$data
    mSet$analSet <<- mSet$storage[[mset_name]]$analysis
  }else{
    mSet$dataSet$cls.name <<- mset_name
    keep.samples <<- mSet$dataSet$covars$sample[which(mSet$dataSet$covars[[input$subset_var]] %in% input$subset_group)]
    mSet$dataSet$covars <<- mSet$dataSet$covars[sample %in% keep.samples]
    keep.log.procr <<- rownames(mSet$dataSet$procr) %in% keep.samples
    mSet$dataSet$procr <<- mSet$dataSet$procr[keep.log.procr,]
    keep.log.norm <<- rownames(mSet$dataSet$norm) %in% keep.samples
    mSet$dataSet$norm <<- mSet$dataSet$norm[keep.log.norm,]
    mSet$dataSet$cls <<- mSet$dataSet$cls[keep.log.norm]
    if("facA" %in% names(mSet$dataSet)){
      mSet$dataSet$facA <<- mSet$dataSet$facA[keep.log.norm]
      mSet$dataSet$facB <<- mSet$dataSet$facB[keep.log.norm]
      mSet$dataSet$time.fac <<- mSet$dataSet$time.fac[keep.log.norm]
      mSet$dataSet$exp.fac <<- mSet$dataSet$exp.fac[keep.log.norm]
    } 
  }
  global$constants$last_mset <<- mset_name
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
  datamanager$reload <- input$statistics # reload current one
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
    # adjust bivariate/multivariate (2, >2)...
    mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
    mSet$dataSet$cls.name <<- mset_name
    mSet$analSet <- NULL
  }
  global$constants$last_mset <<- mset_name
  
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
  
  datamanager$reload <- input$statistics # reload current one
  
})

observeEvent(input$subset_var, {
  lvls = levels(as.factor(mSet$dataSet$covars[[input$subset_var]]))
  updateSelectizeInput(session, "subset_group", choices = lvls)
},ignoreInit = T)