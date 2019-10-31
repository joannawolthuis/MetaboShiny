# triggers when the 'change variable' dropdown menu is filled and button is clicked
observeEvent(input$change_cls, {

  # check if previous analysis storage already exists, if not, make it
  mset_name = mSet$dataSet$cls.name

  # save previous analyses (should be usable in venn diagram later)
  mSet$storage[[mset_name]] <- list(analysis = mSet$analSet)

  lcl$constants$last_mset <- mset_name

  # adjust name of experimental variable
  if(grepl(mSet$dataSet$cls.name, pattern = ":")){
    subset_name <- paste0(":", gsub(mSet$dataSet$cls.name, pattern = ".*:", replacement = ""))
  }else{
    subset_name <- ""
  }
  
  mSet <- switch(input$stats_type,
         "1f"={
           # change current variable of interest to user pick from covars table
           mSet$dataSet$cls <- as.factor(mSet$dataSet$covars[,input$stats_var, with=F][[1]])
           # adjust bivariate/multivariate (2, >2)...
           mSet$dataSet$cls.num <- length(levels(mSet$dataSet$cls))
           mSet$dataSet$cls.name <- paste0(input$stats_var, subset_name)
           # - - - 
           mSet
         },
         "2f"={
           # facB should be time if it is there...
           time.check = grepl("time", input$stats_var)
           is.time = any(time.check)
           if(is.time){
             print("time series potential")
             idx2 = which(time.check)
             idx1 = setdiff(c(1,2), idx2)
           }else{
             idx1=1
             idx2=2
           }
           mSet$dataSet$facA <- as.factor(mSet$dataSet$covars[,input$stats_var, with=F][[idx1]])
           mSet$dataSet$facB <- as.factor(mSet$dataSet$covars[,input$stats_var, with=F][[idx2]])
           mSet$dataSet$facA.lbl <- input$stats_var[idx1]
           mSet$dataSet$facB.lbl <- input$stats_var[idx2]
           mSet$dataSet$cls.name <- paste0(paste0(input$stats_var[c(idx1,idx2)],collapse="+"), subset_name)
           # ONLY ANOVA2
           # - - - 
           mSet
         },
         "t"={
           mSet <- MetaboAnalystR::SetDesignType(mSet, "time0")
           mSet$dataSet$sbj <- as.factor(mSet$dataSet$covars$individual)
           if(!any(duplicated(mSet$dataSet$sbj))){
             print("Won't work, need multiple of the same sample in the 'individual' metadata column!")
             return(NULL)
           }
           mSet$dataSet$time.fac <- as.factor(mSet$dataSet$covars[,input$time_var, with=F][[1]])
           mSet$dataSet$cls.name <- paste0("time", subset_name)
           
           # - - - 
           mSet
         },
         "t1f"={
           mSet <- MetaboAnalystR::SetDesignType(mSet, "time")
           if(!any(duplicated(as.factor(mSet$dataSet$covars$individual)))){
             print("Won't work, need multiple of the same sample in the 'individual' metadata column!")
             return(NULL)
           }
           mSet$dataSet$facA <- as.factor(mSet$dataSet$covars[,input$stats_var, with=F][[1]])
           mSet$dataSet$facB <- as.factor(mSet$dataSet$covars[,input$time_var, with=F][[1]])
           mSet$dataSet$facA.lbl <- input$stats_var
           mSet$dataSet$facB.lbl <- "time"
           mSet$dataSet$exp.fac <- mSet$dataSet$facA
           mSet$dataSet$time.fac <- mSet$dataSet$facB
           mSet$dataSet$cls.name <- paste0(paste0(input$stats_var,"+time"), subset_name)
           # - - - 
           mSet
         })
  

  if(mSet$dataSet$cls.name %in% names(mSet$storage)){
    mSet$analSet <- mSet$storage[[input$stats_var]]$analysis
  }else{
    # remove old analSet
    mSet$analSet <- NULL
  }
  
  mSet$dataSet$exp.type <- input$stats_type

  mSet <<- mSet
  
  interface$mode <- input$stats_type
  datamanager$reload <- "general"
  
  shiny::updateNavbarPage(session, "statistics", selected = "inf")

})
