routines = list(tasks = routine_tt_fc)
routines <- shiny::reactiveValues(tasks = list())

routinemanager <- shiny::reactive({
  if(length(routines$tasks) > 0 & input$start_routine){
    shiny::withProgress(max = length(routines$tasks), {
      
      mSet_before = mSet
      pb = pbapply::startpb(max = length(routines$tasks))
      last_task = list()
      
      for(i in 1:length(routines$tasks)){
        mSet_last <<- mSet
        success = F
        
        task = routines$tasks[[i]]
        last_task <<- task
        covars <<- mSet$dataSet$covars
        settings <<- mSet$settings
        task_type = if(task$type %in% c("subset","change","load")) "update" else "analyse"
        # do task...
        done = switch(task_type,
                      update = {
                        if(task$type == "subset"){
                          input$subset_group <- task$settings$subset_group
                          input$subset_var <- task$settings$subset_var
                        }else if(task$type == "load"){
                          
                          input$storage_choice <- task$settings$storage_choice
                        }else{
                          input$stats_var <- task$settings$stats_var
                          input$stats_type <- task$settings$stats_type
                        }
                        updated <- doUpdate(mSet,
                                            input = input,
                                            lcl = lcl,
                                            do = task$type)
                        mSet <- updated$mSet
                        lcl <- updated$lcl
                      },
                      analyse = {
                        if(!("multirank_yes" %in% names(task))){
                          multirank_yes = data.frame()
                        }else{
                          multirank_yes = task$multirank_yes
                        }
                        #try({
                        results <- runStats(mSet,
                                            input = input,
                                            lcl = lcl,
                                            analysis = task$settings$analysis,
                                            multirank_yes = multirank_yes,
                                            ml_queue = ml_queue,
                                            cl=NULL)
                          mSet <- results$mSet
                        lcl <- results$lcl  
                        #})
                      })
        if("dataSet" %in% names(mSet)){
          success <<- T
        }
        
        if(success){
          print("yay!")
        }else{
          print("nay")
          stop("oh no")
          mSet <<- mSet_last
        }
      } 
        #})
    })
  }
})