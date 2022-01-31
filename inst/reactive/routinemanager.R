# paperRoutine
paper_routine <- list(
  list(type = "change",
       settings = list(stats_var = "group",
                       stats_type = "1f")),
  # --------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold01", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold02", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold03", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold04", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold05", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold06", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold07", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold08", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold09", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold10", subset_group = "yes"
  )),
  list(type = "analyse",
       settings = list(
         analysis = "tt"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "fc"
       )),
  list(type = "analyse",
       settings = list(
         analysis = "combi"
       )),
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )))
routines = list(tasks = paper_routine)

# -------------------------
routines <- shiny::reactiveValues(tasks = paper_routine)

routinemanager <- shiny::reactive({
  if(length(routines$tasks) > 0 & input$start_routine){
    shiny::withProgress(max = length(routines$tasks), {
      
      pb = pbapply::startpb(max = length(routines$tasks))
      
      for(i in 1:length(routines$tasks)){
        print(i)
        task = routines$tasks[[i]]
        print(task$type)
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
                 results <- runStats(mSet,
                                  input = input,
                                  lcl = lcl,
                                  analysis = task$settings$analysis)
                 mSet <- results$mSet
                 lcl <- results$lcl
               })

        # ----------
        #pbapply::setpb(pb, i)
        #try({
        #  shiny::setProgress(value = i)
        #}, silent = T)
      } 
    }
    )
  }
})