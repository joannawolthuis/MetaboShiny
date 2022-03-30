# input$stats_var <- "group"
# input$stats_type <- "1f"
# input$tt_nonpar=F
# input$tt_p_thresh=0.05
# input$tt_eqvar = F
# input$tt_multi_test="fdr"
# input$fc_thresh=1.2
# input$combi_anal1 = "fc"
# input$combi_anal2 = "tt"
# input$combi_anal1_var = "log2(FC)"
# input$combi_anal2_var = "-log10(p)"
# input$combi_dist_metric="multiplication"
# input$ml_use_slurm = TRUE
# input$ml_slurm_job_mem = "2G"

# paperRoutine
paper_routine_10fold <- list(
  list(type = "change",
       settings = list(stats_var = "group",
                       stats_type = "1f")),
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
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

paper_routine_5fold <- list(
  list(type = "change",
       settings = list(stats_var = "group",
                       stats_type = "1f")),
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold1", subset_group = "yes"
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold2", subset_group = "yes"
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold3", subset_group = "yes"
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold4", subset_group = "yes"
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
  # -------------------------
  list(type = "subset", settings = list(
    subset_var = "in_featsel_fold5", subset_group = "yes"
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
  #list(type = "subset", settings = list(
  #  subset_var = "averaged", subset_group = "TRUE"
  #)),
  # -------------------------
  list(type = "load",
       settings = list(
         storage_choice = "group"
       )))
routines = list(tasks = paper_routine_5fold)

# -------------------------
routines <- shiny::reactiveValues(tasks = paper_routine_5fold)

routinemanager <- shiny::reactive({
  if(length(routines$tasks) > 0 & input$start_routine){
    shiny::withProgress(max = length(routines$tasks), {
      
      pb = pbapply::startpb(max = length(routines$tasks))
      
      for(i in 1:length(routines$tasks)){
        task = routines$tasks[[i]]
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