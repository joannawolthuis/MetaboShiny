# create listener for what mode we're currently working in (bivariate, multivariate, time series...)
uimanager <- shiny::reactiveValues()

# pmake pca/plsda
shiny::observe({
  if(is.null(uimanager$refresh)){
    NULL # if not makeing anything, nevermind
  }else{
    if(!is.null(mSet)){
      success = F
      try({
        mSet.old <- mSet
        for(do in uimanager$refresh){
          suppressWarnings({
            switch(do,
                   general = {
                     # make sidebar
                     # make pca, plsda, ml(make uimanager do that)
                     # update select input bars with current variable and covariables defined in excel
                     if(is.null(mSet)){
                       interface$mode <- NULL
                     }else{
                       if(is.null(mSet$settings$exp.type)){
                         mSet$settings$exp.type <- "1fb" # one factor, binary class
                       }  
                       
                       shiny::showNotification("Updating interface...")
                       
                       interface$mode <<- mSet$settings$exp.type
                       
                       output$curr_name <- shiny::renderUI(
                         shiny::HTML({
                           expname = mSet$settings$cls.name
                           calc = stringr::str_match(expname, pattern = "(^.*?):")
                           expname = gsub(expname, pattern = "^.*?:", replacement="") 
                           if(!all(is.na(calc[1,]))){
                             calc = calc[,2]
                             subsets = stringr::str_split(expname, ",:")[[1]]
                             subsetPart = paste0(shiny::icon("filter", "fa-xs"),
                                                 " ",
                                                 subsets, collapse="<br>")
                             statsPart = h2(calc)
                           }else{
                             subsetPart = ""
                             statsPart = h2(expname)
                           }
                           paste0("<br>", 
                                  statsPart, 
                                  subsetPart, 
                                  "<br>",
                                  "<br>")
                         })
                       )
                       
                       if(any(grepl("RT", colnames(mSet$dataSet$norm)))){
                         for(picker in c("add", "iso")){
                           output[[paste0(picker, "_rt_ui")]] <- shiny::renderUI({
                             list(
                               shiny::helpText("Only consider options within a certain retention time window?"),
                               shinyWidgets::switchInput(
                                 inputId = paste0(picker, "_use_rt"),
                                 size = "mini",
                                 onLabel = "Yes", 
                                 offLabel = "No", 
                                 value = TRUE
                               ),
                               shiny::conditionalPanel(gsubfn::fn$paste("input.$picker_use_rt == true"),
                                                       shiny::numericInput(paste0(picker, "_rt_perc") ,
                                                                           label = "Max. retention time error margin:", 
                                                                           value = 0.1, 
                                                                           min = 0,
                                                                           max = 100,
                                                                           width = "30%")
                               ) 
                             )
                           })
                         }
                       }
                       # --- MZ PICKER ---
                       neededForChoices = lapply(colnames(mSet$dataSet$norm), function(mzfull){
                         if(grepl("RT", mzfull)){
                           split = stringr::str_split(mzfull, "RT")[[1]]
                           mz = split[1]
                           rt = split[2]
                         }else{
                           mz = mzfull
                           rt = ""
                         }
                         list(full = mzfull, mz = mz, rt = rt)
                       })
                       
                       allMz = lapply(neededForChoices, function(x) x$full)
                       names(allMz) = lapply(neededForChoices, function(x) paste(x$mz, "m/z"))
                       subtext = lapply(neededForChoices, function(x) x$rt)
                       newOrder = order(unlist(allMz))
                       
                       shinyWidgets::updatePickerInput(session, 
                                                       "curr_mz",
                                                       choices = append(allMz[newOrder], 
                                                                        list("select a m/z")), 
                                                       choicesOpt = list(subtext = c(subtext[newOrder], "select a m/z"),
                                                                         icon = c(rep('',length(subtext)), "fa-cat"),
                                                                         style = c(rep('text-align:center;',length(subtext) + 1))),
                                                       selected = "fa-cat"
                       )
                       
                       shinyWidgets::updatePickerInput(session, 
                                                       "ml_mzs",
                                                       choices = colnames(mSet$dataSet$norm), 
                                                       choicesOpt = list(subtext = c(subtext[newOrder], "select a m/z"),
                                                                         icon = c(rep('',length(subtext)), "fa-cat"),
                                                                         style = c(rep('text-align:center;',length(subtext) + 1))),
                                                       selected = "fa-cat"
                       )
                       
                       # -----------------
                       shiny::updateNavbarPage(session, "statistics", selected = "inf")
                       origcount = mSet$settings$orig.count 
                       output$samp_count <- shiny::renderText({
                         paste0(as.character(nrow(mSet$dataSet$norm)),
                                if(nrow(mSet$dataSet$norm) == origcount) "" else paste0("/",as.character(origcount)))
                       })
                       
                       output$mz_count <- shiny::renderText({
                         as.character(ncol(mSet$dataSet$norm))
                       })
                       
                       if(length(mSet$storage) > 0){
                         storeNames = c(names(mSet$storage)[sapply(mSet$storage, function(x) "settings" %in% names(x))])
                         subtext = sapply(storeNames, function(storeName){
                           sz=format(object.size(mSet$storage[[storeName]]), unit="Mb")
                           if(sz != "0 Mb"){
                             sz
                           }else{
                             ""
                           }
                         })
                       }else{
                         storeNames=c(" ")
                         subtext=c(" ")
                       }
                       shinyWidgets::updatePickerInput(session, 
                                                       "storage_choice",
                                                       selected = if(!is.null(storeNames[[1]])) storeNames else c(" "), 
                                                       choices = if(length(storeNames) > 0) storeNames else c(" "),
                                                       choicesOpt = list(subtext = subtext,
                                                                         style = c(rep('text-align:center;',length(subtext))))
                       )
                       
                       usesCovars <- paste0(c("stats","time","shape","col","txt","subset","fill","ml_mistake"), "_var")
                       lapply(usesCovars, function(inputId){
                         shiny::updateSelectizeInput(session, inputId,  
                                                  choices = c(#"label", 
                                                    colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, 
                                                                                              MARGIN = 2, 
                                                                                              function(col) length(unique(col)) < gbl$constants$max.cols))]),
                                                  selected = mSet$settings$exp.var)
                         
                       })
                     }
                     
                     shiny::updateSelectizeInput(session, "stats_type", 
                                              selected = {
                                                if(grepl(mSet$settings$exp.type, pattern = "^1f\\w")){
                                                  gsub(mSet$settings$exp.type, pattern="^1f\\w", "1f")
                                                }else{
                                                  mSet$settings$exp.type
                                                }
                                              })
                     
                     shiny::updateSelectizeInput(session, "stats_type", 
                                              choices = if(!any(duplicated(mSet$dataSet$covars$individual))){
                                                list("one factor"="1f", 
                                                     "two factors"="2f")
                                              }else{
                                                list("one factor"="1f", 
                                                     "two factors"="2f",
                                                     "time series"="t",
                                                     "time series + one factor"="t1f")
                                              })
                     
                     shiny::setProgress(0.7)
                     
                     if(mSet$metshiParams$prematched){
                       search_button$on <- FALSE
                     }else{
                       search_button$on <- TRUE}
                     
                     ml_queue$jobs = list()
                     
                     shiny::updateNavbarPage(session, "statistics", selected = "inf")
                   },
                   wordcloud = {
                     wordcloud_filters = file.path(lcl$paths$work_dir, "wordcloud")
                     if(dir.exists(wordcloud_filters)){
                       filter_files = list.files(wordcloud_filters, full.names = T)
                       items = lapply(filter_files, function(file){
                         data.table::fread(file)$word
                       })
                       itemNames = gsub(pattern = "\\.csv", replacement = "", x = basename(filter_files))
                       names(items) <- itemNames
                       gbl$vectors$wordcloud$filters <<- append(gbl$vectors$wordcloud$filters, items)
                       shiny::updateTabsetPanel(session, "wordclouds", selected = "plot")
                       shiny::updateSelectizeInput(session, 
                                                inputId = "wordcloud_filter", 
                                                #selected="stopwords",
                                                choices = c("stopwords","metabolomics","default", itemNames))
                     }
                   },
                   statspicker = {
                     output$stats_picker <- shiny::renderUI({
                       if(input$stats_type != "t"){
                         shiny::selectizeInput("stats_var", 
                                               label="Experimental variables:", 
                                               multiple = if(input$stats_type == "2f") T else F,
                                               selected = mSet$settings$exp.var, 
                                               choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, 
                                                                                                              MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]),
                                               options = list(maxItems = 2),
                                               width = "80%") 
                       }else{
                         list()
                       }
                     }) 
                   },
                   combi = {
                     combi_opts = setdiff(names(mSet$analSet)[sapply(mSet$analSet,function(anal) any(grepl("\\.mat", names(anal))))], "combi")
                     if(length(combi_opts)>1){
                       shiny::updateSelectizeInput(session, "combi_anal1",  
                                                choices = combi_opts,
                                                selected = combi_opts[1])
                       shiny::updateSelectizeInput(session, "combi_anal2",  
                                                choices = combi_opts,
                                                selected = combi_opts[2])
                     }
                   },
                   venn = {
                     shiny::updateSelectizeInput(session, 
                                                 "intersect_venn", 
                                                 choices =  names(lcl$vectors$venn_lists))
                   },
                   corr = {
                     output$jqui_ui <- shiny::renderUI(suppressWarnings(shinyjqui::orderInput(inputId = 'corr_seq',
                                                                                              label = 'Drag panels to select pattern for correlation (low-high)', 
                                                                                              items = {
                                                                                                lvls = levels(mSet$dataSet$cls)
                                                                                                numconv = as.numeric(as.character(lvls))
                                                                                                if(all(!is.na(numconv))){
                                                                                                  order = order(numconv)
                                                                                                }else{
                                                                                                  order = order(as.character(lvls))
                                                                                                }
                                                                                                as.character(lvls)[order]
                                                                                              }))
                     )
                     storeNames = c(names(mSet$storage)[sapply(mSet$storage, function(x) "settings" %in% names(x))])
                     shiny::updateSelectizeInput(session, "storage_choice", 
                                              choices = if(!is.null(storeNames[[1]])) storeNames else c()
                     ) 
                   },
                   ml = {
                     shiny::updateSelectizeInput(session, "ml_include_covars", 
                                              choices = c(colnames(mSet$dataSet$covars)[!(colnames(mSet$dataSet$covars) %in% c("label", "sample", "individual"))]))
                     
                     shiny::updateSelectizeInput(session, "ml_batch_covars", 
                                              choices = c(colnames(mSet$dataSet$covars)[!(colnames(mSet$dataSet$covars) %in% c("label", "sample", "individual"))]))
                     # --- ML PCA ---
                     npc = nrow(mSet$dataSet$norm) - 1
                     shiny::updateSliderInput(session, "ml_keep_pcs", max = npc)
                     
                     if(length(mSet$analSet$ml) > 0){
                       
                       shiny::showTab(session = session, 
                                      inputId = "ml2", 
                                      target = "res")
                       
                       ###
                       data = mSet$analSet$ml[[mSet$analSet$ml$last$method]][[mSet$analSet$ml$last$name]]
                       classes = if(!is.null(data$res$prediction)){
                         colnames(data$res$prediction)
                       }else{
                         colnames(data$res[[1]]$prediction)
                       }
                       
                       # try both options and check higher AUC
                       perf.per.posclass = sapply(classes, function(pos.class.test){
                         perf = getMLperformance(ml_res = data$res[[1]], 
                                                 pos.class = pos.class.test,
                                                 x.metric = input$ml_plot_x,
                                                 y.metric = input$ml_plot_y)
                         t = perf$coords
                         pracma::trapz(t[`Test set` == "Test"]$x, 
                                       t[`Test set` == "Test"]$y)
                       })
                       pos.class <- names(which.max(perf.per.posclass))
                       
                       shiny::updateSelectizeInput(session, "ml_plot_posclass", 
                                                choices = classes, 
                                                selected = pos.class)
                       ###
                       
                       choices = c()
                       methods <- setdiff(names(mSet$analSet$ml), "last")
                       for(method in methods){
                         model.names = names(mSet$analSet$ml[[method]])
                         choices <- c(choices, paste0(method, " - ", paste0(model.names)))
                       }
                       
                       shiny::updateSelectizeInput(session, "ml_samp_distr", 
                                                choices = c(" ", choices))
                       
                       shiny::updateSelectizeInput(session, "show_which_ml", 
                                                choices = choices, 
                                                selected = paste0(mSet$analSet$ml$last$method, 
                                                                  " - ", 
                                                                  mSet$analSet$ml$last$name),
                                                server = T)
                       
                     }else{
                       shiny::hideTab(session = session, inputId = "ml2", target = "res")
                     }
                   },
                   tt = {
                     if("V" %in% colnames(mSet$analSet$tt$sig.mat)){
                       shiny::updateCheckboxInput(session, "tt_nonpar", value = T)
                     }else{
                       shiny::updateCheckboxInput(session, "tt_nonpar", value = F)
                     }
                   },
                   heatmap = 
                     {
                       NULL
                     }, 
                   power = {
                     if(grepl(mSet$settings$exp.type, pattern = "1f.")){
                       pairs = MetaboShiny::expand.grid.unique(levels(mSet$dataSet$cls), 
                                                               levels(mSet$dataSet$cls))
                       pairs = paste(pairs[,1], "vs.", pairs[,2])
                       shiny::updateSelectizeInput(session, inputId = "power_comps", choices = pairs, selected = pairs[1])
                     }
                   }
            )
          })
        }
        success = T
      })
      if(!success){
        MetaboShiny::metshiAlert("Changing UI failed!")
        mSet <- mSet.old
      }
    }
    uimanager$refresh <- NULL # set making to 'off'
  }
})