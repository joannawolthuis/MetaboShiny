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
                       if(is.null(mSet$dataSet$exp.type)){
                         mSet$dataSet$exp.type <- "1f" # one factor, binary class
                       }  
                       
                       shiny::showNotification("Updating interface...")
                       
                       interface$mode <<- mSet$dataSet$exp.type
                       
                       output$curr_name <- shiny::renderUI(
                         shiny::HTML({
                           expname = mSet$dataSet$cls.name
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
                       
                       # -----------------
                       shiny::updateNavbarPage(session, "statistics", selected = "inf")
                       origcount = mSet$settings$orig.count 
                       output$samp_count <- shiny::renderText({
                         paste0(as.character(nrow(mSet$dataSet$norm)),
                                if(nrow(mSet$dataSet$norm) == origcount) "" else paste0("/",as.character(origcount)))
                       })
                       storeNames = c(names(mSet$storage)[sapply(mSet$storage, function(x) "settings" %in% names(x))])
                       shiny::updateSelectInput(session, "storage_choice", 
                                                choices = if(!is.null(storeNames[[1]])) storeNames else c()
                       )
                       usesCovars <- paste0(c("stats","time","shape","col","txt","subset"), "_var")
                       lapply(usesCovars, function(inputId){
                         shiny::updateSelectInput(session, inputId,  
                                                  choices = c("label", 
                                                              colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, 
                                                                                                        MARGIN = 2, 
                                                                                                        function(col) length(unique(col)) < gbl$constants$max.cols))]))
                         
                       })
                     }
                     
                     shiny::updateSelectInput(session, "stats_type", 
                                              selected = {
                                                if(grepl(mSet$dataSet$exp.type, pattern = "^1f\\w")){
                                                  gsub(mSet$dataSet$exp.type, pattern="^1f\\w", "1f")
                                                }else{
                                                  mSet$dataSet$exp.type
                                                }
                                              })
                     
                     shiny::updateSelectInput(session, "stats_type", 
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
                       shiny::updateSelectInput(session, 
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
                                               selected = mSet$dataSet$exp.var, 
                                               choices = c("label", colnames(mSet$dataSet$covars)[which(apply(mSet$dataSet$covars, 
                                                                                                              MARGIN = 2, function(col) length(unique(col)) < gbl$constants$max.cols))]),
                                               options = list(maxItems = 2),
                                               width = "80%") 
                       }else{
                         list()
                       }
                     }) 
                   },
                   vennrich = {
                     NULL
                   },
                   pattern = {
                     output$jqui_ui <- shiny::renderUI(suppressWarnings(shinyjqui::orderInput(inputId = 'pattern_seq',
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
                     shiny::updateSelectInput(session, "storage_choice", 
                                              choices = if(!is.null(storeNames[[1]])) storeNames else c()
                     ) 
                   },
                   ml = {
                     shiny::updateSelectInput(session, "ml_include_covars", 
                                              choices = c(colnames(mSet$dataSet$covars)[!(colnames(mSet$dataSet$covars) %in% c("label", "sample", "individual"))]))
                     
                     if("ml" %in% names(mSet$analSet)){
                       shiny::showTab(session = session, 
                                      inputId = "ml2", 
                                      target = "res")
                       
                       choices = c()
                       methods <- setdiff(names(mSet$analSet$ml), "last")
                       for(method in methods){
                         model.names = names(mSet$analSet$ml[[method]])
                         choices <- c(choices, paste0(method, " - ", paste0(model.names)))
                       }
                       shiny::updateSelectInput(session, "show_which_ml", choices = choices, selected = paste0(mSet$analSet$ml$last$method, " - ", mSet$analSet$ml$last$name))
                       
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
                     if(grepl(mSet$dataSet$exp.type, pattern = "1f.")){
                       pairs = MetaboShiny::expand.grid.unique(levels(mSet$dataSet$cls), 
                                                               levels(mSet$dataSet$cls))
                       pairs = paste(pairs[,1], "vs.", pairs[,2])
                       shiny::updateSelectInput(session, inputId = "power_comps", choices = pairs, selected = pairs[1])
                     }
                   }
            )
          })
        }
        success = T
      })
      if(!success){
        MetaboShiny::metshiAlert("Changing UI failed!")
        mSet <<- mSet.old
      }
    }
    uimanager$refresh <- NULL # set making to 'off'
  }
})
