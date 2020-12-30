# print current compound in sidebar
observe({
  shinyWidgets::updatePickerInput(session, 
                                  "curr_mz", 
                                  selected = my_selection$mz)
})

shiny::observeEvent(input$undo_match_filt, {
  result_filters$add <- list(positive = c(), negative = c())
  result_filters$db <- c()
  result_filters$iso <- c()
  search$go <- T
}) 

output$curr_add <- shiny::renderText({
  paste0(unlist(result_filters$add), collapse=", ")
})

output$curr_iso <- shiny::renderText({
  paste0(result_filters$iso, collapse=", ")
})

output$curr_db <- shiny::renderText({
  paste0(result_filters$db, collapse=", ")
})

observe({
  if(my_selection$mz != ""){
    for(pie in c("add", "iso","db")){
      if(pie == "add"){
        mzMode =if(grepl(my_selection$mz, pattern = "-")) "negative" else "positive"
        result_filters$add[[mzMode]] <- result_filters$add[[mzMode]][!is.na(result_filters$add[[mzMode]])]
      }else{
        result_filters[[pie]] <- result_filters[[pie]][!is.na(result_filters[[pie]])]
      }
    }   
  }
})

shiny::observeEvent(input$show_iso_labels, {
  if( my_selection$mz != "" & nrow(shown_matches$forward_full) > 0){
    search$go <- TRUE
  }
})

shiny::observe({
  if(my_selection$mz != ""){
    if(mSet$metshiParams$prematched){
      search$go <- TRUE
      my_selection$name <- ""
    }
    # update star logo
    shinyWidgets::updatePrettyToggle(session, 
                                     "star_mz",
                                     value =  mSet$report$mzStarred[my_selection$mz]$star)
  }
})

shiny::observe({
  if(my_selection$name != ""){
    if(nrow(shown_matches$forward_full) > 0 ){
      subsec = data.table::as.data.table(shown_matches$forward_full)[compoundname == my_selection$name]
      
      if(grepl("SYNONYMS:", x = subsec$description)){
        has_syn = T
        subsec$source <- paste0("synonymSEPERATOR", subsec$source)
        subsec$description <- gsub(subsec$description, pattern = "SYNONYMS: ", replacement="")
      }else{
        has_syn = F
      }
      
      subsec <- subsec[, .(compoundname, 
                           source = strsplit(source, split = "SEPERATOR")[[1]], 
                           structure = structure, 
                           description = strsplit(description, split = "SEPERATOR")[[1]]
      )
      ]
      
      subsec <- aggregate(subsec, by = list(subsec$source), FUN=function(x) paste0(unique(x), collapse="."))
      
      keep <- which(trimws(subsec$description) %not in% c("","Unknown","unknown", " ",
                                                          "Involved in pathways: . More specifically: . Also associated with compound classes:"))
      subsec <- subsec[keep,]
      
      if(has_syn){
        subsec <- subsec[order(as.numeric(grepl(subsec$source, pattern = "synonym")), decreasing = T),]
      }
      
      if(nrow(subsec) > 0){
        
        # render descriptions seperately
        output$desc_ui <- shiny::renderUI({
          
          lapply(1:nrow(subsec), function(i){
            
            row = subsec[i,]
            
            # icon(s) in one row
            db = row$source
            desc_id = paste0("curr_desc_", db)
            desc = row$description
            #output[[desc_id]] <- shiny::renderText({desc})
            
            if(db == "synonym"){
              ui = shiny::fluidRow(align="center", 
                                   tags$h3("Synonyms:"),
                                   helpText(desc),
                                   shiny::br()
              )
            }else{
              id = gbl$constants$db.build.info[[db]]$image_id
              address = unlist(sapply(gbl$constants$images, function(item) if(item$name == id) item$path else NULL))
              ui = shiny::fluidRow(align="center", 
                                   shiny::tags$button(
                                     id = paste0(db, "_copy_id"),
                                     class = "btn btn-default action-button shiny-bound-input",
                                     shiny::img(src = address,
                                                height = "30px"),
                                     style = "vertical-align: middle;border-radius: 0px;border-width: 0px;background-color: #ff000000;"
                                   ),
                                   shiny::br(),
                                   helpText(desc),
                                   shiny::br()
              )
              
              shiny::observeEvent(input[[paste0(db, "_copy_id")]], {
                shiny::showNotification("Saving database identifier to clipboard!")
                dbID = stringr::str_match(desc, "Database ID: (.*?). ")[,2]
                clipr::write_clip(dbID, allow_non_interactive = TRUE)
                shiny::updateTextInput(session,
                                       "wordcloud_searchTerm",
                                       value = dbID)
              })
            }
            return(ui)
          })
        })  
        
      }else{
        output$desc_ui <- shiny::renderUI({
          helpText("No additional info available!")
        })
      }
    }
  }
})

output$curr_struct <- renderPlot({
  width=min(c(300, shiny::reactiveValuesToList(session$clientData)$output_empty4_width))
  plot_mol(my_selection$struct,
           style = "cow",
           width=width, 
           height=width)
},
height=reactive(min(c(300, shiny::reactiveValuesToList(session$clientData)$output_empty4_width))),
width=reactive(min(c(300, shiny::reactiveValuesToList(session$clientData)$output_empty4_width)))
)

output$curr_formula <- shiny::renderUI({
  shiny::tags$div(
    HTML(gsub(my_selection$form,pattern="(\\d+)", replacement="<sub>\\1</sub>",perl = T))
  )})

shiny::observeEvent(input$curr_mz, {
  if(input$curr_mz %in% colnames(mSet$dataSet$norm)){
    my_selection$mz <- input$curr_mz   
    plotmanager$make <- "summary"
  }
})

shiny::observeEvent(input$subset_var, {
  lvls = levels(as.factor(mSet$dataSet$covars[[input$subset_var]]))
  subtext = sapply(lvls, function(lv){
    samp.count = sum(unlist(mSet$dataSet$covars[[input$subset_var]]) == lv)
    paste(samp.count, "samples")
  })
  shinyWidgets::updatePickerInput(session, 
                                  "subset_group",
                                  selected = if(!is.null(lvls[1])) lvls[1] else c(" "), 
                                  choices = if(length(lvls) > 0) lvls else c(" "),
                                  choicesOpt = list(subtext = subtext,
                                                    style = c(rep('text-align:center;',
                                                                  length(subtext))))
  )
},ignoreInit = T)