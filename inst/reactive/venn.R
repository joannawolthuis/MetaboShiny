# nonselected

venn_yes = shiny::reactiveValues(start = data.frame(),
                                 now = data.frame())

venn_no = shiny::reactiveValues(
                        start = data.frame(c("a", "b", "c")),
                        now = data.frame(c("a", "b", "c")))

multirank_yes = shiny::reactiveValues(start = data.frame(),
                                 now = data.frame())

multirank_no = shiny::reactiveValues(
  start = data.frame(c("a", "b", "c")),
  now = data.frame(c("a", "b", "c")))

venn_members <- shiny::reactiveValues(mzvals = list())

shiny::observeEvent(input$venn_add, {
  # add to the 'selected' table
  rows <- input$venn_unselected_rows_selected
  # get members and send to members list
  added = venn_no$now[rows,]
  venn_yes$now <<- rbind(venn_yes$now, added)
  venn_no$now <<- venn_no$now[-rows,]
})

shiny::observeEvent(input$venn_remove, {
  # add to the 'selected' table
  rows <- input$venn_selected_rows_selected
  # get members and send to non members list
  removed = venn_yes$now[rows,]
  venn_no$now <<- rbind(venn_yes$now, removed)
  venn_yes$now <<- data.frame(venn_yes$now[-rows,])
})

shiny::observeEvent(input$multirank_add, {
  # add to the 'selected' table
  rows <- input$multirank_unselected_rows_selected
  # get members and send to members list
  added = multirank_no$now[rows,]
  multirank_yes$now <<- rbind(multirank_yes$now, added)
  multirank_no$now <<- multirank_no$now[-rows,]
})

shiny::observeEvent(input$multirank_remove, {
  # add to the 'selected' table
  rows <- input$multirank_selected_rows_selected
  # get members and send to non members list
  removed = multirank_yes$now[rows,]
  multirank_no$now <<- rbind(multirank_yes$now, removed)
  multirank_yes$now <<- data.frame(multirank_yes$now[-rows,])
})

output$venn_threshold_ui <- shiny::renderUI({
  row = input$venn_selected_row_last_clicked
  if(!is.null(row)){
    if(venn_yes$now[row, "threshold"] != "any"){
      thresh = venn_yes$now[row, "threshold"]#">10"
      sign = stringr::str_extract(thresh,pattern = ">|<|=")
      value = as.numeric(gsub(sign, "", thresh))
    }else{
      sign = "="
      value = 0
    }
    shiny::fluidRow(column(width = 4, shiny::selectInput("venn_threshold_sign",
                                                         label="sign:",
                                                         choices = c("=",">","<"),
                                                         multiple = F,
                                                         selected = sign)), 
                    column(width = 4, shiny::numericInput("venn_threshold_value",
                                                          label = "value:",
                                                          min = 0,
                                                          max = 100,
                                                          step = 0.01,
                                                          value = value)),
                    column(width=4, br(),shinyWidgets::circleButton("venn_threshold_set",
                                                               icon = icon("check"),
                                                               size="sm")))
  }else{"Select analysis in bottom table to set threshold!"}
})

shiny::observeEvent(input$venn_threshold_set, {
  rows = input$venn_selected_rows_selected
  sign = input$venn_threshold_sign
  value = input$venn_threshold_value
  venn_yes$now[rows,"threshold"] <- paste0(sign, value)
})

# the 'non-selected' table
output$venn_unselected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames=c("result","threshold"))
  try({
    res = DT::datatable(venn_no$now, 
                        rownames = FALSE, 
                        colnames=c("result","threshold"), 
                        selection = "multiple")
  })
  res
})

# the 'selected' table
output$venn_selected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames=c("result","threshold"))
  try({
    res = DT::datatable(venn_yes$now,rownames = FALSE, colnames=c("result","threshold"), 
                        selection = "multiple")
  })
  res
})

# the 'non-selected' table
output$multirank_unselected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames=c("result","threshold"))
  try({
    res = DT::datatable(multirank_no$now, 
                        rownames = FALSE, 
                        colnames=c("result","threshold"), 
                        selection = "multiple")
  })
  res
})

# the 'selected' table
output$multirank_selected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames=c("result","threshold"))
  try({
    res = DT::datatable(multirank_yes$now,rownames = FALSE, colnames=c("result","threshold"), 
                        selection = "multiple")
  })
  res
})

# triggers on clicking the 'go' button on the venn diagram sidebar panel
shiny::observeEvent(input$venn_build, {
  plotmanager$make <- "venn"
  uimanager$refresh <- "venn"
})

# triggers when users pick which intersecting hits they want
shiny::observeEvent(input$intersect_venn, {

  if(length(input$intersect_venn) == 0){
    lcl$tables$venn_overlap <<- data.frame()
  }else if(length(input$intersect_venn) == 1){

    l = lcl$vectors$venn_lists
    # Get the combinations of names of list elements
    nms <- combn(names(l), 2, FUN = paste0, collapse = "  ~ ", simplify = FALSE)

    # Make the combinations of list elements
    ll <- combn(l, 2, simplify = FALSE)

    # Intersect the list elements
    out <- lapply(ll , function(x) (intersect(x[[1]],
                                              x[[2]])))
    # Output with names
    intersecties <- unique(unlist(out))

    uniqies = lcl$vectors$venn_lists[[input$intersect_venn]]
    overlap <- setdiff(uniqies,intersecties)
    overlap = gsub("\\.$", "-", overlap)
    mSet$analSet$venn <<- list(mzs = overlap)
    lcl$tables$venn_overlap <<- overlap
  }else{
    overlap <- Reduce("intersect", lcl$vectors$venn_lists[input$intersect_venn])
    overlap = gsub("\\.$", "-", overlap)
    mSet$analSet$venn <<- list(mzs = overlap)
    lcl$tables$venn_overlap <<- overlap 
    }

  lcl$tables$venn_overlap <<- data.frame(mz = lcl$tables$venn_overlap)
  rownames(lcl$tables$venn_overlap) <<- lcl$tables$venn_overlap$mz

  # hypergeometric testing...
  if(nrow(lcl$tables$venn_overlap) > 0){
    try({
      pval = switch(as.character(length(input$intersect_venn)),
                    "2" = 1 - phyper(nrow(lcl$tables$venn_overlap),
                                     length(lcl$vectors$venn_lists[[input$intersect_venn[1]]]),
                                     ncol(mSet$dataSet$norm) - length(lcl$vectors$venn_lists[[1]]),
                                     length(lcl$vectors$venn_lists[[input$intersect_venn[2]]])),
                    "3" = venn_sample_3(nrepl = 100,
                                        intersectn = nrow(lcl$tables$venn_overlap),
                                        n = 1:ncol(mSet$dataSet$norm),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[1]]]),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[2]]]),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[3]]])),
                    "4" = venn_sample_4(nrepl = 100,
                                        intersectn = nrow(lcl$tables$venn_overlap),
                                        n = 1:ncol(mSet$dataSet$norm),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[1]]]),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[2]]]),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[3]]]),
                                        length(lcl$vectors$venn_lists[[input$intersect_venn[4]]])))
      stars = p2stars(pval)
      boxcolor = if(stars == "") "blue" else if(stars == "*") "green" else if(stars == "**") "yellow" else if(stars == "***") "orange" else if (stars == "****") "red"
      output$venn_pval <- shiny::renderUI({
        div(shiny::renderText({
          paste0("significance: ", stars, " (p = ", pval, ")" )
        }), style=paste0("background:", boxcolor))
      })
    })
  }else{
    pval = NULL
    output$venn_pval <- NULL
  }
  output$venn_tab <- DT::renderDataTable({
    # -------------
    DT::datatable(lcl$tables$venn_overlap,
                  selection = 'single',
                  rownames = F,
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
  })
})
