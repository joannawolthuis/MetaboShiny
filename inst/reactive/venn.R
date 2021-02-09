# nonselected

venn_yes = shiny::reactiveValues(start = data.frame(),
                                 now = data.frame())

venn_no = shiny::reactiveValues(
                        start = data.frame(c("a", "b", "c")),
                        now = data.frame(c("a", "b", "c")))

venn_members <- shiny::reactiveValues(mzvals = list())

shiny::observeEvent(input$venn_add, {
  # add to the 'selected' table
  rows <- input$venn_unselected_rows_selected
  # get members and send to members list
  added = venn_no$now[rows,]
  venn_yes$now <<- data.frame(c(unlist(venn_yes$now), added))
  venn_no$now <<- data.frame(venn_no$now[-rows,])
})

shiny::observeEvent(input$venn_remove, {
  # add to the 'selected' table
  rows <- input$venn_selected_rows_selected
  # get members and send to non members list
  removed = venn_yes$now[rows,]
  venn_no$now <<- data.frame(c(unlist(venn_no$now), removed))
  venn_yes$now <<- data.frame(venn_yes$now[-rows,])
})


# the 'non-selected' table
output$venn_unselected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames="excluded", options = list(dom = 'tp'))
  try({
    res = DT::datatable(venn_no$now, rownames = FALSE, colnames="excluded", selection = "multiple", options = list(dom = 'tp'))
  })
  res
})

# the 'selected' table
output$venn_selected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames="included", options = list(dom = 'tp'))
  try({
    res = DT::datatable(venn_yes$now,rownames = FALSE, colnames="included", selection = "multiple", options = list(dom = 'tp'))
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
    lcl$tables$venn_overlap <<- setdiff(uniqies,intersecties)

  }else{
    lcl$tables$venn_overlap <<- Reduce("intersect", lapply(input$intersect_venn, function(x){ # get the intersecting hits for the wanted tables
      lcl$vectors$venn_lists[[x]]
    }))
  }

  lcl$tables$venn_overlap <<- data.frame(mz = lcl$tables$venn_overlap)
  rownames(lcl$tables$venn_overlap) <<- lcl$tables$venn_overlap$mz

  # hypergeometric testing...
  if(length(input$intersect_venn) == 2 & nrow(lcl$tables$venn_overlap) > 0){

    pval = 1 - phyper(nrow(lcl$tables$venn_overlap),
                      length(lcl$vectors$venn_lists[[input$intersect_venn[1]]]),
                      ncol(mSet$dataSet$norm) - length(lcl$vectors$venn_lists[[1]]),
                      length(lcl$vectors$venn_lists[[input$intersect_venn[2]]]))
    stars = p2stars(pval)
    boxcolor = if(stars == "") "blue" else if(stars == "*") "green" else if(stars == "**") "yellow" else if(stars == "***") "orange" else if (stars == "****") "red"
    output$venn_pval <- shiny::renderUI({
      div(shiny::renderText({
        paste0("significance: ", stars, " (p = ", pval, ")" )
      }), style=paste0("background:", boxcolor))
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
