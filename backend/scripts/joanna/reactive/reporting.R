dataModal <- function(failed = FALSE) {
  modalDialog(
    fluidRow(align="center",
             textInput("report_plot_title", "Title", value = switch(input$statistics,
                                                                    tt = paste0("T-test", " - ", local$curr_mz, " m/z"),
                                                                    fc = paste0("Fold-change",  " - ", local$curr_mz, " m/z"),
                                                                    anova = paste0("ANOVA",  " - ", local$curr_mz, " m/z"),
                                                                    pca = "PCA",
                                                                    plsda = "PLS-DA",
                                                                    volc = "Volcano",
                                                                    meba = paste0("MEBA",  " - ", local$curr_mz, " m/z"),
                                                                    asca = paste0("ASCA",  " - ", local$curr_mz, " m/z"),
                                                                    ml = "Machine learning")),
             textAreaInput("report_plot_notes", "Notes", value = "", height = "100px")
    ),
    footer = tagList(
      fluidRow(align="center",
               modalButton("Cancel"),
               actionButton("report_plot", "Add plot to report")
      )

    )
  )
}

observeEvent(input$show_window, {
  showModal(dataModal())
  shinyjqui::jqui_draggable(selector = '.modal-content') # make draggable
})

# Report buttons
reportAppend = function(reportPlot, plotTitle, plotNotes){
  # if no report folder, create
  dir.create(file.path(options$work_dir, "report"))
  dir.create(file.path(options$work_dir, "report",
                       "figures"))

  # Report file name (Rmd)
  reportName <- file.path(options$work_dir, "report", paste("Report_", gsub("[^[:alnum:]]", "_", options$proj_name), ".Rmd", sep = ""))

  # If report file doesn't exist yet
  if (!file.exists(reportName)) {
    reportTmp <- file(reportName, open = "w")
    # Add Rmd base
    cat(paste(
      "---",
      paste("title:", options$proj_name),
      paste("author:", Sys.info()[["user"]]),
      paste("date:", Sys.Date()),
      "output: \n html_document: \n self_contained: no",
      "---", "\n", "```{r setup, include=FALSE}", "knitr::opts_chunk$set(echo = FALSE)", "```", "\n",
      sep = "\n"),
      file = reportTmp)
    close(reportTmp)
  }

  # append plot to file
  reportTmp <- file(reportName, open = "a")
  # unique file name to store figure for report
  tmpFigureName <- paste(tempfile(pattern = "plot", tmpdir = file.path(options$work_dir, "report", "figures")), ".png", sep = "")
  # save plo† as PDF
  # ggsave(tmpFigureName, plot = local$last_plot)
  ggsave(tmpFigureName, plot = reportPlot)
  dev.off()

  # htmlwidgets::saveWidget(local$last_plot, tmpFigureName)
  # Add code to Rmd file
  cat(paste(
    paste("#", plotTitle, sep = " "),

    # Input HTML in HTML (plotly interactive)
    # paste("<iframe width=\"600\" height=\"600\" src=\"", tmpFigureName, "\"></iframe>", sep = ""),

    # input figures in HTML
    paste("<img src = \"", tmpFigureName, "\">", sep = ""),

    # inputs PDF in HTML (zoom and scroll possible)
    # "---", "\n", "```{r, out.width = \"85%\"}", "\n", paste("knitr::include_graphics(\"", tmpFigureName, "\")", sep = ""), "```", "\n",
    plotNotes,
    "\n",
    sep = "\n"),
    file = reportTmp)
  close(reportTmp)

  # create HTML
  currentWD = getwd()
  setwd(file.path(options$work_dir, "report"))
  rmarkdown::render(reportName)
  setwd(currentWD)

  # reset title, caption, and notes
  updateTextInput(session, "plotTitle", label = "Title", value = "")
  updateTextAreaInput(session, "plotNotes", label = "Notes", value = "")
}
# observe report button presses, need one for each button*

observeEvent(input$report_plot, {
  reportAppend(local$last_plot, input$report_plot_title, input$report_plot_notes)
})

# For report tab where user can choose plots to appear in report
# nonselected


report_yes = reactiveValues(start = data.frame(),
                            now = data.frame())

report_no = reactiveValues(start = data.frame(c("a", "b", "c")),
                           now = data.frame(c("a", "b", "c")))

observe({
  if(exists("mSet")){
    if("storage" %in% names(mSet)){
      analyses = names(mSet$storage)
      report_no$start <- rbindlist(lapply(analyses, function(name){
        analysis = mSet$storage[[name]]$analysis
        analysis_names = names(analysis)
        # - - -
        with.subgroups <- intersect(analysis_names, c("ml", "plsr"))
        if(length(with.subgroups) > 0){
          extra_names <- lapply(with.subgroups, function(anal){
            switch(anal,
                   ml = {
                     which.mls <- intersect(c("rf", "ls"), names(analysis$ml))
                     ml.names = sapply(which.mls, function(meth){
                       if(length(analysis$ml[[meth]]) > 0){
                         paste0(meth, " - ", names(analysis$ml[[meth]]))
                       }
                     })
                     unlist(ml.names)
                   },
                   plsr = {
                     c ("plsda - PC1", "plsda - PC2", "plsda - PC3")
                   })
          })
          analysis_names <- c(setdiff(analysis_names, c("ml", "plsr", "plsda")), unlist(extra_names))
        }
        # - - -
        data.frame(
          paste0(analysis_names, " (", name, ")")
        )
      }))
      report_no$now <- report_no$start
    }else{
      report_no$start <- data.frame(names(mSet$analSet))
      report_no$now <- report_no$start
    }
  }
})

report_members <- reactiveValues(mzvals = list())

observeEvent(input$report_add, {
  # add to the 'selected' table
  rows <- input$report_unselected_rows_selected
  # get members and send to members list
  added = report_no$now[rows,]
  report_yes$now <- data.frame(c(unlist(report_yes$now), added))
  report_no$now <- data.frame(report_no$now[-rows,])
})

observeEvent(input$report_remove, {
  # add to the 'selected' table
  rows <- input$report_selected_rows_selected
  # get members and send to non members list
  removed = report_yes$now[rows,]
  report_no$now <- data.frame(c(unlist(report_no$now), removed))
  report_yes$now <- data.frame(report_yes$now[-rows,])
})


# the 'non-selected' table
output$report_unselected <- DT::renderDataTable({
  res = DT::datatable(data.table(), rownames=FALSE, colnames="excluded", options = list(dom = 'tp'))
  try({
    res = DT::datatable(report_no$now,rownames = FALSE, colnames="excluded", selection = "multiple", options = list(dom = 'tp'))
  })
  res
})

# the 'selected' table
output$report_selected <- DT::renderDataTable({
  res = DT::datatable(data.table(), rownames=FALSE, colnames="included", options = list(dom = 'tp'))
  try({
    res = DT::datatable(report_yes$now,rownames = FALSE, colnames="included", selection = "multiple", options = list(dom = 'tp'))
  })
  res
})
