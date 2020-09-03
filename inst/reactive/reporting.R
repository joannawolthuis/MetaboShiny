shiny::observeEvent(input$star_mz, {
    if(!is.null(mSet)){
      if(!("report" %in% names(mSet))){
        mSet$report <<- list(mzStarred = data.table::data.table(mz = colnames(mSet$dataSet$norm),
                                                               star = c(FALSE)))  
        data.table::setkey(mSet$report$mzStarred, mz)
      }
    }
    if(my_selection$mz != ""){
      
      mSet$report$mzStarred[my_selection$mz]$star <<- input$star_mz
      
      try({
        tablemanager$make <- plotmanager$make <- input$statistics
      })
    }
})

# nonselected

report_yes = shiny::reactiveValues(start = data.frame(),
                                 now = data.frame())

report_no = shiny::reactiveValues(
  start = data.frame(c("a", "b", "c")),
  now = data.frame(c("a", "b", "c")))

report_members <- shiny::reactiveValues(mzvals = list())

shiny::observeEvent(input$report_add, {
  # add to the 'selected' table
  rows <- input$report_unselected_rows_selected
  # get members and send to members list
  added = report_no$now[rows,]
  report_yes$now <- data.frame(c(unlist(report_yes$now), added))
  report_no$now <- data.frame(report_no$now[-rows,])
})

shiny::observeEvent(input$report_remove, {
  # add to the 'selected' table
  rows <- input$report_selected_rows_selected
  # get members and send to non members list
  removed = report_yes$now[rows,]
  report_no$now <- data.frame(c(unlist(report_no$now), removed))
  report_yes$now <- data.frame(report_yes$now[-rows,])
})

# the 'non-selected' table
output$report_unselected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames="Not in report", options = list(dom = 'tp'))
  try({
    res = DT::datatable(report_no$now, rownames = FALSE, colnames="Not in report", selection = "multiple", options = list(dom = 'tp'))
  })
  res
})

# the 'selected' table
output$report_selected <- DT::renderDataTable({
  res = DT::datatable(data.table::data.table(), rownames=FALSE, colnames="In report", options = list(dom = 'tp'))
  try({
    res = DT::datatable(report_yes$now,rownames = FALSE, colnames="In report", selection = "multiple", options = list(dom = 'tp'))
  })
  res
})

# if(interactive()){
#   # GENERAL IDEA
#   "
# PREP
# create mSet$report list (mReport)
# example:
# list(analyses=list(tt = list(
#                    plots = list(),
#                    tables = list()),
#                    mz = c(list(mz = .., matchRow = ..., note = ..., db = ...)))
# 
# RUN ANALYSES OF INTEREST 
# example: t-test
# STAR specific plot or tables
# ADD to mReport
# 
# SELECT MZ AND DESC OF INTEREST
# sift through mz of interest
# if interesting, STAR that m/z 
# OR star interesting description (savd)
# 
# create template for report
# 
# CREATING REPORT
# WRITE HEADER
# - logo - 
# - title - 
# - date - 
# - user experiment description -
# LOOP through analysis types
# FOR EACH EXPERIMENT
#   - ROW 1:(overview plot(s)) (starred mz highlighted)
#   - ROW 2:
#     - table with starred m/zs 
#     - OR top x results
#   - ROW 3: 
#     - top x summary figures
#     - OR FOR EACH m/z starred:
#       match table rows (name, adduct, iso, description (only if starred))
# 
# "
#   
#   # SERVER
#   
#   ## rmd example ##
#   # load data 
#   set.seed(500)
#   Score <- rnorm(40, 100, 15)
#   Criteria1<-rnorm(40, 10, 5)
#   Criteria2<-rnorm(40, 20, 5)
#   ID <- sample(1:1000,8,replace=T)
#   df <- data.frame(ID,Score,Criteria1,Criteria2)
#   
#   library("rmarkdown")
#   
#   experimentTitle = "myExperiment"
#   fullpath = system.file("www/gemmy_rainbow.png", package="MetaboShiny")
#   template = "
# ```{r, echo=FALSE, fig.align='center',out.width = '150px'}
# knitr::include_graphics(fullpath)
# ```      
# <center> <h1>`r experimentTitle`</h1> </center>
# <center> <h3>analysed with MetaboShiny</h3> </center>
# 
# Test OWO haha.
# 
# ```{r, echo=FALSE}
# #Report Analysis
# summary(subgroup)
# ```"
#   tempLoc = "hewwo.Rmd"
#   writeLines(template, tempLoc)
#   
#   id = unique(df$ID)[1]
#   subgroup <- df[df$ID == id,]
#   
#   render(input = tempLoc,
#          output_dir = getwd(),
#          output_file = paste0('report.', id, '.html'))    
#   
# }