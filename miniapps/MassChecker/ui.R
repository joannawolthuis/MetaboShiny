shinyUI(fluidPage(theme = "metaboshiny.css",
                  navbarPage(title=h1("MassChecker"), tags$head(tags$script(src = "javascript.js")),
                             tabPanel(icon("battery-0"), 
                                      sidebarLayout(
                                        sidebarPanel(
                                          h1("File processing"),
                                          hr(),
                                          shinyDirButton("dir", "Choose directory", "Pick"),
                                          hr(),
                                          helpText("Convert .raw to .mzXML"),
                                          actionButton("convert", "Convert"),
                                          hr(),
                                          helpText("Check corrupt .mzXML"),
                                          actionButton("check", "Check"),
                                          hr()
                                        ),
                                        mainPanel(
                                          h4("Directory"),
                                          verbatimTextOutput("dir"), br(),
                                          h4("Files (.raw & .mzXML)"),
                                          verbatimTextOutput("files")
                                        )
                                      )),
                             tabPanel(icon("battery-1"),
                                      sidebarLayout(
                                        sidebarPanel(
                                          wellPanel(id = "lPanel",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                                    
                                          h1("Signal checking"),
                                          hr(),
                                          helpText("Generate plots"),
                                          actionButton("plot", "Plot"),
                                          hr(),
                                          helpText("Submit sorting results"),
                                          actionButton("submit", "Submit"),
                                          shinySaveButton("ticreport_save_loc", "Save report", title = "Choose a file name",filetype = ".csv"),
                                          hr(),
                                          div(DT::dataTableOutput('ticTable'),style='font-size:80%')
                                        )),
                                        mainPanel(
                                          wellPanel(id = "rPanel",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                            uiOutput("tics")
                                          )
                                        )
                                      )),
                             # tabPanel(icon("battery-2"),
                             #          sidebarLayout(
                             #            sidebarPanel(
                             #              wellPanel(id = "lPanel2",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                             #                        textInput("report_title", label = "Formal project name"),
                             #                        hr(),
                             #                        shinyFilesButton("sample_sheet_loc",label = "Upload sample sheet", title="Open excel sheet", multiple = F),
                             #                        verbatimTextOutput("sample_sheet"),
                             #                        hr(),
                             #                        actionButton("make_report", "Generate report"),
                             #                        hr(),
                             #                        shinySaveButton("report_save_loc", "Save report", title = "Choose a file name",filetype = ".csv")
                             #              )),
                             #            mainPanel(
                             #              wellPanel(id = "rPanel2",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                             #                        div(DT::dataTableOutput('final_report'),style='font-size:80%')
                             #              )
                             #            )
                             #          )),
                             tabPanel(icon("battery-2"),
                                      sidebarLayout(
                                        sidebarPanel(
                                          wellPanel(id = "lPanel2", style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                                    fileInput("inj_loc", "Upload injection order"),
                                                    actionButton("create_sample_names", "Create sample sheet")
                                        )),
                                        mainPanel(
                                          wellPanel(id = "rPanel2",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                                    dataTableOutput("sample_names")
                                          )
                                        )
                                      )),
                             tabPanel(icon("battery-3"),
                                      sidebarLayout(
                                        sidebarPanel(
                                          wellPanel(id = "lPanel2",style = "#overflow-y:scroll; max-height: 1100px; background-color: #ffffff;",
                                                    fluidRow(column(12,
                                                                    align="center",
                                                                    icon("magic","fa-2x"),
                                                                    h4("Pipeline settings"),
                                                                    sliderInput("nrepl", "Replicates", 0, 10, 3),
                                                                    hr(),
                                                                    sliderInput("thresh_pos", "Positive threshold", 1000, 10000, 2000),
                                                                    sliderInput("thresh_neg", "Negative threshold", 1000, 10000, 2000),
                                                                    sliderInput("dimsThresh", "DIMS threshold", 0, 2000, 100),
                                                                    sliderInput("trim", "Trim", 0, 0.5, 0.1),
                                                                    sliderInput("resol", "Resolution", 100000, 200000, 140000),
                                                                    sliderInput("ppm", "Error margin", 1, 3, 2, step=0.5, post = " ppm"),
                                                                    sliderInput("cores", 
                                                                                "Cores used", 
                                                                                1, 
                                                                                parallel::detectCores(), 
                                                                                1,
                                                                                step = 1),
                                                                    tags$div(HTML('<div id="peak_calling" class="form-group shiny-input-radiogroup shiny-input-container shiny-input-container-inline">
                                                                                  <label class="control-label" for="peak_calling">Peak calling style</label>
                                                                                  <div class="shiny-options-group">
                                                                                  <label class="radio-inline">
                                                                                  <input type="radio" name="peak_calling" value="gaussian" checked="checked"/>
                                                                                  <span><img src="gaussian.jpg" height="42" width="42"/></span>
                                                                                  </label>
                                                                                  <label class="radio-inline">
                                                                                  <input type="radio" name="peak_calling" value="wavelet"/>
                                                                                  <span><img src="wavelet.png" height="42" width="42"/></span>
                                                                                  </label>
                                                                                  </div>
                                                                                  </div>')),
                                                                    radioButtons("peak_grouping",label = "Peak grouping style",choices = list("Mean shifting"="meanshift",
                                                                                                                                              "Hierarchical clustering"="hclust"),inline = T, 
                                                                                 selected="meanshift"),
                                                                    actionButton('start_pipeline', 'Go', icon=icon("paw"),style='font-size:150%')
                                                                    )
                                                             )
                                          )),
                                        mainPanel(
                                          wellPanel(id = "rPanel2",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                                    bsCollapse(id = "pipeline",
                                                               bsCollapsePanel("Initial setup",
                                                                               actionButton("do_step_1", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("DIMS data extraction", 
                                                                               actionButton("do_step_2", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("Average technical replicates", 
                                                                               actionButton("do_step_3", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("Peak finding", 
                                                                               actionButton("do_step_4", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("Collect samples (I)", 
                                                                               actionButton("do_step_5", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("Peak grouping", 
                                                                               actionButton("do_step_6", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("Fill missing values", 
                                                                               actionButton("do_step_7", 'Retry'),
                                                                               style = "warning"),
                                                               bsCollapsePanel("Collect samples (II)", 
                                                                               actionButton("do_step_8", 'Retry'),
                                                                               style = "warning")
                                                    )
                                          )
                                          )
                                        )
                                      )
                             ),
                  div(class="spinnylocation1",
                      div(class="plus", img(class="imagetop", src="program.svg", width="120px", height="120px")),
                      div(class="minus", img(class="imagebottom", src="program.svg", width="120px", height="120px"))),
                  div(class="line"))                    
        ) 
                       
