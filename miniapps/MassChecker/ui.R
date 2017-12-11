shinyUI(fluidPage(theme = "metaboshiny.css",
                  navbarPage(title=h1("MassChecker"),
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
                             tabPanel(icon("battery-2"),
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
                             tabPanel(icon("battery-4"),
                                      sidebarLayout(
                                        sidebarPanel(
                                          wellPanel(id = "lPanel2",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                                    textInput("report_title", label = "Formal project name"),
                                                    hr(),
                                                    shinyFilesButton("sample_sheet_loc",label = "Upload sample sheet", title="Open excel sheet", multiple = F),
                                                    verbatimTextOutput("sample_sheet"),
                                                    hr(),
                                                    actionButton("make_report", "Generate report"),
                                                    hr(),
                                                    shinySaveButton("report_save_loc", "Save report", title = "Choose a file name",filetype = ".csv")
                                          )),
                                        mainPanel(
                                          wellPanel(id = "rPanel2",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                                    div(DT::dataTableOutput('final_report'),style='font-size:80%')
                                          )
                                        )
                                      ))
                             ),
                  div(class="spinnylocation1",
                      div(class="plus", img(class="imagetop", src="program.svg", width="120px", height="120px")),
                      div(class="minus", img(class="imagebottom", src="program.svg", width="120px", height="120px"))),
                  div(class="line"))                    
        ) 
                       
