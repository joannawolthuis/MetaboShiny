

shinyUI(fluidPage(theme = "metaboshiny.css",
                  navbarPage(title=h1("MassChecker"),
                             tabPanel(icon("battery-0"), 
                                      sidebarLayout(
                                        sidebarPanel(
                                          h1("File processing"),
                                          hr(),
                                          shinyDirButton("dir", "Chose directory", "Pick"),
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
                                          h1("Signal checking"),
                                          hr(),
                                          helpText("Generate plots"),
                                          actionButton("plot", "Plot"),
                                          hr(),
                                          helpText("Submit sorting results"),
                                          actionButton("submit", "Submit")
                                        ),
                                        mainPanel(
                                          wellPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 800px; background-color: #ffffff;",
                                            uiOutput("tics")
                                          )
                                        )
                                      )),
                             tabPanel(icon("battery-4"),
                                      helpText("lalala"))
                             ),
                  div(class="spinnylocation1",
                      div(class="plus", img(class="imagetop", src="program.svg", width="120px", height="120px")),
                      div(class="minus", img(class="imagebottom", src="program.svg", width="120px", height="120px"))),
                  div(class="line"))                    
        ) 
                       

