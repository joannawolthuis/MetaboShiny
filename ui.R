# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).
library(datasets)

asca.table <- read.csv("backend/appdata/asca_sigAB.csv")
mode <- "time"
# Use a fluid Bootstrap layout

navbarPage("MetaboShiny", id="nav_general",
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Database check", 
                    helpText("placeholder")),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Experiment import", 
                    helpText("placeholder")),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Match finder", 
                    helpText("placeholder")),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Analysis", # multiple options here :-)
                    if(mode == "time"){
                                        navbarPage("Time Series", id="nav_time",
                                                   tabPanel("iPCA", helpText("placeholder")),
                                                   # =================================================================================
                                                   tabPanel("Heatmap", helpText("placeholder")),
                                                   # =================================================================================
                                                   tabPanel("MEBA", helpText("placeholder")),
                                                   # =================================================================================
                                                   tabPanel("ASCA",
                                                            navbarPage("", 
                                                              tabPanel("Plots & ID",
                                                                       
                                                            fluidRow(column(4,div(DT::dataTableOutput('asca_tab'),style='font-size:80%')),
                                                                     column(8, plotOutput('asca_plot'))
                                                                      ),
                                                            fluidRow(
                                                              column(2, checkboxGroupInput("checkGroup", 
                                                                                           label = h4("Search mz in..."), 
                                                                                           choices = list("Internal DB" = "matches_internal", 
                                                                                                          "Noise DB" = "matches_noise", 
                                                                                                          "HMDB" = "matches_hmdb",
                                                                                                          "ChEBI" = "matches_chebi"),
                                                                                           selected = 1),actionButton("do", "Update Table")),
                                                              column(3, verbatimTextOutput("sql_results")),
                                                              column(4,div(DT::dataTableOutput('match_tab'),style='font-size:80%'))
                                                                    )
                                                              )
                                                          )
                                                          )
                                                   # =================================================================================
                                                   )
                                      }
           # --------------------------------------------------------------------------------------------------------------------------------
))