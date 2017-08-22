# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).

navbarPage("MetaboShiny", id="nav_general",windowTitle = "MetaboShiny",
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Databases" , icon = icon("database"),
                    # -- header row ---
                    fluidRow(column(11, align="center",
                                    h3("")
                             ),
                    # --- db check cols ---
                    fluidRow(column(4, align="center",
                                    h2("UMC"),
                                    helpText("Internal commonly known metabolites and common pollutants found in DIMS."),
                                    br(),
                                    imageOutput("umc_logo",inline = T),
                                    br(),br(),br(),
                                    actionButton("check_umc", "Check", icon = icon("check")),
                                    actionButton("build_umc", "Build", icon = icon("wrench")),
                                    br(),br(),
                                    imageOutput("umc_check",inline = T)
                                    ),
                             column(3,  align="center",
                                    h2("HMDB"),
                                    helpText("Metabolites commonly found in human biological samples."),
                                    br(),
                                    imageOutput("hmdb_logo",inline = T),
                                    br(),br(),br(),
                                    actionButton("check_hmdb", "Check", icon = icon("check")),
                                    actionButton("build_hmdb", "Build", icon = icon("wrench")),
                                    br(),br(),
                                    imageOutput("hmdb_check",inline = T)
                                    ),
                             column(4,  align="center",
                                    h2("ChEBI"),
                                    helpText("A broad database with known chemicals of biological interest."),
                                    br(),
                                    imageOutput("chebi_logo",inline = T),
                                    br(),br(),br(),
                                    actionButton("check_chebi", "Check", icon = icon("check")),
                                    actionButton("build_chebi", "Build", icon = icon("wrench")),
                                    br(),br(),
                                    imageOutput("chebi_check",inline = T)
                                    )
                             )
                    )
                    ),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Import", icon = icon("upload"), 
                    br(),br(),br(),
                    fluidRow(column(3,  align="center",
                                    imageOutput("pos_icon",inline = T),
                                    br(),br(),
                                    fileInput("outlist_pos", "Positive outlist",
                                                  multiple = F,
                                                  accept = c(".RData"))),
                             column(3,  align="center",
                                    imageOutput("excel_icon",inline = T),
                                    br(),br(),
                                    fileInput("excel", "Excel workbook",
                                              multiple = F,
                                              accept = c(".xslx"))),
                             column(3,  align="center",
                                    imageOutput("neg_icon",inline = T),
                                    br(),br(),
                                    fileInput("outlist_neg", "Negative outlist",
                                                  multiple = F,
                                                  accept = c(".RData")))
                             ),
                    fluidRow(column(9, align="center",
                                    actionButton("create_db", "Create database", icon = icon("magic")),
                                    br(),br(),hr(),
                                    imageOutput("db_icon",inline = T),
                                    br(),br(),
                                    fileInput("pat_db", "Import previous experiment",
                                              multiple = F,
                                              accept = c(".db")),
                                    actionButton("import_db", "Import", icon = icon("hand-peace-o")),
                                    imageOutput("db_upload_check",inline = T)
                                    ))
                    ),
           tabPanel("Document", icon=icon("file-text-o"),
                    fluidRow(column(4, align="center",
                                    selectInput('exp_var', 'Which experimental variable do you want to look at?', choices = c("")),
                                    actionButton("check_excel", "Get variables", icon=icon("refresh")),
                                    br(),br(),
                                    selectInput('exp_type', 'What type of analysis do you want to do?', choices = c("Time series", 
                                                                                                                    "Bivariate comparison",
                                                                                                                    "Multivariate comparison")),
                                    actionButton("create_csv", "Create CSV", icon=icon("file-text-o")),
                                    br(),hr(),
                                    imageOutput("csv_icon",inline = T),
                                    br(),br(),
                                    fileInput("pat_csv", "Import CSV",
                                              multiple = F,
                                              accept = c(".csv")),
                                    actionButton("import_csv", "Import", icon = icon("hand-peace-o")),
                                    imageOutput("csv_upload_check",inline = T)
                                    ),
                             column(7, 
                                    fluidRow(div(DT::dataTableOutput('csv_tab'),style='font-size:80%'))
                                    )
                             )
                    ),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Clean",  icon = icon("filter"),
                    column(3, aligh="center",
                           selectInput('filt_type', 'How will you filter your m/z values?', choices = c("Interquantile range")),
                           selectInput('norm_type', 'What type of normalization do you want to do?', choices = c("Quantile normalization")),
                           selectInput('trans_type', 'How will you transform your data?', choices = c("Log transform")),
                           selectInput('scale_type', 'How will you scale your data?', choices = c("Autoscale")),
                           actionButton("initialize", "Go", icon=icon("hand-o-right"))
                           ), column(9,
                                     navbarPage("Explore",
                                                tabPanel("Variables", icon=icon("braille"),
                                                         plotOutput('var_norm_plot', width = '800px', height='800px')),
                                                tabPanel("Samples", icon=icon("tint"),
                                                         plotOutput('samp_norm_plot', width = '800px', height='800px'))
                                                ))

           ),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Analysis",  icon = icon("bar-chart"),# multiple options here :-)
                    column(8, 
                    if(mode == "time"){
                                        navbarPage("Time Series", id="nav_time",
                                                   tabPanel("iPCA", 
                                                            plotlyOutput("plot_ipca" ),
                                                            selectInput("ipca_factor", label = "Color based on:", choices =list("Time"="facA",
                                                                                                                                "Experimental group"="facB"))
                                                            ),
                                                   # =================================================================================
                                                   #tabPanel("Heatmap",
                                                   #         plotOutput("time_heat_plot", height='600px', width='600px')
                                                   #        ),
                                                   # =================================================================================
                                                   tabPanel("MEBA", 
                                                            fluidRow(plotOutput('meba_plot')),
                                                            fluidRow(div(DT::dataTableOutput('meba_tab'),style='font-size:80%'))),
                                                   # =================================================================================
                                                   tabPanel("ASCA",
                                                            navbarPage("Explore", 
                                                              tabPanel("Overview", icon=icon("eye"), helpText("...")),
                                                              tabPanel("Plots", icon=icon("bar-chart-o"),
                                                                     fluidRow(plotOutput('asca_plot')),
                                                                     fluidRow(div(DT::dataTableOutput('asca_tab'),style='font-size:80%'))
                                                              )
                                                          )
                                                          )
                                                   # =================================================================================
                                                   )
                                      }), column(4, align="center",
                                                 imageOutput("find_mol_icon",inline = T),
                                                 div(checkboxGroupInput("checkGroup", 
                                                                              label = h4("Find selected m/z in:"), 
                                                                              choices = list("Internal DB" = "matches_internal", 
                                                                                              "Noise DB" = "matches_noise", 
                                                                                              "HMDB" = "matches_hmdb",
                                                                                              "ChEBI" = "matches_chebi"),
                                                                              selected = 1,inline = T),style='font-size:90%'),
                                                     actionButton("search_mz", "Search", icon=icon("search")),
                                                 hr(),
                                                   div(DT::dataTableOutput('match_tab'),style='font-size:80%'),
                                                 hr(),
                                                   div(textOutput("curr_definition"))
                                                 )
           # --------------------------------------------------------------------------------------------------------------------------------
           ),tabPanel("",  icon = icon("cog"), 
                    h2("General settings for this application"),
                    hr(),
                    shinyDirButton("get_work_dir", "Choose a working directory" ,
                                   title = "Browse",
                                   buttonType = "default", class = NULL),
                    helpText("Your results will be stored here for later access."),
                    textOutput("exp_dir"),
                    hr(),
                    textInput(inputId="proj_name", label="Project name", value = ''),
                    actionButton("set_proj_name", label="Apply"),
                    helpText("This name will be used in all save files."),
                    textOutput("proj_name")
           )
           )