# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).

navbarPage("MetaboShiny", id="nav_general",
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
                             ),
                    fluidRow(column(11, align="center",
                                    br(),br(),hr(),
                                    imageOutput("db_icon2",inline = T),
                                    br(),br(),
                                    fileInput("own_db", "Build your own",
                                              multiple = F,
                                              accept = c(".csv")),
                                    actionButton("import_files", "Import", icon = icon("hand-peace-o"))
                    ))
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
                                    actionButton("import_files", "Create database", icon = icon("magic")),
                                    br(),br(),hr(),
                                    imageOutput("db_icon",inline = T),
                                    br(),br(),
                                    fileInput("pat_db", "Import experiment database",
                                              multiple = F,
                                              accept = c(".db")),
                                    actionButton("import_files", "Import", icon = icon("hand-peace-o"))
                                    ))
                    ),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("Match", icon = icon("search"), 
                    fluidRow(
                      column(3, align="center",
                             checkboxGroupInput("match_group",
                                                label = h4("Match m/z values with..."),
                                                choices = list("Internal DB" = "Internal", 
                                                               "Noise DB" = "Noise",
                                                               "HMDB" = "HMDB",
                                                               "ChEBI" = "ChEBI"),
                                                selected = 1),
                             actionButton("match_start", "Let's match!", icon = icon("search"))),
                      column(7, helpText("statistics on matches here"))
                      )
           ),
           tabPanel("Setup",  icon = icon("sliders"),
                    column(3, aligh="center",
                           helpText("Use this window to setup your experiment."),
                           selectInput('exp_var', 'Which experimental variable do you want to look at?', choices = c("")),
                           actionButton("check_excel", "Get variables", icon=icon("refresh")),
                           br(),br(),
                           selectInput('exp_type', 'What type of analysis do you want to do?', choices = c("Time series", 
                                                                                                           "Bivariate comparison",
                                                                                                           "Multivariate comparison")), 
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
                    column(7, 
                    if(mode == "time"){
                                        navbarPage("Time Series", id="nav_time",
                                                   tabPanel("iPCA", helpText("placeholder")),
                                                   # =================================================================================
                                                   tabPanel("Heatmap", helpText("placeholder")),
                                                   # =================================================================================
                                                   tabPanel("MEBA", helpText("placeholder")),
                                                   # =================================================================================
                                                   tabPanel("ASCA",
                                                            navbarPage("Explore", 
                                                              tabPanel("Overview", icon=icon("eye"), helpText("...")),
                                                              tabPanel("Plots", icon=icon("bar-chart-o"),
                                                                     fluidRow(plotOutput('asca_plot',height = '300px', width='700px')),
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
                                                 hr()
                                                 )
           # --------------------------------------------------------------------------------------------------------------------------------
           ),
           tabPanel("",  icon = icon("cog"), 
                    helpText("placeholder"))# multiple options here :-))
                    
           )