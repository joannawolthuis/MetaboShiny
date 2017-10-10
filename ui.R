# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).

source("./Rsource/SwitchButton.R")
sardine <- function(content) div(style="display: inline-block;vertical-align:top;", content)

shinyUI(fluidPage(
  tags$head(
    tags$style(HTML("
                    @import url('//fonts.googleapis.com/css?family=Bungee Outline');
                    
                    h1 {
                    margin: 2px;
                    font-family: 'Bungee Outline';
                    font-weight: 15;
                    line-height: 0.5;
                    }
                    
                    "))
    ),theme = "button.css",
                  if(options$packages_installed != "Y"){
                    navbarPage(title=h1("MetaboSetup"), 
                               id="nav_setup",
                               windowTitle = "MetaboShiny",
                               tabPanel("", icon = icon("wrench"), value="setup",
                                        # --- db check cols ---
                                        fluidRow(column(width=2),column(width=5, align="center",
                                                                        h2("Setup"),
                                                                        br(),
                                                                        imageOutput("cute_package",inline = T),
                                                                        hr(),
                                                                        helpText("Needed packages:"),
                                                                        div(DT::dataTableOutput('package_tab'),style='font-size:80%'),
                                                                        br(),                                  
                                                                        actionButton("install_packages", "Install", icon = icon("wrench")),
                                                                        br(),br(),
                                                                        imageOutput("package_check")
                                        )))
                    )  
                  } else{
                    navbarPage(title=h1("MetaboShiny"),
                               id="nav_general",
                               windowTitle = "MetaboShiny",
                               tabPanel("", icon = icon("share-alt"), value="setup",
                                        # --- db check cols ---
                                        fluidRow(column(width=2),column(width=5, align="center",
                                                                        h2("Setup"),
                                                                        br(),
                                                                        imageOutput("cute_package",inline = T),
                                                                        hr(),
                                                                        helpText("Installed packages:"),
                                                                        div(DT::dataTableOutput('package_tab'),style='font-size:80%'),
                                                                        br(),                                  
                                                                        #                                        actionButton("install_packages", "Install", icon = icon("wrench")),
                                                                        actionButton("update_packages", "Update", icon = icon("star")),
                                                                        br(),br(),
                                                                        imageOutput("package_check")
                                        ))),
                               tabPanel("",  icon = icon("save"), value="save",
                                        helpText("Save your dataset")
                                        #fluidRow(column(width=1,fadeImageButton("test", img.path = "cutemolecule.png")),
                                        #         column(width=1,fadeImageButton("test2", img.path = "cutemolecule.png"))
                                        #         )
                                        ),
                               # --------------------------------------------------------------------------------------------------------------------------------
                               tabPanel("", icon = icon("database"), value="database",
                                        # -- header row ---
                                        fluidRow(column(12, align="center",
                                                        h3("")
                                        ),
                                        # --- db check cols --- common pollutants found in DIMS.
                                        fluidRow(column(3, align="center",
                                                        h2("UMC Internal"),
                                                        helpText("Internal commonly known metabolites."),
                                                        br(),
                                                        imageOutput("umc_logo_int",inline = T),
                                                        br(),br(),br(),
                                                        actionButton("check_internal", "Check", icon = icon("check")),
                                                        actionButton("build_internal", "Build", icon = icon("wrench")),
                                                        br(),br(),
                                                        imageOutput("internal_check",inline = T)
                                        ),column(3, align="center",
                                                 h2("UMC Noise"),
                                                 helpText("Internal common pollutants found in DIMS using local method."),
                                                 br(),
                                                 imageOutput("umc_logo_noise",inline = T),
                                                 br(),br(),br(),
                                                 actionButton("check_noise", "Check", icon = icon("check")),
                                                 actionButton("build_noise", "Build", icon = icon("wrench")),
                                                 br(),br(),
                                                 imageOutput("noise_check",inline = T)
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
                                        )),
                                        fluidRow(column(3,  align="center",
                                                        h2("ChEBI"),
                                                        helpText("A broad database with known chemicals of biological interest."),
                                                        br(),
                                                        imageOutput("chebi_logo",inline = T),
                                                        br(),br(),br(),
                                                        actionButton("check_chebi", "Check", icon = icon("check")),
                                                        actionButton("build_chebi", "Build", icon = icon("wrench")),
                                                        br(),br(),
                                                        imageOutput("chebi_check",inline = T)
                                        ),column(3,  align="center",
                                                        h2("WikiPathways"),
                                                        helpText("Compounds associated with open source biological pathways."),
                                                        imageOutput("wikipath_logo",inline = T),
                                                        br(),br(),
                                                        actionButton("check_wikipathways", "Check", icon = icon("check")),
                                                        actionButton("build_wikipathways", "Build", icon = icon("wrench")),
                                                        br(),br(),
                                                        imageOutput("wikipathways_check",inline = T)
                                        ),column(3,  align="center",
                                                        h2("KEGG"),
                                                        helpText("Compounds associated with biological pathways."),
                                                        imageOutput("kegg_logo",inline = T),
                                                        br(),br(),
                                                        actionButton("check_kegg", "Check", icon = icon("check")),
                                                        actionButton("build_kegg", "Build", icon = icon("wrench")),
                                                        br(),br(),
                                                        imageOutput("kegg_check",inline = T)
                                                        )
                                        # ,
                                        #          column(3,  align="center",
                                        #                 h2("PubChem"),
                                        #                 helpText("A huge database with known pretty much all known chemicals."),
                                        #                 br(),
                                        #                 imageOutput("pubchem_logo",inline = T),
                                        #                 br(),br(),br(),
                                        #                 actionButton("check_pubchem", "Check", icon = icon("check")),
                                        #                 actionButton("build_pubchem", "Build", icon = icon("wrench")),
                                        #                 br(),br(),
                                        #                 imageOutput("pubchem_check",inline = T)
                                        #                 )
                                        )
                                        )
                               ),
                               # --------------------------------------------------------------------------------------------------------------------------------
                               tabPanel("", icon = icon("upload"), value="upload", ## this guy gives error???
                                        fluidRow(column(9, align="center", h4("Create project"), br())),
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
                                                        actionButton("create_db", "Go", icon = icon("magic")),
                                                        br(),br(),hr(),
                                                        h4("Import project"),
                                                        br(),
                                                        imageOutput("db_icon",inline = T),
                                                        br(),br(),
                                                        fileInput("pat_db", "",
                                                                  multiple = F,
                                                                  accept = c(".db")),
                                                        actionButton("import_db", "Import", icon = icon("hand-peace-o")),
                                                        imageOutput("db_upload_check",inline = T)
                                        ))
                               ),
                               tabPanel("", icon=icon("file-text-o"), value="document",
                                        fluidRow(column(4, align="center",
                                                        selectInput('exp_var', 'Which experimental variable do you want to look at?', choices = c("")),
                                                        actionButton("check_excel", "Get variables", icon=icon("refresh")),
                                                        br(),br(),
                                                        selectInput('exp_type', 'What type of analysis do you want to do?', choices = list("Standard" = "stat.standard",
                                                                                                                                           "Time series - standard" = "time.standard", 
                                                                                                                                           "Time series - custom end" = "time.custom",
                                                                                                                                           "Time series - subtract" = "time.subtract"
                                                        )),
                                                        uiOutput("exp_opt"),
                                                        br(),br(),
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
                               tabPanel("",  icon = icon("shower"), value="filter",
                                        fluidRow(column(3, aligh="center",
                                               selectInput('filt_type', 'How will you filter your m/z values?', choices = list("Interquantile range"="iqr",
                                                                                                                               "Relative stdev"="rsd",
                                                                                                                               "Non-parametric relative stdev"="nrsd",
                                                                                                                               "Mean"="mean",
                                                                                                                               "Standard deviation"="sd",
                                                                                                                               "Median absolute deviation"="mad",
                                                                                                                               "Median"="median")),
                                               selectInput('norm_type', 'What type of normalization do you want to do?', choices = list("Quantile normalization"="QuantileNorm",
                                                                                                                                        "By reference feature"="CompNorm",
                                                                                                                                        "Sum"="SumNorm",
                                                                                                                                        "Median"="MedianNorm")),
                                               uiOutput("ref_select"),
                                               selectInput('trans_type', 'How will you transform your data?', choices = list("Log transform"="LogNorm",
                                                                                                                             "Cubic root transform"="CrNorm",
                                                                                                                             "None"="N/A")),
                                               selectInput('scale_type', 'How will you scale your data?', choices = list("Autoscale"="AutoNorm",
                                                                                                                         "Mean-center"="MeanCenter",
                                                                                                                         "Pareto Scaling"="ParetoNorm",
                                                                                                                         "Range scaling"="RangeNorm",
                                                                                                                         "None"="N/A")),
                                               actionButton("initialize", "Go", icon=icon("hand-o-right")),
                                               hr(),
                                               imageOutput("dataset_icon",inline = T),
                                               fileInput("pat_dataset", "Import dataset",
                                                         multiple = F,
                                                         accept = c(".RData")),
                                               actionButton("import_dataset", "Import", icon = icon("hand-peace-o")),
                                               imageOutput("dataset_upload_check",inline = T)
                                        ), column(9,
                                                  navbarPage("Explore",
                                                             tabPanel("Variables", icon=icon("braille"),
                                                                      plotOutput('var_norm_plot', width = '800px', height='800px')),
                                                             tabPanel("Samples", icon=icon("tint"),
                                                                      plotOutput('samp_norm_plot', width = '800px', height='800px'))
                                                  ))
                                        )),
                               # --------------------------------------------------------------------------------------------------------------------------------
                               tabPanel("",  icon = icon("bar-chart"), value = "analysis",# multiple options here :-)
                                        sidebarLayout(position="right",
                                                      mainPanel = mainPanel(
                                                        uiOutput("analUI")
                                                      ),
                                                      sidebarPanel = sidebarPanel(align="center",
                                                        fluidRow(imageOutput("find_mol_icon",inline = T),
                                                                 h4("Selected database(s):"),br(),
                                                                 sardine(fadeImageButton("search_internal", img.path = "umcinternal.png")),
                                                                 sardine(fadeImageButton("search_noise", img.path = "umcnoise.png")),
                                                                 sardine(fadeImageButton("search_hmdb", img.path = "hmdblogo.png")),br(),
                                                                 sardine(fadeImageButton("search_chebi", img.path = "chebilogo.png")),
                                                                 sardine(fadeImageButton("search_wikipathways", img.path = "wikipathways.png")),
                                                                 sardine(fadeImageButton("search_kegg", img.path = "kegglogo.gif")),
                                                                 #sardine(fadeImageButton("search_pubchem", img.path = "pubchemlogo.png")),
                                                                 br(),
                                                                 sardine(switchButton(inputId = "autosearch",
                                                                              label = "Autosearch", 
                                                                              value = FALSE, col = "BW", type = "OO"))
                                                                 ),
                                                        hr(),
                                                        h4("Current compound:"),
                                                        verbatimTextOutput("curr_mz"), hr(),
                                                        fluidRow(navbarPage("Search", id="tab_iden",
                                                                            tabPanel("Current", icon=icon("sort-numeric-asc"),
                                                                                     actionButton("search_mz", "Find hits", icon=icon("search")),
                                                                                     hr(),
                                                                                     div(DT::dataTableOutput('match_tab'),style='font-size:80%'),
                                                                                     hr(),
                                                                                     div(textOutput("curr_definition"))
                                                                            ),
                                                                            tabPanel("Target", icon=icon("sort-alpha-asc"),
                                                                                     actionButton("browse_db", "Browse compounds", icon=icon("eye")),
                                                                                     hr(),
                                                                                     div(DT::dataTableOutput('browse_tab'),style='font-size:80%'),
                                                                                     hr(),
                                                                                     div(textOutput("browse_definition"),style='font-size:80%'),
                                                                                     hr(),
                                                                                     actionButton("search_cpd", "Find hits", icon=icon("search")),
                                                                                     hr(),
                                                                                     div(DT::dataTableOutput('hits_tab'),style='font-size:80%')
                                                                            )))
                                                      )
                                                      )),
                               tabPanel("",  icon = icon("cog"), value="options", 
                                        navbarPage("Settings", id="tab_settings",
                                                   tabPanel("Project", icon=icon("gift"),
                                                            textInput(inputId="proj_name", label="Project name", value = ''),
                                                            actionButton("set_proj_name", label="Apply"),
                                                            helpText("This name will be used in all save files."),
                                                            textOutput("proj_name")
                                                   ),
                                                   tabPanel("Storage", icon=icon("folder-open-o"),
                                                            shinyDirButton("get_db_dir", "Choose a database directory" ,
                                                                           title = "Browse",
                                                                           buttonType = "default", class = NULL),
                                                            helpText("Your databases will be stored here. 500GB recommended for all (without pubchem ±3GB)"),
                                                            textOutput("curr_db_dir"),
                                                            hr(),
                                                            shinyDirButton("get_work_dir", "Choose a working directory" ,
                                                                           title = "Browse",
                                                                           buttonType = "default", class = NULL),
                                                            helpText("Your results will be stored here for later access."),
                                                            textOutput("exp_dir")
                                                   ),
                                                   tabPanel("Graphics", icon=icon("paint-brush"),
                                                            h2("Continuous data"),
                                                            selectInput("color_ramp", label = "Color scheme", choices = list("RAINBOW!"="rb",
                                                                                                                             "Yellow - blue"="y2b",
                                                                                                                             "Matlab 1"="ml1",
                                                                                                                             "Matlab 2 "="ml2",
                                                                                                                             "Magenta - Green"="m2g",
                                                                                                                             "Cyan - yellow"="c2y",
                                                                                                                             "Blue - yellow"="b2y",
                                                                                                                             "Green - red"="g2r",
                                                                                                                             "Blue - green"="b2g",
                                                                                                                             "Blue - red"="b2r",
                                                                                                                             "Blue - pink (pastel)"="b2p",
                                                                                                                             "Blue - green - yellow"="bgy",
                                                                                                                             "Green - yellow - white"="gyw",
                                                                                                                             "Red - yellow - white"="ryw")
                                                            ),
                                                            fluidRow(plotlyOutput("ramp_plot",inline = T, width = "400px", height="200px")),
                                                            h2("Discrete data"),
                                                            uiOutput("colourPickers")
                                                   ),
                                                   tabPanel("Statistics", icon=icon("bar-chart-o"), 
                                                            textInput(inputId="ppm", label="PPM", value = '', width = "50px"),
                                                            actionButton("set_ppm", label="Apply"),
                                                            helpText("Parts per million window for matching. Please re-import outlists and excel file after changing this!"),
                                                            textOutput("ppm"))
                                        )
                                        
                                        
                               ))}
))
                  