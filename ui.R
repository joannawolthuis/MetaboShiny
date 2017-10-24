# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).

shinyUI(fluidPage(theme = "metaboshiny.css",
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
                                                                        DT::dataTableOutput('package_tab',width="100%"),style='font-size:80%'),
                                                                        br(),                                  
                                                                        actionButton("install_packages", "Install", icon = icon("wrench")),
                                                                        br(),br(),
                                                                        imageOutput("package_check")
                                        ))
                               )
                  } else{
                    navbarPage(title=h1("MetaboShiny"), inverse = T,
                               id="nav_general",
                               windowTitle = "MetaboShiny",
                               tabPanel("", icon = icon("share-alt"), value="setup",
                                        # --- db check cols ---
                                        fluidRow(column(width=2),column(width=8, align="center",
                                                                        h2("Setup"),
                                                                        br(),
                                                                        imageOutput("cute_package",inline = T),
                                                                        hr(),
                                                                        helpText("Installed packages:"),
                                                                        div(DT::dataTableOutput('package_tab', width="100%"),style='font-size:80%'),
                                                                        br(),                                  
                                                                        #                                        actionButton("install_packages", "Install", icon = icon("wrench")),
                                                                        actionButton("update_packages", "Update", icon = icon("star")),
                                                                        br(),br(),
                                                                        imageOutput("package_check")
                                        ))),
                               # --------------------------------------------------------------------------------------------------------------------------------
                               tabPanel("", icon = icon("database"), value="database",
                                        # -- header row ---
                                        fluidRow(column(12, align="center",
                                                        h3("")
                                        )),
                                        # --- db check cols --- common pollutants found in DIMS.
                                        fluidRow(column(3, align="center",
                                                        h2("UMC Internal"),
                                                        helpText("Internal commonly known metabolites.")
                                        ),column(3, align="center",
                                                 h2("UMC Noise"),
                                                 helpText("Internal common pollutants found in DIMS using local method.")
                                        ),
                                        column(3,  align="center",
                                               h2("HMDB"),
                                               helpText("Metabolites commonly found in human biological samples.")
                                        ),column(3, align="center",
                                                h2("ChEBI"),
                                                helpText("A broad database with known chemicals of biological interest.")
                                        )), br(),
                                        fluidRow(column(3, align="center",
                                                        imageOutput("umc_logo_int",inline = T),
                                                        br(),br()
                                        ),column(3, align="center",
                                                 imageOutput("umc_logo_noise",inline = T),
                                                 br(),br()
                                        ),
                                        column(3,  align="center",
                                               imageOutput("hmdb_logo",inline = T),
                                               br(),br()
                                               ),
                                        column(3, align="center",
                                               imageOutput("chebi_logo",inline = T),
                                               br(),br()
                                        )),
                                        fluidRow(column(3, align="center",
                                                        actionButton("check_internal", "Check", icon = icon("check")),
                                                        actionButton("build_internal", "Build", icon = icon("wrench")),
                                                        br(),br(),
                                                        imageOutput("internal_check",inline = T)
                                        ),column(3, align="center",
                                                 actionButton("check_noise", "Check", icon = icon("check")),
                                                 actionButton("build_noise", "Build", icon = icon("wrench")),
                                                 br(),br(),
                                                 imageOutput("noise_check",inline = T)
                                        ),
                                        column(3,  align="center",
                                               actionButton("check_hmdb", "Check", icon = icon("check")),
                                               actionButton("build_hmdb", "Build", icon = icon("wrench")),
                                               br(),br(),
                                               imageOutput("hmdb_check",inline = T)),
                                        column(3, align="center",
                                               actionButton("check_chebi", "Check", icon = icon("check")),
                                               actionButton("build_chebi", "Build", icon = icon("wrench")),
                                               br(),br(),
                                               imageOutput("chebi_check",inline = T)
                                        )),
                                        br(),br(),
                                        # --- db check cols --- common pollutants found in DIMS.
                                        fluidRow(column(3, align="center",
                                                 h2("WikiPathways"),
                                                 helpText("Open source biological pathway database. Currently only partially available."),
                                                 h3(" - REQUIRES CHEBI - ")
                                        ),
                                        column(3,  align="center",
                                               h2("KEGG"),
                                               helpText("Large pathway database with info on pathways in various organisms, involved enzymes, and connected disease phenotypes.")
                                        ),
                                        column(3,  align="center",
                                               h2("SMPDB"),
                                               helpText("The Small Molecule Pathway Database is a database containing > 700 small molecule pathways found in humans."),
                                               h3(" - REQUIRES HMDB - ")
                                        )), br(),
                                        fluidRow(column(3, align="center",
                                                 imageOutput("wikipath_logo",inline = T),
                                                 br(),br()
                                        ),
                                        column(3,  align="center",
                                               imageOutput("kegg_logo",inline = T),
                                               br(),br()
                                               ),
                                        column(3,  align="center",
                                               imageOutput("smpdb_logo",inline = T),
                                               br(),br()
                                        )),
                                        fluidRow(column(3, align="center",
                                                 actionButton("check_wikipathways", "Check", icon = icon("check")),
                                                 actionButton("build_wikipathways", "Build", icon = icon("wrench")),
                                                 br(),br(),
                                                 imageOutput("wikipathways_check",inline = T)
                                        ),
                                        column(3,  align="center",
                                               actionButton("check_kegg", "Check", icon = icon("check")),
                                               actionButton("build_kegg", "Build", icon = icon("wrench")),
                                               br(),br(),
                                               imageOutput("kegg_check",inline = T)),
                                        column(3,  align="center",
                                               actionButton("check_smpdb", "Check", icon = icon("check")),
                                               actionButton("build_smpdb", "Build", icon = icon("wrench")),
                                               br(),br(),
                                               imageOutput("smpdb_check",inline = T)))
                                       # ----------------------------------------------------
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
                                        sidebarLayout(position="left",
                                                      sidebarPanel = sidebarPanel(align="center",
                                                        selectInput('exp_var', 'Which experimental variable do you want to look at?', choices = c("")),
                                                        actionButton("check_excel", "Get variables", icon=icon("refresh")),
                                                        br(),br(),
                                                        selectInput('exp_type', 'What type of analysis do you want to do?', choices = list("Standard" = "stat.standard",
                                                                                                                                           "Time series - standard" = "time.standard", 
                                                                                                                                           "Time series - custom end" = "time.custom",
                                                                                                                                           "Time series - subtract" = "time.subtract"
                                                        )),
                                                        uiOutput("exp_opt"),
                                                        # ------------------
                                                          tabsetPanel( id = "adductSettings", selected="db",
                                                                       tabPanel(icon("database"), value="db",
                                                                                br(),
                                                                                tags$i("Select database(s)"),
                                                                                br(),
                                                                                fluidRow(
                                                                                  sardine(fadeImageButton("add_internal", img.path = "umcinternal.png")),
                                                                                  sardine(fadeImageButton("add_noise", img.path = "umcnoise.png")),
                                                                                  sardine(fadeImageButton("add_hmdb", img.path = "hmdblogo.png")),
                                                                                  sardine(fadeImageButton("add_chebi", img.path = "chebilogo.png")),
                                                                                  sardine(fadeImageButton("add_smpdb", img.path = "smpdb_logo_adj.png")),
                                                                                  sardine(fadeImageButton("add_wikipathways", img.path = "wikipathways.png")),
                                                                                  sardine(fadeImageButton("add_kegg", img.path = "kegglogo.gif", value = T))
                                                                                )),
                                                                       tabPanel(icon("id-card-o"), value = "identifier",
                                                                                br(),
                                                                                tags$i("Select identifier"),
                                                                                radioButtons(inputId = "group_by", label = NULL, choices = 
                                                                                               list("Pathway ID" = "pathway",
                                                                                                    "Database ID" = "identifier", 
                                                                                                    "Compound name" = "compoundname",
                                                                                                    "Molecular formula" = "baseformula",
                                                                                                    "Mass/charge" = "mz"), 
                                                                                             selected = "baseformula",
                                                                                             width="100%")                             ), 
                                                                       tabPanel(icon("plus-square"), value="adducts",
                                                                                br(),
                                                                                tags$i("Select adduct(s)"),
                                                                                fluidRow(column(width=6, div(style="font-size:120%",icon("search-plus"))), 
                                                                                         column(width=6,div(style="font-size:120%",icon("search-minus")))
                                                                                         ),
                                                                                fluidRow(column(width=6,div(DT::dataTableOutput('pos_add_tab',width="100%"),style='font-size:70%')),
                                                                                         column(width=6,div(DT::dataTableOutput('neg_add_tab',width="100%"),style='font-size:70%'))
                                                                                        )
                                                                       )
                                                        ),
                                                        # ------------------
                                                        actionButton("create_csv", "Create CSV", icon=icon("file-text-o"))
                                                        # ,imageOutput("csv_icon",inline = T),
                                                        # br(),br(),
                                                        # fileInput("pat_csv", "Import CSV",
                                                        #           multiple = F,
                                                        #           accept = c(".csv")),
                                                        # actionButton("import_csv", "Import", icon = icon("hand-peace-o")),
                                                        # imageOutput("csv_upload_check",inline = T)
                                                      ),
                                                      mainPanel = mainPanel(align="center",
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
                                                                                                                                        "By reference feature"="ProbNorm",
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
                                                                                  bsCollapse(id = "dbSelect", open = "Settings",
                                                                                             bsCollapsePanel(h3("Databases"), "",style = "info",
                                                                                               fluidRow(#imageOutput("find_mol_icon",inline = T),
                                                                                                 sardine(fadeImageButton("search_internal", img.path = "umcinternal.png")),
                                                                                                 sardine(fadeImageButton("search_noise", img.path = "umcnoise.png")),
                                                                                                 sardine(fadeImageButton("search_hmdb", img.path = "hmdblogo.png")),
                                                                                                 sardine(fadeImageButton("search_chebi", img.path = "chebilogo.png")),
                                                                                                 sardine(fadeImageButton("search_smpdb", img.path = "smpdb_logo_adj.png")),
                                                                                                 sardine(fadeImageButton("search_wikipathways", img.path = "wikipathways.png")),
                                                                                                 sardine(fadeImageButton("search_kegg", img.path = "kegglogo.gif")),
                                                                                                 #sardine(fadeImageButton("search_pubchem", img.path = "pubchemlogo.png")),
                                                                                                 br(),
                                                                                                 sardine(switchButton(inputId = "autosearch",
                                                                                                                      label = "Autosearch", 
                                                                                                                      value = FALSE, col = "BW", type = "OO"))
                                                                                               ),
                                                                                               h4("Current compound:"),
                                                                                               verbatimTextOutput("curr_cpd") 
                                                                                             )
                                                                                             ),
                                                        hr(),
                                                        fluidRow(navbarPage("Search", id="tab_iden",
                                                                            tabPanel("Current", icon=icon("sort-numeric-asc"),
                                                                                     actionButton("search_cpd", "Find hits", icon=icon("search")),
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
                                                                                     actionButton("revsearch_cpd", "Find hits", icon=icon("search")),
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
                                                            fluidRow(plotlyOutput("ramp_plot",inline = T, width="100%")),
                                                            h2("Discrete data"),
                                                            uiOutput("colourPickers")
                                                   ),
                                                   tabPanel("Statistics", icon=icon("bar-chart-o"), 
                                                            textInput(inputId="ppm", label="PPM", value = '', width = "50px"),
                                                            actionButton("set_ppm", label="Apply"),
                                                            helpText("Parts per million window for matching. Please re-import outlists and excel file after changing this!"),
                                                            textOutput("ppm")),
                                                   tabPanel("Adducts", icon=icon("plus-square"),
                                                            h3("Current adduct table:"),
                                                            rHandsontableOutput("adduct_tab", width=800, height=600),
                                                            shinySaveButton("save_adducts", 
                                                                            "Save changed table", 
                                                                            "Save file as ...", 
                                                                            filetype=list(RData="RData", csv="csv")
                                                                            ),
                                                            hr(),
                                                            fileInput("add_tab", "Import adduct table",
                                                                      multiple = F,
                                                                      accept = c(".RData", ".csv")),
                                                            sardine(actionButton("import_adducts", "Import", icon = icon("hand-peace-o"))),
                                                            sardine(imageOutput("adduct_upload_check",inline = T))
                                                   )
                                        )
                               ), div(class="spinnylocation1", 
                                      div(class="plus", img(class="imagetop", src="heart-research.png", width="120px", height="120px")),
                                      div(class="minus", img(class="imagebottom", src="heart-research.png", width="120px", height="120px"))
                                      ),
                               div(class="line")
                    )
                    }
                  )
        )
                  