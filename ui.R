#TODO: zoom font size

shinyUI(fluidPage(theme = "metaboshiny.css",tags$head(tags$script(src = "sparkle.js"), 
                                                      tags$style(type="text/css", bar.css)
),
shinyalert::useShinyalert(),
navbarPage(inverse=TRUE,title=div(h1("MetaboShiny"), tags$head(tags$style(type="text/css", font.css)),
                                  class="sparkley"),
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
                                                    div(DT::dataTableOutput('package_tab', 
                                                                            width="100%"),
                                                        style='font-size:80%'),
                                                    br(),                                  
                                                    actionButton("update_packages", "Install missing packages", icon = icon("heart")),
                                                    br(),br(),
                                                    imageOutput("package_check")
                    ))),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("", icon = icon("database"), value="database",
                    # -- header row ---
                    fluidRow(column(12, align="center",
                                    h3("")
                    )), # FIRST ROW
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
                    fluidRow(column(3, align="center", # SECOND ROW
                                    h2("WikiPathways"),
                                    helpText("Open source biological pathway database. Currently only partially available. Requires CHEBI to be built.")
                    ),
                    column(3,  align="center",
                           h2("KEGG"),
                           helpText("Large pathway database with info on pathways in various organisms, involved enzymes, and connected disease phenotypes.")
                    ),
                    column(3,  align="center",
                           h2("SMPDB"),
                           helpText("Small molecule pathway database. Compounds overlap with HMDB.")
                    ),column(3,  align="center",
                             h2("MetaCyc"),
                             helpText("Large pathway database with over 10 000 available compounds. Spans several organisms!"),
                             helpText(a("Download spreadsheet to /backend/db/metacyc_source.", href = 'https://metacyc.org/group?id=biocyc17-31223-3729417004'))
                    )
                    ), br(),
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
                    )
                    ,column(3,  align="center",
                            imageOutput("metacyc_logo",inline = T),
                            br(),br()
                    )
                    ),
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
                           imageOutput("kegg_check",inline = T))
                    ,
                    column(3,  align="center",
                           actionButton("check_smpdb", "Check", icon = icon("check")),
                           actionButton("build_smpdb", "Build", icon = icon("wrench")),
                           br(),br(),
                           imageOutput("smpdb_check",inline = T))
                    ,
                    column(3,  align="center",
                           actionButton("check_metacyc", "Check", icon = icon("check")),
                           actionButton("build_metacyc", "Build", icon = icon("wrench")),
                           br(),br(),
                           imageOutput("metacyc_check",inline = T))
                    ),br(),# THIRD ROW
                    fluidRow(column(3, align="center",
                                    h2("DIMEdb"),
                                    helpText("A direct infusion database of biologically relevant metabolite structures and annotations.")
                    ),column(3, align="center",
                             h2("Wikidata"),
                             helpText("Central storage for the data of its Wikimedia sister projects including Wikipedia, Wikivoyage, Wikisource, and others.")
                    ),column(3, align="center",
                             h2("VMH"),
                             helpText("Virtual Metabolic Human (VMH) hosts ReconMap, an extensive network of human metabolism, and bacterial metabolites.")
                    )
                    
                    ),
                    fluidRow(column(3, align="center",
                                    imageOutput("dimedb_logo",inline = T),
                                    br(),br()
                    ),column(3, align="center",
                             imageOutput("wikidata_logo",inline = T),
                             br(),br()
                    ),column(3, align="center",
                             imageOutput("vmh_logo",inline = T),
                             br(),br()
                    )),
                    fluidRow(column(3, align="center",
                                    actionButton("check_dimedb", "Check", icon = icon("check")),
                                    actionButton("build_dimedb", "Build", icon = icon("wrench")),
                                    br(),br(),
                                    imageOutput("dimedb_check",inline = T)
                    ),column(3, align="center",
                             actionButton("check_wikidata", "Check", icon = icon("check")),
                             actionButton("build_wikidata", "Build", icon = icon("wrench")),
                             br(),br(),
                             imageOutput("wikidata_check",inline = T)
                    ),
                    column(3, align="center",
                           actionButton("check_vmh", "Check", icon = icon("check")),
                           actionButton("build_vmh", "Build", icon = icon("wrench")),
                           br(),br(),
                           imageOutput("vmh_check",inline = T)
                    ))
                    # ----------------------------------------------------
           ),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("", icon = icon("upload"), value="upload", ## this guy gives error???
                    fluidRow(column(9, align="center", 
                                    h2("Create project"))),
                    hr(),
                    tabsetPanel(id="new_proj",selected = "From DB",
                               tabPanel(id="db", title="From DB",
                                        br(),br(),
                                        fluidRow(
                                          column(3, align="center",
                                                        imageOutput("db_icon", inline = T),
                                                        br(),br()
                                                ),
                                                column(3, align="center",
                                                       br(),br(),
                                                        imageOutput("excel_icon",inline = T)
                                                        )
                                                 ), 
                                        fluidRow(
                                          column(3, align="center",
                                                 shinyFilesButton('database', 'DATABASE', 'Please select an database file', FALSE)
                                          ),
                                          column(3, align="center",
                                                 shinyFilesButton('excel', 'METADATA', 'Please select an excel file', FALSE)
                                          )
                                        )
                                        ),
                               tabPanel(id="csv", title="From CSV",
                                        br(),br(),
                                        fluidRow(column(3,  align="center",
                                                        imageOutput("pos_icon",inline = T),
                                                        br(),br(),
                                                        shinyFilesButton('outlist_pos', 'POSITIVE PEAKS', 'Please select a csv file', FALSE)
                                        ),
                                        column(3,  align="center",
                                               imageOutput("excel_icon_2",inline = T),
                                               br(),br(),
                                               shinyFilesButton('excel', 'METADATA', 'Please select an excel file', FALSE)
                                        ),
                                        column(3,  align="center",
                                               imageOutput("neg_icon",inline = T),
                                               br(),br(),
                                               shinyFilesButton('outlist_neg', 'NEGATIVE PEAKS', 'Please select a csv file', FALSE)
                                        )
                                        )
                                        )
                    ),
                    hr(),
                    fluidRow(column(9, align="center",
                                    actionButton("create_db", "Go", icon = icon("magic"))
                    ))
           ),
           tabPanel("", icon=icon("file-text-o"), value="document",
                    sidebarLayout(position="left",
                                  sidebarPanel = sidebarPanel(align="center",
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
                                                                                      sardine(fadeImageButton("add_metacyc", img.path = "metacyc.png")),
                                                                                      sardine(fadeImageButton("add_wikipathways", img.path = "wikipathways.png")),
                                                                                      sardine(fadeImageButton("add_kegg", img.path = "kegglogo.gif", value = T)),
                                                                                      sardine(fadeImageButton("add_dimedb", img.path = "dimedb.png")),
                                                                                      sardine(fadeImageButton("add_wikidata", img.path = "wikidata.png")),
                                                                                      sardine(fadeImageButton("add_vmh", img.path = "vmh.png"))
                                                                                    )),
                                                                           tabPanel(icon("id-card-o"), value = "identifier",
                                                                                    br(),
                                                                                    tags$i("Select identifier"),
                                                                                    radioButtons(inputId = "group_by", label = NULL, choices = 
                                                                                                   list("Molecular formula" = "baseformula",
                                                                                                        "Mass/charge" = "mz"), 
                                                                                                 selected = "mz",
                                                                                                 width="100%")                             ), 
                                                                           tabPanel(icon("plus-square"), value="adducts",
                                                                                    br(),
                                                                                    tags$i("Select adduct(s)"),
                                                                                    fluidRow(column(width=6, div(style="font-size:120%",icon("search-plus"))), 
                                                                                             column(width=6,div(style="font-size:120%",icon("search-minus")))
                                                                                    ),
                                                                                    fluidRow(column(width=6,div(DT::dataTableOutput('pos_add_tab',width="100%"),style='font-size:70%')),
                                                                                             column(width=6,div(DT::dataTableOutput('neg_add_tab',width="100%"),style='font-size:70%'))
                                                                                    ),
                                                                                    fluidRow(sardine(div(actionButton(inputId = "sel_all_adducts",
                                                                                                                      label = "", icon=icon("circle")),style='font-size:70%')),
                                                                                             sardine(div(actionButton(inputId = "sel_comm_adducts",
                                                                                                                      label = "", icon=icon("check-circle-o")),style='font-size:70%')),
                                                                                             sardine(div(actionButton(inputId = "sel_no_adducts",
                                                                                                                      label = "", icon=icon("circle-o")),style='font-size:70%'))
                                                                                    ),br()
                                                                                    
                                                                           )
                                                              ),
                                                              # ------------------
                                                              actionButton("create_csv", "Create CSV", icon=icon("file-text-o")),
                                                              hr(),
                                                              imageOutput("csv_icon",inline = T),
                                                              br(),br(),
                                                              fileInput("pat_csv", "Import CSV",
                                                                        multiple = F,
                                                                        accept = c(".csv")),
                                                              actionButton("import_csv", "Import", icon = icon("hand-peace-o")),
                                                              imageOutput("csv_upload_check",inline = T)
                                  ),
                                  mainPanel = mainPanel(align="center",
                                                        fluidRow(div(DT::dataTableOutput('csv_tab'),style='font-size:80%'))
                                  )
                    )
           ),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("",  icon = icon("shower"), value="filter",
                    fluidRow(column(3, aligh="center",
                                    selectInput('samp_var', 'Which variable represents sample amount/concentration?', choices = c("")), #TODO: only show this when normalize by sample specific factor (specnorm) is selected
                                    selectizeInput('batch_var', 'What are your batch variables?', choices = c("batch"), multiple=TRUE, options = list(maxItems = 2L)),
                                    actionButton("check_csv", "Get options", icon=icon("refresh")),
                                    hr(),
                                    #sliderInput("perc_limit", label = "Max. missing feature percent", min = 0, 
                                    #            max = 100, value = 20),
                                    shinyWidgets::sliderTextInput("perc_limit","Max. missing feature percent:",
                                                                  choices=c(0, 0.0001, 0.001, 0.01, 0.1, seq(1, 100, 1)),
                                                                  selected=0, grid = T),
                                    selectInput('filt_type', 'How will you filter your m/z values?', choices = list("Interquantile range"="iqr",
                                                                                                                    "Relative stdev"="rsd",
                                                                                                                    "Non-parametric relative stdev"="nrsd",
                                                                                                                    "Mean"="mean",
                                                                                                                    "Standard deviation"="sd",
                                                                                                                    "Median absolute deviation"="mad",
                                                                                                                    "Median"="median",
                                                                                                                    "None"="none"),
                                                selected = "none"),
                                    selectInput('norm_type', 'What type of normalization do you want to do?', choices = list("Quantile normalization"="QuantileNorm",
                                                                                                                             "By reference feature"="ProbNorm",
                                                                                                                             "By reference compound"="CompNorm",
                                                                                                                             "By sample specific factor"="SpecNorm",
                                                                                                                             "Sum"="SumNorm",
                                                                                                                             "Median"="MedianNorm",
                                                                                                                             "None"="NULL")),
                                    uiOutput("ref_select"),
                                    selectInput('trans_type', 'How will you transform your data?', choices = list("Log transform"="LogNorm",
                                                                                                                  "Cubic root transform"="CrNorm",
                                                                                                                  "None"="NULL")),
                                    selectInput('scale_type', 'How will you scale your data?', choices = list("Autoscale"="AutoNorm",
                                                                                                              "Mean-center"="MeanCenter",
                                                                                                              "Pareto Scaling"="ParetoNorm",
                                                                                                              "Range scaling"="RangeNorm",
                                                                                                              "None"="NULL")),
                                    selectInput('miss_type', 'How to deal with missing values?', choices = list("Feature minimum"="colmin",
                                                                                                                "Total minimum"="min",
                                                                                                                "Random forest"="rf",
                                                                                                                #"Impute w/ regression"="regr",
                                                                                                                "KNN imputation"="knn",
                                                                                                                "SVD imputation"="svdImpute",
                                                                                                                "BPCA imputation"="bpca",
                                                                                                                "PPCA imputation"="ppca",
                                                                                                                "Median"="median",
                                                                                                                "Mean"="mean",
                                                                                                                "Leave them out"="exclude",
                                                                                                                "Leave them alone"="none"), 
                                                selected = "knn"),
                                    switchButton(inputId = "remove_outliers",
                                                                          label = "Exclude outliers?", 
                                                                          value = TRUE, col = "BW", type = "YN"),
                                    actionButton("initialize", "Go", icon=icon("hand-o-right")),
                                    hr(),
                                    imageOutput("dataset_icon",inline = T),
                                    fileInput("pat_dataset", "Import dataset",
                                              multiple = F,
                                              accept = c(".RData")),
                                    actionButton("import_dataset", "Import", icon = icon("hand-peace-o")),
                                    imageOutput("dataset_upload_check",inline = T)
                    ), column(9,
                              navbarPage(inverse=TRUE,"Explore",
                                         tabPanel("Variables", icon=icon("braille"),
                                                  fluidRow(column(6,plotOutput("var1",height='300px')),
                                                           column(6,plotOutput("var3", height='300px'))
                                                  ),
                                                  fluidRow(column(6,plotOutput("var2", height='500px')),
                                                           column(6,plotOutput("var4", height='500px')))
                                         ),
                                         tabPanel("Samples", icon=icon("tint"),
                                                  fluidRow(column(6,plotOutput("samp1",height='300px')),
                                                           column(6,plotOutput("samp3", height='300px'))
                                                  ),
                                                  fluidRow(column(6,plotOutput("samp2", height='500px')),
                                                           column(6,plotOutput("samp4", height='500px')))
                                         )
                              ))
                    )),
           # --------------------------------------------------------------------------------------------------------------------------------
           tabPanel("",  icon = icon("bar-chart"), value = "analysis",# multiple options here :-)
                    sidebarLayout(position="right",
                                  mainPanel = mainPanel(width = 8,
                                                        navbarPage(inverse=F,h2("Statistics"), id="statistics",
                                                                   # TODO: T-SNE
                                                                   tabPanel(h3("Info"), value = "inf",
                                                                            fluidRow(column(width=12, align="center",
                                                                                            br(),br(),br(),
                                                                                            icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg"),
                                                                                            br(),br(),
                                                                                            h2("Please select a variable of interest in the sidebar!"),
                                                                                            br(),br(),
                                                                                            icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg")
                                                                            ))),
                                                                   tabPanel(h3("PCA"), value = "pca", #icon=icon("cube"),
                                                                            fluidRow(column(10,align="center",plotly::plotlyOutput("plot_pca",height = "600px", width="600px"))
                                                                            ),
                                                                            fluidRow(column(12,align="center",
                                                                                            switchButton("pca_2d3d", label = "", col = "BW", type = "2d3d"))),
                                                                            hr(),
                                                                            fluidRow(column(3,
                                                                                            selectInput("pca_x", label = "X axis:", choices = paste0("PC",1:20),selected = "PC1",width="100%"),
                                                                                            selectInput("pca_y", label = "Y axis:", choices = paste0("PC",1:20),selected = "PC2",width="100%"),
                                                                                            selectInput("pca_z", label = "Z axis:", choices = paste0("PC",1:20),selected = "PC3",width="100%")),
                                                                                     column(9,
                                                                                            tabsetPanel(id="pca_2", 
                                                                                                        tabPanel(title="Table", 
                                                                                                                 div(DT::dataTableOutput('pca_tab',width="100%"),style='font-size:80%')),
                                                                                                        tabPanel(title="Scree",
                                                                                                                 plotOutput("pca_scree")
                                                                                                        ),
                                                                                                        tabPanel(title="Loadings", 
                                                                                                                 div(DT::dataTableOutput('pca_load_tab',width="100%"),style='font-size:80%'))
                                                                                            ))
                                                                            )
                                                                   ),
                                                                   tabPanel(h3("PLSDA"), value = "plsda", 
                                                                            fluidRow(column(12,align="center",plotly::plotlyOutput("plot_plsda",height = "500px", width="500px"))),
                                                                            fluidRow(column(12,align="center",switchButton("plsda_2d3d", label = "", col = "BW", type = "2d3d"))),
                                                                            hr(),
                                                                            fluidRow(column(3,
                                                                                            div(style="display:inline-block",
                                                                                                selectInput("plsda_type", 
                                                                                                            label="Type:", 
                                                                                                            choices=list("Normal"="normal",
                                                                                                                         "Orthogonal"="ortho",
                                                                                                                         "Sparse"="sparse"), 
                                                                                                            width = '100px',
                                                                                                            selected=1)),
                                                                                            div(style="display:inline-block",
                                                                                                shinyWidgets::circleButton("do_plsda", icon = icon("hand-pointer-o"), size = "sm")
                                                                                            ),
                                                                                            selectInput("plsda_x", label = "X axis:", choices = paste0("PC",1:8),selected = "PC1",width="100%"),
                                                                                            selectInput("plsda_y", label = "Y axis:", choices = paste0("PC",1:8),selected = "PC2",width="100%"),
                                                                                            selectInput("plsda_z", label = "Z axis:", choices = paste0("PC",1:8),selected = "PC3",width="100%")),
                                                                                     column(9,
                                                                                            tabsetPanel(id="plsda_2", 
                                                                                                        tabPanel(title="Cross-validation", 
                                                                                                                 plotOutput("plsda_cv_plot")),
                                                                                                        tabPanel(title="Permutation", 
                                                                                                                 plotOutput("plsda_perm_plot")),
                                                                                                        tabPanel(title="Table", 
                                                                                                                 div(DT::dataTableOutput('plsda_tab',width="100%"),style='font-size:80%')),
                                                                                                        tabPanel(title="Loadings", 
                                                                                                                 div(DT::dataTableOutput('plsda_load_tab',width="100%"),style='font-size:80%'))
                                                                                            ))
                                                                            )
                                                                   ),
                                                                   tabPanel(h3("T-test"), value="tt", 
                                                                            fluidRow(plotly::plotlyOutput('tt_specific_plot',width="100%")),
                                                                            navbarPage(inverse=F,"",
                                                                                       tabPanel("", icon=icon("table"),
                                                                                                div(DT::dataTableOutput('tt_tab',width="100%"),style='font-size:80%'))
                                                                                       ,tabPanel("", icon=icon("area-chart"),
                                                                                                 plotly::plotlyOutput('tt_overview_plot',height="300px")
                                                                                       )
                                                                            )),
                                                                   tabPanel(h3("ANOVA"), value="aov",
                                                                            fluidRow(plotly::plotlyOutput('aov_specific_plot',width="100%")),
                                                                            navbarPage(inverse=F,"",
                                                                                       tabPanel("", icon=icon("table"),
                                                                                                div(DT::dataTableOutput('aov_tab',width="100%"),style='font-size:80%'))
                                                                                       ,tabPanel("", icon=icon("area-chart"),
                                                                                                 plotly::plotlyOutput('aov_overview_plot',height="300px")
                                                                                       )
                                                                            )),
                                                                   tabPanel(h3("Fold-change"), value="fc",
                                                                            fluidRow(plotly::plotlyOutput('fc_specific_plot',width="100%")),
                                                                            navbarPage(inverse=F,"",
                                                                                       tabPanel("", icon=icon("table"),
                                                                                                div(DT::dataTableOutput('fc_tab',width="100%"),style='font-size:80%'))
                                                                                       ,tabPanel("", icon=icon("area-chart"),
                                                                                                 plotly::plotlyOutput('fc_overview_plot',height="300px")
                                                                                       ))
                                                                   ),
                                                                   tabPanel(h3("Volcano"), value="volc",
                                                                            fluidRow(plotly::plotlyOutput('volc_plot',width="100%",height="600px"))
                                                                   ),
                                                                   tabPanel(h3("MEBA"), value="meba", 
                                                                            fluidRow(plotly::plotlyOutput('meba_specific_plot'),height="600px"),
                                                                            fluidRow(div(DT::dataTableOutput('meba_tab', width="100%"),style='font-size:80%'))
                                                                   ),
                                                                   # =================================================================================
                                                                   tabPanel(h3("ASCA"), value="asca",
                                                                            fluidRow(plotly::plotlyOutput('asca_specific_plot', height="600px")),
                                                                            fluidRow(div(DT::dataTableOutput('asca_tab',width="100%"),style='font-size:80%'))
                                                                   ),
                                                                   tabPanel(h3("Heatmap"), value="heatmap",
                                                                            plotly::plotlyOutput("heatmap",width="100%",height="700px"),
                                                                            br(),
                                                                            fluidRow(column(align="center",
                                                                                            width=12,
                                                                                            sliderInput("heatmap_topn", "Use top ... from table:", value=100, min = 10, max = 200))
                                                                            ),
                                                                            fluidRow(column(align="center",
                                                                                            width=12,
                                                                                            uiOutput("heatbutton"))
                                                                            )
                                                                   )
                                                                   # ,tabPanel(h3("Enrichment"), value="enrich_biv",
                                                                   #          sidebarLayout(position="left",
                                                                   #                        sidebarPanel = sidebarPanel(align="center",
                                                                   #                                                    fluidRow(
                                                                   #                                                      tags$b("Pathway DBs"),br(),
                                                                   #                                                      sardine(fadeImageButton("enrich_metacyc", 
                                                                   #                                                                              img.path = "metacyc.png", 
                                                                   #                                                                              value = T)),
                                                                   #                                                      sardine(fadeImageButton("enrich_wikipathways", 
                                                                   #                                                                              img.path = "wikipathways.png", 
                                                                   #                                                                              value = T)),
                                                                   #                                                      sardine(fadeImageButton("enrich_kegg", 
                                                                   #                                                                              img.path = "kegglogo.gif", 
                                                                   #                                                                              value = T))
                                                                   #                                                    ),
                                                                   #                                                    selectInput('enrich_stats', 
                                                                   #                                                                'Score source',
                                                                   #                                                                choices = c("T-test"="tt", 
                                                                   #                                                                            "Fold-change"="fc",
                                                                   #                                                                            "PLS-DA"="plsda",
                                                                   #                                                                            "RandomForest"="rf")
                                                                   #                                                    ),
                                                                   #                                                    selectInput('enrich_vals', 
                                                                   #                                                                'Score threshold',
                                                                   #                                                                choices = c("Significant"="sig", 
                                                                   #                                                                            "Top 50"="t50", 
                                                                   #                                                                            "Top 100"="t100", 
                                                                   #                                                                            "Top 200"="t200", 
                                                                   #                                                                            "Top 500"="t500")),
                                                                   #                                                    hr(),
                                                                   #                                                    actionButton("go_enrich", "Analyse", icon=icon("binoculars"))
                                                                   #                        ),
                                                                   #                        # ------------------
                                                                   #                        mainPanel = mainPanel(align="center",
                                                                   #                                              fluidRow(div(DT::dataTableOutput('enriched'),style='font-size:80%')),
                                                                   #                                              fluidRow(div(DT::dataTableOutput('enrich_pw_tab'),style='font-size:80%')),
                                                                   #                                              fluidRow(plotly::plotlyOutput('enrich_pw_specific_plot'),height="200px")
                                                                   #                                              
                                                                   #                        ))
                                                                   # )
                                                                   ,tabPanel(h3("ML"), value = "ml",
                                                                             fluidRow(
                                                                               column(width=3,align="center",
                                                                                      selectInput("ml_method", 
                                                                                                  label = "Type:", 
                                                                                                  selected = "rf", 
                                                                                                  choices = list("Random Forest" = "rf",
                                                                                                                 "Lasso" = "ls",
                                                                                                                 "Group lasso" = "gls"),
                                                                                                  multiple = F),
                                                                                      sliderInput("ml_train_perc", 
                                                                                                  label = "Percentage in training", 
                                                                                                  min = 1,
                                                                                                  max = 100,
                                                                                                  step = 1,
                                                                                                  value = 60, 
                                                                                                  post = "%"),
                                                                                      selectInput("ml_folds", label="Fold CV",choices = c("5", 
                                                                                                                                          "10", 
                                                                                                                                          "20", 
                                                                                                                                          "50", 
                                                                                                                                          "LOOCV"),
                                                                                                  multiple = F),
                                                                                      # - - - - - - - - - -
                                                                                      style="z-index:1002;"
                                                                               ),
                                                                               column(width=6,align="center",
                                                                                      sliderInput("ml_attempts", 
                                                                                                  label = "Attempts", 
                                                                                                  min = 1,
                                                                                                  max = 100,
                                                                                                  step = 1,
                                                                                                  value = 20, 
                                                                                                  post = "x"),
                                                                                      br(),
                                                                                      br(),
                                                                                      shinyWidgets::circleButton("do_ml",icon = h3(paste("Go"), icon("hand-pointer-o", "fa-lg")), status = "default", size = "lg")
                                                                                      #actionButton("do_ml",label=h2("Go"),width = "150px", height="150px")
                                                                               ),
                                                                               column(width=3,align="center",
                                                                                      textInput("ml_train_regex", label = "Regex for train:"),
                                                                                      textInput("ml_test_regex", label = "Regex for test:"),
                                                                                      tags$a(img(src="help.png"), href="https://regex101.com"),
                                                                                      textInput("ml_name", label="Name:"))
                                                                             ),
                                                                             hr(),
                                                                             navbarPage(title="Results",id="ml_results",inverse=F,
                                                                                        tabPanel(title = "ROC",value = "roc",icon=icon("area-chart"),
                                                                                                 plotlyOutput("ml_roc",height = "600px"),
                                                                                                 div(DT::dataTableOutput("ml_tab",width="100%"),style='font-size:80%')),
                                                                                        tabPanel("Model",value= "bar",icon=icon("table"),
                                                                                                 fluidRow(plotlyOutput("ml_bar", width = "100%", height="300px")),
                                                                                                 fluidRow(
                                                                                                   column(5, sliderInput("ml_top_x",
                                                                                                                         label = "Show top:",
                                                                                                                         min = 10,
                                                                                                                         max = 100,
                                                                                                                         step=10,
                                                                                                                         value=20)),
                                                                                                   column(7, plotlyOutput("ml_specific_plot", 
                                                                                                                          height="300px"))
                                                                                                 )
                                                                                        )
                                                                             )
                                                                   ))
                                  ),
                                  sidebarPanel = 
                                    sidebarPanel(align="center",width = 4,
                                                 tabsetPanel(id = "search", #type = "pills",
                                                             tabPanel(NULL, icon=icon("exchange"),
                                                                      br()
                                                                      ,h2("Change variable of interest")
                                                                      ,selectInput("first_var", label="Do statistics on:", choices = c("label"))
                                                                      ,shinyWidgets::circleButton("change_cls", icon = icon("hand-pointer-o"), size = "sm")
                                                                      ,fluidRow(column(12, align="center", uiOutput("timebutton")))
                                                                      ,h2("Shape")
                                                                      ,selectInput("second_var", label="Marker shape based on:", choices = c("label"))
                                                                      ,hr()
                                                                      ,h2("Change sample subset"),
                                                                      textInput("subset_regex", "Subset data based on regex:", value = "...")
                                                                      ,selectInput("subset_var", label="Subset data based on:", choices = c("label"))
                                                                      
                                                             ),
                                                             tabPanel(title=NULL, icon=icon("search"),
                                                                      br(),
                                                                      bsCollapse(id = "dbSelect", open = "Settings",
                                                                                 bsCollapsePanel(h2("Databases"), "",#style = "info",
                                                                                                 fluidRow(#imageOutput("find_mol_icon",inline = T),
                                                                                                   sardine(fadeImageButton("search_internal", img.path = "umcinternal.png")),
                                                                                                   sardine(fadeImageButton("search_noise", img.path = "umcnoise.png")),
                                                                                                   sardine(fadeImageButton("search_hmdb", img.path = "hmdblogo.png")),
                                                                                                   sardine(fadeImageButton("search_chebi", img.path = "chebilogo.png")),
                                                                                                   sardine(fadeImageButton("search_smpdb", img.path = "smpdb_logo_adj.png")),
                                                                                                   sardine(fadeImageButton("search_metacyc", img.path = "metacyc.png")),
                                                                                                   sardine(fadeImageButton("search_wikipathways", img.path = "wikipathways.png")),
                                                                                                   sardine(fadeImageButton("search_kegg", img.path = "kegglogo.gif")),
                                                                                                   sardine(fadeImageButton("search_dimedb", img.path = "dimedb.png")),
                                                                                                   sardine(fadeImageButton("search_wikidata", img.path = "wikidata.png")),
                                                                                                   sardine(fadeImageButton("search_vmh", img.path = "vmh.png"))
                                                                                                 )),
                                                                                 bsCollapsePanel(h2("Miniplot"), "",#style = "info",
                                                                                                 plotly::plotlyOutput("curr_plot", height="300px", width="100%"))
                                                                                 ),
                                                                      tabsetPanel(id="tab_iden",                                                             
                                                                                  tabPanel(title="mz > molecule",
                                                                                           br(),
                                                                                           bsCollapse(id="isoSelect",
                                                                                                      bsCollapsePanel(h2("Isotope scoring"), "",#style = "info",
                                                                                                                      selectInput("iso_score_method", "Which method used to score compounds of same weight?", selected="mscore", 
                                                                                                                                  choices=list("M-score"="mscore",
                                                                                                                                               "Chi-square"="chisq",
                                                                                                                                               "Mean absolute percentage error"="mape",
                                                                                                                                               "SIRIUS"="sirius",
                                                                                                                                               "Network-based"="network"))
                                                                                                      )
                                                                                           ),
                                                                                           hr(),
                                                                                           div(h2(textOutput("curr_cpd"),style="padding:10px;"),
                                                                                               style="background-color:white;
                                                                                               height:55px;
                                                                                               width:115%;
                                                                                               position:relative;
                                                                                               right:30px;
                                                                                               border-top: 1px solid #DFDCDC;
                                                                                               border-bottom: 1px solid #DFDCDC;
                                                                                               "),
                                                                                           fluidRow(
                                                                                             actionButton("search_cpd", "S E A R C H", icon=icon("search"),width = "40%"),
                                                                                             actionButton("score_iso", "S C O R E", icon=icon("balance-scale"), width="40%"),
                                                                                             #,actionButton("score_net", "", icon=icon("arrows-alt")),
                                                                                             align="center"),br(),
                                                                                           fluidRow(
                                                                                             hr(),
                                                                                             bsCollapsePanel(h2("Structure viewer"),"",
                                                                                                             textOutput("curr_formula"),
                                                                                                             plotOutput("curr_struct", height="310px")
                                                                                             ),
                                                                                             # scoring buttons here..
                                                                                             div(DT::dataTableOutput('match_tab', width="100%"),style='font-size:80%'),
                                                                                             hr(),
                                                                                             div(textOutput("curr_definition"))
                                                                                           )),
                                                                                  tabPanel(title="molecule > mz",
                                                                                           br(),
                                                                                           actionButton("browse_db", "Browse compounds", icon=icon("eye")),
                                                                                           hr(),
                                                                                           div(DT::dataTableOutput('browse_tab'),style='font-size:80%'),
                                                                                           hr(),
                                                                                           div(textOutput("browse_definition"),style='font-size:80%'),
                                                                                           hr(),
                                                                                           actionButton("revsearch_cpd", "Find hits", icon=icon("search")),
                                                                                           hr(),
                                                                                           div(DT::dataTableOutput('hits_tab'),style='font-size:80%')
                                                                                  ))
                                                                      ),
                                                             tabPanel(title=NULL, icon=icon("link"),
                                                                      br(),
                                                                      helpText("You can use this tab to find overlapping hits between analyses."),
                                                                      # - - - select candidates - - -
                                                                      
                                                                      checkboxGroupInput("venn_members", label = "Please choose which to combine", choices = list(
                                                                        "T-test" = "tt",
                                                                        "Fold-change" = "fc",
                                                                        "Lasso" = "ls",
                                                                        "Random Forest" = "rf",
                                                                        "Volcano" = "volc",
                                                                        "PLS-DA" = "plsda"
                                                                      ), selected = c("tt", "fc"), inline = T),
                                                                      
                                                                      uiOutput("venn_ml_ui"),
                                                                      
                                                                      sliderInput("venn_tophits", label = "Only include top", min = 1, max = 1000, post = " hits", value=20),
                                                                      
                                                                      # - - - - - - - - - - - - - - -
                                                                      actionButton("build_venn", label = "Find overlap"),
                                                                      plotly::plotlyOutput("venn_plot",inline = F),
                                                                      selectInput("intersect_venn", label = "Show overlap between:", selected = 1,choices = "",multiple = T),
                                                                      div(DT::dataTableOutput('venn_tab'),style='font-size:80%')
                                                                      ),
                                                             tabPanel(NULL, icon=icon("paint-brush"),
                                                                      h2("Plot theme"),
                                                                      selectInput("ggplot_theme", label = "Theme", choices = list("Grid, white bg"="bw",
                                                                                                                                  "No grid, white bg"="classic",
                                                                                                                                  "Grid, gray bg"="gray",
                                                                                                                                  "Minimal"="min",
                                                                                                                                  "Grid, black bg"="dark",
                                                                                                                                  "Grid, white bg, gray axes"="light",
                                                                                                                                  "Line drawing"="line"),
                                                                                  selected = getOptions("user_options.txt")$gtheme),
                                                                      fluidRow(plotOutput("ggplot_theme_example",inline = F, width="100%")),
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
                                                                                                                                       "Red - yellow - white"="ryw",
                                                                                                                                       "Grayscale"="bw",
                                                                                                                                       "Blues (brew)" = "Blues", 
                                                                                                                                       "Blue - green (brew)" = "BuGn", 
                                                                                                                                       "Blue - purple (brew)" = "BuPu", 
                                                                                                                                       "Green - blue (brew)" = "GnBu", 
                                                                                                                                       "Greens (brew)" = "Greens", 
                                                                                                                                       "Grayscale (brew)" = "Greys",
                                                                                                                                       "Oranges (brew)" = "Oranges", 
                                                                                                                                       "Orange - red (brew)" = "OrRd", 
                                                                                                                                       "Purple - blue (brew)" = "PuBu", 
                                                                                                                                       "Purple - blue - green (brew)" = "PuBuGn",
                                                                                                                                       "Purple - red (brew)" = "PuRd", 
                                                                                                                                       "Purples (brew)" = "Purples", 
                                                                                                                                       "Red - purple (brew)" = "RdPu", 
                                                                                                                                       "Reds (brew)" = "Reds", 
                                                                                                                                       "Yellow - green (brew)" = "YlGn", 
                                                                                                                                       "Yellow - green - blue (brew)" = "YlGnBu", 
                                                                                                                                       "Yellow - orange - brown (brew)" = "YlOrBr", 
                                                                                                                                       "Yellow - orange - red (brew)"="YlOrRd", 
                                                                                                                                       "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", #TODO: add descriptions (or remove all?)
                                                                                                                                       "RdGy", "RdYlBu", "RdYlGn", "Spectral", 
                                                                                                                                       "Accent", "Dark2", "Paired", "Pastel1", 
                                                                                                                                       "Pastel2", "Set1", "Set2", "Set3"),selected = getOptions("user_options.txt")$gspec
                                                                      ),
                                                                      fluidRow(plotly::plotlyOutput("ramp_plot",inline = T, width="100%")),
                                                                      h2("Discrete data"),
                                                                      uiOutput("colorPickers")
                                                             ))
                                    )
                    )),
           tabPanel("",  icon = icon("cog"), value="options", 
                    navbarPage(inverse=TRUE,"Settings", id="tab_settings",
                               tabPanel("Project", icon=icon("gift"),
                                        #textInput(inputId="proj_name", label="Project name", value = ''),
                                        selectizeInput(inputId="proj_name", 
                                                       label="Project name", 
                                                       choices=global$vectors$project_names, 
                                                       selected = getOptions("user_options.txt")$proj_name,
                                                       options=list(create = TRUE)),
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
                                        textOutput("curr_exp_dir")
                               ),
                               tabPanel("Adducts", icon=icon("plus-square"),
                                        h3("Current adduct table:"),
                                        rhandsontable::rHandsontableOutput("adduct_tab", width=800, height=600),
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
                               ),
                               tabPanel("Aesthetic", icon=icon("child"),
                                        h3("Change app settings"),
                                        hr(),
                                        h2("Navigation bar colours"),
                                        colourpicker::colourInput(inputId = "bar.col.1", 
                                                                  label = paste("Active background"), 
                                                                  value = options$col1,
                                                                  allowTransparent = FALSE),
                                        colourpicker::colourInput(inputId = "bar.col.2", 
                                                                  label = paste("Inactive background"), 
                                                                  value = options$col2,
                                                                  allowTransparent = FALSE),
                                        colourpicker::colourInput(inputId = "bar.col.3", 
                                                                  label = paste("Active tab"), 
                                                                  value = options$col3,
                                                                  allowTransparent = FALSE),
                                        colourpicker::colourInput(inputId = "bar.col.4", 
                                                                  label = paste("Inactive tab"), 
                                                                  value = options$col4,
                                                                  allowTransparent = FALSE),
                                        br(),
                                        h2("Fonts (Google fonts)"),
                                        textInput(inputId="font.1", label="h1", value = options$font1),
                                        textInput(inputId="font.2", label="h2", value = options$font2),
                                        textInput(inputId="font.3", label="h3", value = options$font3),
                                        textInput(inputId="font.4", label="body", value = options$font4),
                                        br(), # TODO: font size modifier slider
                                        h2("Font size"),
                                        sliderInput("size.1", label="h1", value=as.numeric(options$size1),min = 5, max=50),
                                        sliderInput("size.2", label="h2", value=as.numeric(options$size2),min = 5, max=50),
                                        sliderInput("size.3", label="h3", value=as.numeric(options$size3),min = 5, max=50),
                                        sliderInput("size.4", label="body", value=as.numeric(options$size4),min = 5, max=50),
                                        br(),
                                        h3("Taskbar image"),
                                        div(imageOutput("taskbar_image",inline = T)),
                                        shinyFilesButton('taskbar_image_path', 
                                                         'Select image', 
                                                         'Please select an image file', 
                                                         FALSE), 
                                        hr(),
                                        # - - - -
                                        actionButton("change_css", "Save settings (restart to apply)")
                               )
                    )
           ), 
           #actionButton("test", "test"),
           tabPanel(title = "Quit", value="stop", icon = icon("times-circle")),
           div(class="spinnylocation1",
               div(class="plus", img(class="imagetop", src=getOptions("user_options.txt")$taskbar_image, width="120px", height="120px")),
               div(class="minus", img(class="imagebottom", src=getOptions("user_options.txt")$taskbar_image, width="120px", height="120px"))
           ),
           div(class="line")
)
)
)