shinyServer(function(input, output, session) {
  
  # ================================= DEFAULTS ===================================
  
  tbl <<- NA
  tables <- list()
  db_search_list <- c()
  theme_update(plot.title = element_text(hjust = 0.5))
  color.function <- rainbow
  patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
  csv_loc <<- file.path(options$work_dir, paste0(options$proj_name, ".csv"))
  shinyOptions(progress.style="old")
  ppm <<- 2
  cf <<- rainbow
  nvars <<- 2
  ncores <<- parallel::detectCores() - 1
  session_cl <<- parallel::makeCluster(ncores)
  
  # -----------------------------------------------
  
  spinnyimg <- reactiveVal("www/electron.png")
  
  output$spinny <- renderText({spinnyimg()})
  
  output$adductSettings <- renderUI({
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
                            sardine(fadeImageButton("add_wikidata", img.path = "wikidata.png"))
                          )),
                 tabPanel(icon("id-card-o"), value = "identifier",
                          br(),
                          tags$i("Select identifier"),
                          radioButtons(inputId = "group_by", label = NULL, choices = 
                                         list(#"Pathway ID" = "pathway",
                                           #"Database ID" = "identifier", 
                                           #"Compound name" = "compoundname",
                                           "Molecular formula" = "baseformula",
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
    )
  })
  
  packages <<- c("data.table", "DBI", "RSQLite", "ggplot2", "minval", "enviPat",
                 "plotly", "parallel", "shinyFiles", "curl", "httr", "pbapply", 
                 "sqldf", "plyr", "ChemmineR", "gsubfn", "stringr", "plotly", "heatmaply",
                 "reshape2", "XML", "xlsx", "colourpicker", "DT","Rserve", "ellipse", 
                 "scatterplot3d","pls", "caret", "lattice", "compiler",
                 "Cairo", "randomForest", "e1071","gplots", "som", "xtable",
                 "RColorBrewer", "xcms","impute", "pcaMethods","siggenes",
                 "globaltest", "GlobalAncova", "Rgraphviz","KEGGgraph",
                 "preprocessCore", "genefilter", "pheatmap", "igraph",
                 "RJSONIO", "SSPA", "caTools", "ROCR", "pROC", "sva", "rJava",
                 "colorRamps", "grDevices", "KEGGREST", "manhattanly")
  
  options <- getOptions(".conf")
  
  default.text <- list(list(name='work_dir',text=options$work_dir),
                       list(name='curr_db_dir',text=options$db_dir),
                       list(name='ppm',text=options$ppm),
                       list(name='analUI',text="Please choose an analysis mode!"),
                       list(name='proj_name',text=options$proj_name),
                       list(name="curr_cpd", text="...")
  )
  
  db_list <<- c("internal", 
                "noise", 
                "hmdb", 
                "chebi", 
                #"pubchem", 
                "kegg",
                "metacyc",
                "wikipathways"
                ,"smpdb",
                "dimedb",
                "wikidata"
  )
  
  images <<- list(list(name = 'cute_package', path = 'www/new-product.png', dimensions = c(80, 80)),
                  list(name = 'umc_logo_int', path = 'www/umcinternal.png', dimensions = c(120, 120)),
                  list(name = 'umc_logo_noise', path = 'www/umcnoise.png', dimensions = c(120, 120)),
                  list(name = 'hmdb_logo', path = 'www/hmdblogo.png', dimensions = c(150, 100)),
                  list(name = 'metacyc_logo', path = 'www/metacyc.png', dimensions = c(300, 80)),
                  list(name = 'chebi_logo', path = 'www/chebilogo.png', dimensions = c(140, 140)),
                  list(name = 'wikipath_logo', path = 'www/wikipathways.png', dimensions = c(130, 150)),
                  list(name = 'kegg_logo', path = 'www/kegglogo.gif', dimensions = c(200, 150)),
                  list(name = 'pubchem_logo', path = 'www/pubchemlogo.png', dimensions = c(145, 90)),
                  list(name = 'smpdb_logo', path = 'www/smpdb_logo_adj.png', dimensions = c(200, 160)),
                  list(name = 'dimedb_logo', path = 'www/dimedb_logo.png', dimensions = c(310, 120)),
                  list(name = 'wikidata_logo', path = 'www/wikidata.png', dimensions = c(250, 200)),
                  list(name = 'pos_icon', path = 'www/handpos.png', dimensions = c(120, 120)),
                  list(name = 'neg_icon', path = 'www/handneg.png', dimensions = c(120, 120)),
                  list(name = 'excel_icon', path = 'www/excel.png', dimensions = c(120, 120)),
                  list(name = 'db_icon', path = 'www/office.png', dimensions = c(100, 100)),
                  list(name = 'csv_icon', path = 'www/office.png', dimensions = c(100, 100)),
                  list(name = 'dataset_icon', path = 'www/office.png', dimensions = c(100, 100))
  )
  
  options$db_dir <<- options$db_dir
  
  # --- render text ---
  
  lapply(default.text, FUN=function(default){
    output[[default$name]] = renderText(default$text)
  })
  
  
  stat.ui.bivar <- reactive({
    navbarPage(inverse=F,h2("Standard analysis"), id="tab_stat",
               tabPanel(h3("PCA"), value = "pca", #icon=icon("cube"),
                        fluidRow(column(12,align="center",plotly::plotlyOutput("plot_pca",height = "600px", width="600px"))),
                        fluidRow(column(3,
                                        selectInput("pca_x", label = "X axis:", choices = paste0("PC",1:30),selected = "PC1",width="100%"),
                                        selectInput("pca_y", label = "Y axis:", choices = paste0("PC",1:30),selected = "PC2",width="100%"),
                                        selectInput("pca_z", label = "Z axis:", choices = paste0("PC",1:30),selected = "PC3",width="100%")
                        ), 
                        column(9, 
                               div(DT::dataTableOutput('pca_tab',width="100%"),style='font-size:80%')
                        )
                        )
               ),
               tabPanel(h3("PLSDA"), value = "plsda",
                        fluidRow(
                          selectInput("plsda_type", label="PLSDA subtype:", choices=list("Normal"="normal",
                                                                                         "Orthogonal"="ortho",
                                                                                         "Sparse"="sparse"), selected=1),
                          actionButton("do_plsda",label="Go")
                        ),
                        hr(),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("globe"),
                                            plotly::plotlyOutput("plot_plsda_3d",height = "600px", width="600px")
                                            # ,fluidRow(column(3,
                                            #                 selectInput("plsda_x", label = "X axis:", choices = paste0("PC",1:3),selected = "PC1",width="100%"),
                                            #                 selectInput("plsda_y", label = "Y axis:", choices = paste0("PC",1:3),selected = "PC2",width="100%"),
                                            #                 selectInput("plsda_z", label = "Z axis:", choices = paste0("PC",1:3),selected = "PC3",width="100%")
                                            # )
                                            # ,
                                            # column(9,
                                            #        div(DT::dataTableOutput('plsda_tab',width="100%"),style='font-size:80%')
                                            # ))
                                   ),
                                   tabPanel("", icon=icon("area-chart"), 
                                            plotOutput("plsda_cv_plot"),
                                            plotOutput("plsda_perm_plot")
                                            
                                   ),
                                   tabPanel("", icon=icon("star-o"), 
                                            plotly::plotlyOutput("plsda_vip_specific_plot",width="100%"),
                                            fluidRow(column(3,
                                                            selectInput("plsda_vip_cmp", label = "Compounds from:", choices = paste0("PC",1:5),selected = "PC1",width="100%")
                                            ), 
                                            column(9, 
                                                   div(DT::dataTableOutput('plsda_vip_tab',width="100%"),style='font-size:80%')
                                            ))
                                            
                                   )
                        )),
               tabPanel(h3("T-test"), value="tt", 
                        fluidRow(plotly::plotlyOutput('tt_specific_plot',width="100%")),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('tt_tab',width="100%"),style='font-size:80%'))
                                   ,tabPanel("", icon=icon("area-chart"),
                                             plotly::plotlyOutput('tt_overview_plot',height="300px")
                                   )
                        )),
               tabPanel(h3("Fold-change"), value="fc",
                        fluidRow(plotly::plotlyOutput('fc_specific_plot',width="100%")),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('fc_tab',width="100%"),style='font-size:80%'))
                                   ,tabPanel("", icon=icon("area-chart"),
                                             plotly::plotlyOutput('fc_overview_plot',height="300px")
                                   )
                        )),
               # =================================================================================
               tabPanel(h3("RandomForest"), value="rf",
                        fluidRow(plotly::plotlyOutput('rf_specific_plot',width="100%")),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('rf_tab',width="100%"),style='font-size:80%'))
                        )
               ),
               tabPanel(h3("Heatmap"), value="heatmap_biv",
                        plotly::plotlyOutput("heatmap",width="110%",height="700px"),
                        br(),
                        fluidRow(column(align="center",
                                        width=12,switchButton(inputId = "heatmode",
                                                              label = "Use data from:", 
                                                              value = TRUE, col = "BW", type = "TTFC"))
                        ),fluidRow(plotly::plotlyOutput('heatmap_specific_plot'),height="200px")
               ),
               tabPanel(h3("Volcano"), value="volc",
                        fluidRow(plotly::plotlyOutput('volc_plot',width="100%",height="600px")),
                        fluidRow(plotly::plotlyOutput('volc_specific_plot'),height="200px")
               ),
               tabPanel(h3("Enrichment"), value="enrich_biv",
                        sidebarLayout(position="left",
                                      sidebarPanel = sidebarPanel(align="center",
                                                                  fluidRow(
                                                                    tags$b("Pathway DBs"),br(),
                                                                    sardine(fadeImageButton("enrich_metacyc", 
                                                                                            img.path = "metacyc.png", 
                                                                                            value = T)),
                                                                    sardine(fadeImageButton("enrich_wikipathways", 
                                                                                            img.path = "wikipathways.png", 
                                                                                            value = T)),
                                                                    sardine(fadeImageButton("enrich_kegg", 
                                                                                            img.path = "kegglogo.gif", 
                                                                                            value = T))
                                                                  ),
                                                                  selectInput('enrich_stats', 
                                                                              'Score source',
                                                                              choices = c("T-test"="tt", 
                                                                                          "Fold-change"="fc",
                                                                                          "PLS-DA"="plsda",
                                                                                          "RandomForest"="rf")
                                                                  ),
                                                                  selectInput('enrich_vals', 
                                                                              'Score threshold',
                                                                              choices = c("Significant"="sig", 
                                                                                          "Top 50"="t50", 
                                                                                          "Top 100"="t100", 
                                                                                          "Top 200"="t200", 
                                                                                          "Top 500"="t500")),
                                                                  hr(),
                                                                  actionButton("go_enrich", "Analyse", icon=icon("binoculars"))
                                      ),
                                      # ------------------
                                      mainPanel = mainPanel(align="center",
                                                            fluidRow(div(DT::dataTableOutput('enriched'),style='font-size:80%')),
                                                            fluidRow(div(DT::dataTableOutput('enrich_pw_tab'),style='font-size:80%')),
                                                            fluidRow(plotly::plotlyOutput('enrich_pw_specific_plot'),height="200px")
                                                            
                                      ))
               )
               ,tabPanel(h3("LasNet"), value = "lasnet",
                         fluidRow(
                           column(width=1, h2("Ridge"),style = "margin-top: 25px;", align="right"),
                           column(width=9,checkboxGroupInput("lasnet_alpha", "Alpha values",
                                                             choiceNames = 
                                                               as.character(seq(0,1,0.1)),
                                                             choiceValues =
                                                               seq(0,1,0.1),
                                                             selected = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                                             inline = TRUE), align="center"
                           ),
                           column(width=1, h2("Lasso"), style = "margin-top: 25px;")),
                         fluidRow(
                           column(width=5,sliderInput("lasnet_trainfrac", 
                                                      label = "Percentage in training", 
                                                      min = 0,
                                                      max = 100,
                                                      step = 5,
                                                      value = 60, 
                                                      post = "%")),
                           column(width=2, actionButton("do_lasnet",label="Go",width = "50px"),style = "margin-top: 35px;", align="left"),
                           column(width=5,sliderInput("lasnet_folds", 
                                                      label = "Fold validation", 
                                                      min = 0,
                                                      max = 20,
                                                      step = 1,
                                                      value = 5, 
                                                      post = "x")
                           )),
                         
                         hr()
                         ,
                         navbarPage(title="Results",id="lasnet",inverse=F,
                                    tabPanel(title = "ROC",value = "",icon=icon("area-chart"),
                                             plotOutput("lasnet_roc_plot",height = "600px")),
                                    tabPanel("Model",value= "",icon=icon("table"),
                                             plotlyOutput("lasnet_specific_plot"),
                                             uiOutput("lasnet_table_ui"))
                         )
               ))
  })
  
  stat.ui.multivar <- reactive({
    navbarPage(inverse=F,h2("Standard analysis"), id="tab_stat",
               tabPanel("", value = "intro", icon=icon("comment-o"),
                        helpText("Info text here")
               ), # pca_legend
               tabPanel(h3("PCA"), value = "pca", #icon=icon("cube"),
                        fluidRow(column(10, plotly::plotlyOutput("plot_pca",height = "600px")),
                                 column(2, br(),br(), br(),plotly::plotlyOutput("pca_legend",height = "400px"))
                        ),
                        fluidRow(column(3,
                                        selectInput("pca_x", label = "X axis:", choices = paste0("PC",1:30),selected = "PC1",width="100%"),
                                        selectInput("pca_y", label = "Y axis:", choices = paste0("PC",1:30),selected = "PC2",width="100%"),
                                        selectInput("pca_z", label = "Z axis:", choices = paste0("PC",1:30),selected = "PC3",width="100%")
                        ), 
                        column(9, 
                               div(DT::dataTableOutput('pca_tab',width="100%"),style='font-size:80%')
                        )
                        )
               ),
               tabPanel(h3("PLSDA"), value = "plsda",
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("globe"),
                                            plotly::plotlyOutput("plot_plsda",width="100%"),
                                            fluidRow(column(3,
                                                            selectInput("plsda_x", label = "X axis:", choices = paste0("PC",1:8),selected = "PC1",width="100%"),
                                                            selectInput("plsda_y", label = "Y axis:", choices = paste0("PC",1:8),selected = "PC2",width="100%"),
                                                            selectInput("plsda_z", label = "Z axis:", choices = paste0("PC",1:8),selected = "PC3",width="100%")
                                            ), 
                                            column(9, 
                                                   div(DT::dataTableOutput('plsda_tab',width="100%"),style='font-size:80%')
                                            ))),
                                   tabPanel("", icon=icon("star-o"), 
                                            plotly::plotlyOutput("plsda_vip_specific_plot",width="100%"),
                                            fluidRow(column(3,
                                                            selectInput("plsda_vip_cmp", label = "Compounds from:", choices = paste0("PC",1:5),selected = "PC1",width="100%")
                                            ), 
                                            column(9, 
                                                   div(DT::dataTableOutput('plsda_vip_tab',width="100%"),style='font-size:80%')
                                            ))
                                            
                                   )
                        )
                        
               ),
               tabPanel(h3("ANOVA"), value="aov",
                        fluidRow(plotly::plotlyOutput('aov_specific_plot',width="100%")),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('aov_tab',width="100%"),style='font-size:80%'))
                                   ,tabPanel("", icon=icon("area-chart"),
                                             plotly::plotlyOutput('aov_overview_plot',height="300px")
                                   )
                        )),
               # =================================================================================
               tabPanel(h3("RandomForest"), value="rf",
                        fluidRow(plotly::plotlyOutput('rf_specific_plot',width="100%")),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('rf_tab',width="100%"),style='font-size:80%'))
                        )
               ),
               tabPanel(h3("Heatmap"), value="heatmap_mult",
                        plotly::plotlyOutput("heatmap",width="110%",height="700px"),
                        fluidRow(plotly::plotlyOutput('heatmap_specific_plot'),height="200px")
               ),
               tabPanel(h3("Enrichment"), value="enrich_multi",
                        sidebarLayout(position="left",
                                      sidebarPanel = sidebarPanel(align="center",
                                                                  fluidRow(
                                                                    tags$b("Pathway DBs"),br(),
                                                                    sardine(fadeImageButton("enrich_metacyc", 
                                                                                            img.path = "metacyc.png", 
                                                                                            value = T)),
                                                                    sardine(fadeImageButton("enrich_wikipathways", 
                                                                                            img.path = "wikipathways.png", 
                                                                                            value = T)),
                                                                    sardine(fadeImageButton("enrich_kegg", 
                                                                                            img.path = "kegglogo.gif", 
                                                                                            value = T))
                                                                  ),
                                                                  selectInput('enrich_stats', 
                                                                              'Score source',
                                                                              choices = c("ANOVA"="aov",
                                                                                          "PLS-DA"="plsda",
                                                                                          "RandomForest"="rf")
                                                                  ),
                                                                  selectInput('enrich_vals', 
                                                                              'Score threshold',
                                                                              choices = c("Significant"="sig", 
                                                                                          "Top 50"="t50", 
                                                                                          "Top 100"="t100", 
                                                                                          "Top 200"="t200", 
                                                                                          "Top 500"="t500")),
                                                                  hr(),
                                                                  actionButton("go_enrich", "Analyse", icon=icon("binoculars"))
                                      ),
                                      # ------------------
                                      mainPanel = mainPanel(align="center",
                                                            fluidRow(div(DT::dataTableOutput('enriched'),style='font-size:80%')),
                                                            fluidRow(div(DT::dataTableOutput('enrich_pw_tab'),style='font-size:80%')),
                                                            fluidRow(plotly::plotlyOutput('enrich_specific_plot'),height="200px")
                                                            
                                      ))
               )
    )
    
  })
  
  
  time.ui <- reactive({
    navbarPage(inverse=F,h2("Time Series"), id="tab_time",
               tabPanel(h3("iPCA"), value = "ipca", 
                        plotly::plotlyOutput("plot_ipca",height="600px"),
                        selectInput("ipca_factor", label = "Color based on:", choices =list("Time"="facA",
                                                                                            "Experimental group"="facB"),width="100%"),
                        fluidRow(column(3,
                                        selectInput("ipca_x", label = "X axis:", choices = paste0("PC",1:90),selected = "PC1",width="100%"),
                                        selectInput("ipca_y", label = "Y axis:", choices = paste0("PC",1:90),selected = "PC2",width="100%"),
                                        selectInput("ipca_z", label = "Z axis:", choices = paste0("PC",1:90),selected = "PC3",width="100%")
                        ), 
                        column(9, 
                               div(DT::dataTableOutput('ipca_tab',width="100%"),style='font-size:80%')
                        ))
               ),
               # =================================================================================
               tabPanel(h3("MEBA"), value="meba", 
                        fluidRow(plotly::plotlyOutput('meba_specific_plot'),height="600px"),
                        fluidRow(div(DT::dataTableOutput('meba_tab', width="100%"),style='font-size:80%'))
               ),
               # =================================================================================
               tabPanel(h3("ASCA"), value="asca",
                        fluidRow(plotly::plotlyOutput('asca_specific_plot', height="600px")),
                        fluidRow(div(DT::dataTableOutput('asca_tab',width="100%"),style='font-size:80%'))
               )
               # =================================================================================
    )
  })
  
  # -----------------
  
  update.UI <- function(){
    if(!exists("mSet")){
      return(NULL)
    }else{
      nvars <<- length(levels(mSet$dataSet$cls))
      which.ui <- if(nvars > 2) stat.ui.multivar else stat.ui.bivar
      # CHANGE UI
      switch(input$exp_type,
             stat = {
               output$analUI <- renderUI({which.ui()})
               output$exp_opt <- renderUI({})
             },
             time_std = {
               output$analUI <- renderUI({time.ui()})
               output$exp_opt <- renderUI({})
             },
             time_fin = {
               output$analUI <- renderUI({which.ui()})
               output$exp_opt <- renderUI({optUI()})
             },
             time_min = {
               output$analUI <- renderUI({which.ui()})
               output$exp_opt <- renderUI({optUI()})
             })
      # set colour pickers (generate)
      output$colourPickers <- renderUI({color.pickers()})
    }
  }
  
  optUI <- function(){		
    selectInput('your.time', 'What end time do you want to pick?', choices = get_times(patdb))
  }
  
  observeEvent(input$exp_type,{
    req(input$exp_type)
    # ---------------------
    update.UI()
  })
  
  # -----------------
  
  output$pos_add_tab <-DT::renderDataTable({
    # -------------
    DT::datatable(pos_adducts,
                  selection = list(mode = 'multiple', selected = c(1:3, nrow(pos_adducts)), target = 'row'),
                  options = list(pageLength = 5, dom = 'tp'), 
                  rownames = F)
  })
  
  output$neg_add_tab <-DT::renderDataTable({
    # -------------
    DT::datatable(neg_adducts,
                  selection = list(mode = 'multiple', selected = c(1, 2, 14, 15, nrow(neg_adducts)), target = 'row'),
                  options = list(pageLength = 5, dom = 'tp'), 
                  rownames = F)
  })
  
  observeEvent(input$sel_all_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(pos_adducts,
                    selection = list(mode = 'multiple', selected = c(1:nrow(pos_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(neg_adducts,
                    selection = list(mode = 'multiple', selected = c(1:nrow(neg_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  observeEvent(input$sel_no_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(pos_adducts,
                    selection = list(mode = 'multiple', selected = c(0), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(neg_adducts,
                    selection = list(mode = 'multiple', selected = c(0), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  
  observeEvent(input$sel_comm_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(pos_adducts,
                    selection = list(mode = 'multiple', selected = c(1:3, nrow(pos_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(neg_adducts,
                    selection = list(mode = 'multiple', selected = c(1, 2, 14:15, nrow(neg_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  
  
  # ===== PACKAGE LOADING ====
  
  observeEvent(input$nav_general, {
    #load.necessarities(input$nav_general)
  })
  
  # =================================================
  
  lapply(images, FUN=function(image){
    output[[image$name]] <- renderImage({
      filename <- normalizePath(image$path)
      # Return a list containing the filename and alt text
      list(src = filename, 
           width = image$dimensions[1],
           height = image$dimensions[2])
    }, deleteFile = FALSE)
  })
  
  # ========================================
  
  observeEvent(input$nav_general, {
    output$package_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(get.package.table(input$nav_general),
                    selection = 'none',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(10, 20, 30), pageLength = 10)
                    , rownames = F)
    })
  }
  )
  
  observeEvent(input$update_packages, {
    req(packages)
    # ----------
    p_load(char = packages, update = T, character.only = T)
    # ---------- 
    firstRun <<- F
    setOption(".conf", "packages_installed", "Y")
    output$package_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 70,
           height = 70)
    }, deleteFile = FALSE)
    # --- restart ---??
    stopApp()
    runApp(".")
  })
  
  
  if(options$packages_installed == "N") return(NULL) # BREAK!!
  
  observe({  
    shinyDirChoose(input, "get_db_dir", roots = c(home = '~'), session = session)
    if(is.null(input$get_db_dir)) return()
    dir <- reactive(input$get_db_dir)
    home <- normalizePath("~")
    given_dir <- file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
    if(is.null(given_dir)) return()
    options$db_dir <<- given_dir
    output$curr_db_dir <- renderText(options$db_dir)
    # --- connect ---
    setOption(".conf", "db_dir", options$db_dir)
    options <<- getOptions(".conf")
  })
  
  observe({  
    shinyDirChoose(input, "get_work_dir", roots = c(home = '~'), session = session)
    if(is.null(input$get_work_dir)) return()
    dir <- reactive(input$get_work_dir)
    home <- normalizePath("~")
    given_dir <- file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
    if(is.null(given_dir)) return()
    options$work_dir <<- given_dir
    output$exp_dir <- renderText(options$work_dir)
    patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
    if(!dir.exists(options$work_dir)) dir.create(options$work_dir)
    # --- connect ---
    setOption(".conf", "work_dir", options$work_dir)
    options <<- getOptions(".conf")
  })
  # ----------------------------------
  
  observeEvent(input$set_proj_name, {
    proj_name <<- input$proj_name
    patdb <<- file.path(options$work_dir, paste0(proj_name,".db", sep=""))
    output$proj_name <<- renderText(proj_name)
    # --- connect ---
    setOption(".conf", "proj_name", proj_name)
    options <<- getOptions(".conf")
  })
  
  observeEvent(input$set_ppm, {
    ppm <<- input$ppm
    output$ppm <<- renderText(ppm)
    # --- connect ---
    setOption(".conf", "ppm", ppm)
    options <<- getOptions(".conf")
  })
  
  observeEvent(input$color_ramp,{
    output$ramp_plot <- plotly::renderPlotly({
      color.function <<- switch(input$color_ramp,
                                "rb"=rainbow,
                                "y2b"=ygobb,
                                "ml1"=matlab.like2,
                                "ml2"=matlab.like,
                                "m2g"=magenta2green,
                                "c2y"=cyan2yellow,
                                "b2y"=blue2yellow,
                                "g2r"=green2red,
                                "b2g"=blue2green,
                                "b2r"=blue2red,
                                "b2p"=cm.colors,
                                "bgy"=topo.colors,
                                "gyw"=terrain.colors,
                                "ryw"=heat.colors,
                                "bw"=blackwhite.colors)
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        titlefont = list(size = 20)
      )
      # --- render ---
      plotly::plot_ly(z = volcano, 
                      colors = color.function(100), 
                      type = "heatmap",
                      showscale=FALSE)  %>%
        layout(xaxis = ax, yaxis = ax)
    })
  })
  
  observeEvent(input$ggplot_theme,{
    plot.theme <<- switch(input$ggplot_theme,
                          bw=ggplot2::theme_bw,
                          classic=ggplot2::theme_classic,
                          gray=ggplot2::theme_gray,
                          min=ggplot2::theme_minimal,
                          dark=ggplot2::theme_dark,
                          light=ggplot2::theme_light,
                          line=ggplot2::theme_linedraw)
  })
  
  observeEvent(input$change_css, {
    # input <- list(bar.col.1 = "ffc8c2", 
    #               bar.col.2 = "e1897f", 
    #               bar.col.3 = "ffffff", 
    #               bar.col.4 = "ffffff",
    #               font.1 = "Press Start 2P",
    #               font.2 = "Raleway bold",
    #               font.3 = "Raleway",
    #               font.4 = "Raleway serif")
    # --- connect ---
    setOption(".conf", "col1", input$bar.col.1)
    setOption(".conf", "col2", input$bar.col.2)
    setOption(".conf", "col3", input$bar.col.3)
    setOption(".conf", "col4", input$bar.col.4)
    
    setOption(".conf", "font1", input$font.1)
    setOption(".conf", "font2", input$font.2)
    setOption(".conf", "font3", input$font.3)
    setOption(".conf", "font4", input$font.4)
  })
  
  color.pickers <- reactive({
    req(mSet$dataSet)
    # -------------
    switch(input$exp_type,
           stat = {
             facs <- levels(mSet$dataSet$cls)
           },
           time_std = {
             lbl.fac <- if(mSet$dataSet$facA.lbl == "Time") "facB" else "facA"
             facs <- levels(mSet$dataSet[[lbl.fac]])
           },
           time_fin = {
             facs <- levels(mSet$dataSet$cls)
           },
           time_min = {           
             facs <- levels(mSet$dataSet$cls)
           })
    default.colours <- rainbow(length(facs))
    # -------------
    lapply(seq_along(facs), function(i) {
      colourpicker::colourInput(inputId = paste("col", i, sep="_"), 
                                label = paste("Choose colour for", facs[i]), 
                                value = default.colours[i],
                                allowTransparent = T) 
    })
  })
  
  color.vec <- reactive({
    req(mSet$dataSet)
    # ----------
    switch(input$exp_type,
           stat = {
             facs <- mSet$dataSet$cls
           },
           time_std = {
             if("facA" %not in% names(mSet$dataSet) & "facB" %not in% names(mSet$dataSet)) return(c("Blue", "Pink"))
             lbl.fac <- if(mSet$dataSet$facA.lbl == "Time") "facB" else "facA"
             facs <- mSet$dataSet[[lbl.fac]]
           },
           time_fin = {
             facs <- mSet$dataSet$cls
           },
           time_min = {           
             facs <- mSet$dataSet$cls
           })
    default.colours <- rainbow(length(facs))
    # -------------
    unlist(lapply(seq_along(facs), function(i) {
      input[[paste("col", i, sep="_")]]
    }))
  })
  
  # --- adduct table editing ---
  
  values = reactiveValues()
  
  observeEvent(input$import_adducts, {
    req(input$add_tab)
    # ----------------
    DF = fread(input$add_tab$datapath)
    output$adduct_tab <- rhandsontable::renderRHandsontable({
      if (!is.null(DF))
        rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
    })
    output$adduct_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  
  adduct_tab_data <- reactive({
    if (!is.null(input$adduct_tab)) {
      DF = hot_to_r(input$adduct_tab)
    } else {
      if (is.null(values[["DF"]]))
        DF = adducts
      else
        DF = values[["DF"]]
    }
    values[["DF"]] = DF
    # ---------------
    DF
  })
  
  output$adduct_tab <- rhandsontable::renderRHandsontable({
    DF = adduct_tab_data()
    if (!is.null(DF))
      rhandsontable::rhandsontable(DF, stretchH = "all", useTypes = TRUE)
  })
  
  observe({
    DF = adduct_tab_data()
    shinyFileSave(input, "save_adducts", roots = c(home = '~'), session=session)
    fileinfo <- parseSavePath(roots = c(home = '~'), input$save_adducts)
    if (nrow(fileinfo) > 0) {
      switch(fileinfo$type,
             csv = fwrite(file = fileinfo$datapath, x = DF)
      )}
  })
  
  # ======================== DB CHECK ============================
  
  # --- check for db files ---
  
  lapply(db_list, FUN=function(db){
    observeEvent(input[[paste0("check_", db)]],{
      db_folder_files <- list.files(options$db_dir)
      is.present <- paste0(db, ".full.db") %in% db_folder_files
      check_pic <- if(is.present) "yes.png" else "no.png"
      output[[paste0(db,"_check")]] <- renderImage({
        # When input$n is 3, filename is ./images/image3.jpeg
        filename <- normalizePath(file.path('www', check_pic))
        # Return a list containing the filename and alt text
        list(src = filename, width = 70,
             height = 70)
      }, deleteFile = FALSE)
    })
  })
  
  # --- build db ---
  
  lapply(db_list, FUN=function(db){
    observeEvent(input[[paste0("build_", db)]], {
      # ---------------------------
      library(RCurl)
      library(XML)
      library(SPARQL)
      # ---------------------------
      withProgress({
        
        parallel::clusterEvalQ(session_cl, library(data.table))
        parallel::clusterEvalQ(session_cl, library(enviPat))
        parallel::clusterEvalQ(session_cl, library(KEGGREST))
        parallel::clusterEvalQ(session_cl, library(RCurl))
        parallel::clusterEvalQ(session_cl, library(SPARQL))
        parallel::clusterEvalQ(session_cl, library(XML))
        
        parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
          "isotopes",
          "subform.joanna", 
          "mergeform.joanna",
          "multiform.joanna",
          "check.ded.joanna",
          #"data.table",
          #"rbindlist",
          #"isopattern",
          #"keggFind",
          #"keggGet",
          "kegg.charge",
          #,"regexpr",
          #"regmatches",
          "xmlParse",
          "getURL"
        ))
        
        # pkgs = c("data.table", "enviPat", "KEGGREST", "XML", "SPARQL", "RCurl")
        # parallel::clusterCall(session_cl, function(pkgs) {
        #   for (req in pkgs) {
        #     require(req, character.only = TRUE)
        #   }
        # }, pkgs = pkgs)
        #setProgress(message = "Working...")
        #build.base.db(db, outfolder=options$db_dir, cl=session_cl)
        setProgress(0.5)
        build.extended.db(db, 
                          outfolder=options$db_dir,
                          adduct.table = adducts, 
                          cl=session_cl, 
                          fetch.limit=200)
      })
    })
  })
  
  # ================== DATA IMPORT ===========================
  
  observeEvent(input$create_db,{
    # --------------------
    patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
    # --------------------
    withProgress({
      setProgress(.1,message = "Loading outlists into memory...")
      req(input$outlist_neg, input$outlist_pos, input$excel)
      
      build.pat.db(patdb,
                   ppm = ppm,
                   pospath = input$outlist_pos$datapath,
                   negpath = input$outlist_neg$datapath,
                   overwrite = T)
      
      setProgress(.95,message = "Adding excel sheets to database...")
      
      exp_vars <<- load.excel(input$excel$datapath, patdb)})
  })
  
  observeEvent(input$import_db, {
    req(input$pat_db)
    patdb <<- input$pat_db$datapath
    output$db_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  
  observeEvent(input$import_csv, {
    req(input$pat_csv)
    csv_loc <<- input$pat_csv$datapath
    output$csv_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
  })
  
  # ==================== CREATE CSV =======================
  
  observeEvent(input$create_csv, {
    req(options$proj_name)
    req(options$work_dir)
    # ---------
    withProgress({
      setProgress(1/4)
      # create csv
      # input = list(group_by = "mz", exp_type = "stat")
      #print(if(input$broadvars) "individual_data" else "setup")
      tbl <- get.csv(patdb,
                     time.series = if(input$exp_type == "time_std") T else F,
                     #exp.condition = input$exp_var,
                     group_adducts = if(length(add_search_list) == 0) F else T,
                     group_by = input$group_by,
                     which_dbs = add_search_list,
                     which_adducts = selected_adduct_list
                     #,var_table = if(input$broadvars) "individual_data" else "setup",
                     #batches = input$batch_var
      )
      
      # --- experimental stuff ---
      switch(input$exp_type,
             stat = {
               tbl.adj <- tbl[,-"Time"]
             },
             time_std = {
               tbl.adj <- tbl
             },
             time_fin = {
               tbl.adj <<- unique(tbl[Time == input$your.time, -"Time"])
               mainmode <<- "stat"
               tbl.adj$Sample <- gsub(tbl.adj$Sample,pattern = "_T.*", replacement = "")
             },
             time_min = {           
               uniq.samples <- unique(tbl$Sample)
               table.base <- tbl[Time == input$your.time,c(1,3)][on=uniq.samples]
               time.end <- tbl[Time == input$your.time,][,4:ncol(tbl),on=uniq.samples]
               time.begin <- tbl[Time == min(as.numeric(tbl$Time)),][,4:ncol(tbl),on=uniq.samples]
               tbl.adj <<- cbind(table.base, 
                                 time.end - time.begin) 
               tbl.adj$Sample <- gsub(tbl.adj$Sample,pattern = "_T.*", replacement = "")
             }
      )
      # -------------------------
      # save csv
      setProgress(2/4)
      csv_loc <<- file.path(options$work_dir, paste0(options$proj_name,".csv"))
      fwrite(tbl.adj, csv_loc, sep="\t")
      # --- overview table ---
      nvars <- length(tbl.adj[,1:(which(colnames(tbl.adj) == "$"))])
      
      setProgress(3/4)
      
      output$csv_tab <-DT::renderDataTable({
        overview_tab <- if(input$exp_type == "time_std"){
          t(data.table(keep.rownames = F,
                       Identifiers = ncol(tbl.adj) - nvars,
                       Samples = nrow(tbl.adj),
                       Times = length(unique(tbl.adj$Time))
          ))
        } else{
          t(data.table(keep.rownames = F,
                       Identifiers = ncol(tbl.adj) - nvars,
                       Samples = nrow(tbl.adj)
          )) 
        }
        # --- render ---
        DT::datatable(overview_tab, 
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE))
      })
    })
    update.UI()
  })
  
  # load previous dataset
  
  observeEvent(input$import_dataset, {
    req(input$pat_dataset)
    # -----------------------------------
    data_loc <<- input$pat_dataset$datapath
    output$dataset_upload_check <- renderImage({
      # When input$n is 3, filename is ./images/image3.jpeg
      filename <- normalizePath('www/yes.png')
      # Return a list containing the filename and alt text
      list(src = filename, width = 20,
           height = 20)
    }, deleteFile = FALSE)
    # -----------
    load(data_loc, envir = .GlobalEnv)
    # -----------
    output$var_norm_plot <- renderPlot(PlotNormSummary())
    output$samp_norm_plot <- renderPlot(PlotSampleNormSummary())
    # ===========
    switch(mSet$dataSet$shinymode,
           stat = {
             output$exp_opt <- renderUI({optUI()})
             output$analUI <- renderUI({stat.ui()})
           },
           time = {
             output$exp_opt <- renderUI({})
             output$analUI <- renderUI({time.ui()})
           })
    output$colourPickers <- renderUI({color.pickers()})
  })
  
  # ===================== METABOSTART ========================
  
  observeEvent(input$check_csv, {
    # ----------------------
    
    csv <- fread(csv_loc, 
                 header = T)
    
    nvars <- length(csv[,1:(which(colnames(csv) == "$"))])
    
    opts <<- colnames(csv[,1:(nvars - 1)])
    
    batch <<- which(sapply(1:(nvars - 1), function(x) length(unique(csv[,..x][[1]])) < nrow(csv)))
    
    numi <<- which(sapply(1:(nvars - 1), function(x) is.numeric(csv[,..x][[1]])))
    
    # get excel table stuff.
    updateSelectInput(session, "exp_var",
                      choices = opts[batch])
    updateSelectInput(session, "samp_var",
                      choices = opts[numi])
    updateSelectizeInput(session, "batch_var",
                         choices = opts[batch],
                         options = list(maxItems = 3L - (length(input$batch_var)))
    )
  })
  
  # observeEvent(input$batch_var, {
  #   maxI = ifelse(any(grepl(pattern = "Batch", x = input$batch_var)), 3L, 2L)
  #   print(input$batch_var)
  #   print(maxI)
  #   updateSelectizeInput(session, "batch_var",
  #                        choices = opts[batch],
  #                        selected = input$batch_var,
  #                        options = list(maxItems = maxI - (length(input$batch_var))))
  # })
  
  ref.selector <- reactive({
    # -------------
    if(input$norm_type == "ProbNorm" | input$norm_type == "CompNorm"){
      fluidRow(
        hr(),
        selectInput('ref_var', 
                    'What is your reference condition?', 
                    choices = c("")),
        actionButton("check_csv", 
                     "Get options", 
                     icon=icon("search")),
        hr()
      )
    }
  })
  
  observeEvent(input$check_csv, {
    req(csv_loc)
    switch(input$norm_type,
           ProbNorm=updateSelectInput(session, "ref_var",
                                      choices = get_ref_vars(fac = "Label") # please add options for different times later, not difficult
           ),
           CompNorm=updateSelectInput(session, "ref_var",
                                      choices = get_ref_cpds() # please add options for different times later, not difficult
           ))
    # get excel table stuff.
    
  })
  
  output$ref_select <- renderUI({ref.selector()})
  
  observeEvent(input$initialize,{
    req(input$filt_type)
    req(input$norm_type)
    req(input$trans_type)
    req(input$scale_type)
    #req(input$batch_corr)
    # match
    withProgress({
      # get curr values from: input$ exp_type, filt_type, norm_type, scale_type, trans_type (later)
      setProgress(.1)
      # ------------------------------
      #Below is your R command history: 
      mSet <<- switch(input$exp_type,
                      stat = {
                        #  mSet<-
                        InitDataObjects("pktable", 
                                        "stat", 
                                        FALSE)
                      },
                      time_std = {
                        InitDataObjects("pktable", 
                                        "ts", 
                                        FALSE)
                      },
                      time_fin = {
                        InitDataObjects("pktable", 
                                        "stat", 
                                        FALSE)
                      },
                      time_min = {   
                        InitDataObjects("pktable", 
                                        "stat", 
                                        FALSE)
                      }
      )
      
      # print(names(mSet))
      # ------- load and re-save csv --------
      #csv_loc = "/Users/jwolthuis/Analysis/BR/BR_FirstEight.csv"
      #csv_loc = "/Users/jwolthuis/Analysis/SP/Spain1.csv"
      #csv_loc <- "/Users/jwolthuis/Analysis/SP/Brazil1.csv"
      #csv_loc <- "/Users/jwolthuis/Documents/umc/data/Data/BrSpIt/MZXML/results/specpks_grouped_mdq/grouped_pos.csv"
      #csv_loc <- "/Users/jwolthuis/Analysis/SP/BrazilAndSpain.csv"
      #csv_loc <- "~/Analysis/SP/Gilbert_W.csv"
      
      #f = fread(csv_loc, header = TRUE)
      
      csv_orig <- fread(csv_loc, 
                        data.table = TRUE,
                        header = T)
      
      nvars <- length(csv_orig[,1:(which(colnames(csv_orig) == "$"))])
      
      # --- batches ---
      
      batches <- input$batch_var
      condition <- input$exp_var
      
      if(is.null(batches)) batches <- ""
      
      batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE
      
      print(batch_corr)
      
      if("Batch" %in% batches){ batches = c(batches, "Injection") }
      
      first_part <<- csv_orig[,1:(nvars-1),with=FALSE]
      first_part[first_part == "" | is.null(first_part)] <- "Unknown"
      
      csv_temp <- cbind(first_part[,!duplicated(names(first_part)),with=FALSE],
                        "Label" = first_part[,..condition][[1]],
                        "$" = c(0),
                        csv_orig[,-c(1:nvars),with=FALSE])
      
      nvars <- length(csv_temp[,1:(which(colnames(csv_temp) == "$"))])
      
      # --- remove the rest ---
      
      remove = setdiff(colnames(csv_temp)[1:nvars], c(batches, "Label","Batch","Sample","Time"))
      remain = intersect(colnames(csv_temp)[1:nvars], c(batches, "Label","Batch","Sample","Time"))
      
      csv_subset <- csv_temp[, !remove, with=FALSE]
      
      # --- remove outliers? ---
      print("Removing outliers...")
      
      sums <- rowSums(csv_subset[, -c(1:length(remain)),with=FALSE],na.rm = TRUE)
      names(sums) <- csv_subset$Sample
      outliers = c(car::Boxplot(as.data.frame(sums)))
      
      print(paste("Outliers removed (based on total signal):", paste(outliers, collapse=" and ")))
      csv_temp_no_out <- csv_subset[!(Sample %in% outliers)]
      
      batchview = if(condition == "Batch") TRUE else FALSE
      
      if(any(grepl("QC", csv_temp_no_out$Sample))){
        print("hereee")
        batch_table <- csv_temp_no_out[,c("Sample","Batch"),with=FALSE]
        samps <- which(!grepl(csv_temp_no_out$Sample, pattern = "QC"))
        batchnum <- unique(csv_temp_no_out[samps, "Batch"][[1]])
        keep_samps <- batch_table[which(batch_table$Batch %in% batchnum),"Sample"][[1]]
        batch_table <- batch_table[which(batch_table$Batch %in% batchnum)]
        if("Batch" %in% batches){
          csv_temp_filt <- csv_temp_no_out[which(csv_temp_no_out$Sample %in% keep_samps),]
        }else{
          csv_temp_filt <- csv_temp_no_out[which(csv_temp_no_out$Sample %in% keep_samps), -"Batch"]
        }
      }else{
        csv_temp_filt <- csv_temp_no_out[, -"Batch"]
      }
      
      csv_loc_no_out <- gsub(pattern = "\\.csv", replacement = "_no_out.csv", x = csv_loc)
      
      if(batch_corr){
        batch_table <- csv_temp_filt[,c("Sample",batches),with=FALSE]
        fwrite(csv_temp_filt[,-batches, with=FALSE], 
               csv_loc_no_out)  
      }else{
        
        batch_table <- csv_temp_filt[,c("Sample"),with=FALSE]
        fwrite(csv_temp_filt, 
               csv_loc_no_out)  
      }
      
      rownames(batch_table) <- batch_table$Sample
      
      # -------------------------------------
      
      mSet <<- Read.TextData(mSet, 
                             filePath = csv_loc_no_out, 
                             "rowu")
      
      mSet$dataSet$batches <<- batch_table
      
      # ----------- sanity check ------------
      
      mSet <<- SanityCheckData(mSet)
      
      # input = list(perc_limit = 50, miss_type = "colmin", filt_type = "iqr", norm_type = "QuantileNorm", trans_type = "LogNorm", scale_type = "AutoNorm", ref_var=NULL)
      
      mSet <<- RemoveMissingPercent(mSet, 
                                    percent = input$perc_limit/100)
      
      if(input$miss_type != "none"){
        mSet <<- ImputeVar(mSet,
                           method = input$miss_type)
      }
      
      # ------------------------
      
      if(input$filt_type != "none"){
        print("filtering...")
        mSet <<- FilterVariable(mSet,
                                filter = input$filt_type,
                                qcFilter = "F",
                                rsd = 25)
      }
      # ------------------------------------
      
      setProgress(.2)
      
      keep <- intersect(first_part$Sample, rownames(mSet$dataSet$preproc))
      
      first_part_no_out <<- first_part[Sample %in% keep]
      
      if(input$norm_type == "SpecNorm"){
        norm.vec <<- first_part_no_out[match(first_part_no_out$Sample,
                                             rownames(mSet$dataSet$preproc)),][[input$samp_var]]
        norm.vec <<- scale(x = norm.vec,center = 1)
        print(norm.vec)
      }else{
        norm.vec <<- rep(1, length(mSet$dataSet$cls))
      }
      
      mSet <<- Normalization(mSet,
                             rowNorm = input$norm_type,
                             transNorm = input$trans_type,
                             scaleNorm = input$scale_type,
                             ref = input$ref_var
      )
      
      
      second_part_no_out <- mSet$dataSet$norm[match(rownames(mSet$dataSet$norm),
                                                    first_part_no_out$Sample),]
      csv_filled = cbind(first_part_no_out, 
                         second_part_no_out)
      
      csv_loc_filled <- gsub(pattern = "\\.csv", replacement = "_filled.csv", x = csv_loc)
      
      fwrite(csv_filled, file = csv_loc_filled)
      #print(ncol(mSet$dataSet$orig))
      #print(ncol(mSet$dataSet$norm))
      
      # mSet <- Normalization(mSet = mSet,
      #                       rowNorm = "QuantileNorm",
      #                       transNorm = "LogNorm",
      #                       scaleNorm = "AutoNorm"
      #                       # ref = "C15H13N1O2"
      #                       #ratio=FALSE,
      #                       #ratioNum=20
      #                       )
      # ====
      
      setProgress(.4)
      
      #mSet <- mSet_beforebatch
      
      matr <- mSet$dataSet$norm
      
      smps <- rownames(mSet$dataSet$norm)
      
      qc_rows <- which(grepl(pattern = "QC", x = smps))
      if(length(qc_rows) == 0) qc_rows = 1:length(smps)
      
      smpnames = smps[qc_rows]
      # mSet$dataSet$qc <<- mSet$dataSet$norm[qc_rows,]
      
      if(batch_corr){
        
        if("Batch" %in% colnames(mSet$dataSet$batches)){
          
          # library(BatchCorrMetabolomics)
          
          mSet_before_qc <<- mSet
          
          # --- QC CORRECTION ---
          
          print("Correcting batch effect w/ QC...")
          
          corr_cols <- pbapply::pblapply(1:ncol(matr), FUN=function(i){
            vec = matr[,i]
            print(vec)
            corr_vec = BatchCorrMetabolomics::doBC(Xvec = as.numeric(vec), 
                                                   ref.idx = as.numeric(qc_rows), 
                                                   batch.idx = as.numeric(mSet$dataSet$batches[match(smps, rownames(mSet$dataSet$batches)),"Batch"][[1]]),
                                                   seq.idx = as.numeric(mSet$dataSet$batches[match(smps, rownames(mSet$dataSet$batches)),"Injection"][[1]]))
            # ---------
            corr_vec
          })
          
          qc_corr_matrix <- as.data.frame(do.call(cbind, corr_cols))
          
          colnames(qc_corr_matrix) <- colnames(matr)
          rownames(qc_corr_matrix) <- rownames(matr)
          
          mSet$dataSet$norm <<- as.data.frame(qc_corr_matrix)
          
        }
        
        if(!batchview){
          mSet$dataSet$norm <<- mSet$dataSet$norm[!grepl(rownames(mSet$dataSet$norm),pattern= "QC"),]
          mSet$dataSet$cls <<- mSet$dataSet$cls[which(!grepl(rownames(mSet$dataSet$norm),pattern= "QC")), drop = TRUE]
          mSet$dataSet$batches <<- mSet$dataSet$batches[!grepl(mSet$dataSet$batches$Sample,pattern= "QC"),]
          mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
        }
        
        left_batch_vars <- grep(colnames(mSet$dataSet$batches), 
                                pattern =  "Batch|Injection|Sample",
                                value = T,
                                invert = T)
        
        print(left_batch_vars)
        
        if(length(left_batch_vars) > 2){ 
          print("Can only correct for 2 batches...")
        } else if(length(left_batch_vars) == 0){
          print("No vars usable for ComBat...")
          NULL
        } else{
          
          smp <- rownames(mSet$dataSet$norm)
          exp_lbl <- mSet$dataSet$cls
          
          print(paste("Batches:", left_batch_vars))
          
          csv <- as.data.table(cbind(Sample = smp, 
                                     Label = mSet$dataSet$cls,
                                     mSet$dataSet$norm))
          
          csv_edata <-t(csv[,!c(1,2)])
          colnames(csv_edata) <- csv$Sample
          
          if(length(left_batch_vars) == 1){
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$batches[match(smp, mSet$dataSet$batches$Sample),left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = c(0),
                                    outcome = as.factor(exp_lbl))     
            batch_normalized= t(sva::ComBat(dat=csv_edata,
                                            batch=csv_pheno$batch1
                                            #mod=mod.pheno,
                                            #par.prior=TRUE
            ))
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }else{
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$batches[match(smp, mSet$dataSet$batches$Sample), left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = mSet$dataSet$batches[match(smp, mSet$dataSet$batches$Sample), left_batch_vars[2], with=FALSE][[1]],
                                    outcome = as.factor(exp_lbl))
            batch_normalized = t(limma::removeBatchEffect(x = csv_edata,
                                                          #design = mod.pheno,
                                                          batch = csv_pheno$batch1,
                                                          batch2 = csv_pheno$batch2))
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }
          
          print(head(csv_pheno))
          
          #mod.batch = model.matrix(~ batch1 + batch2 + outcome, data=csv_pheno)
          
          mSet$dataSet$norm <<- as.data.frame(batch_normalized)
        }
      } else{
        if(!batchview){
          mSet$dataSet$norm <<- mSet$dataSet$norm[!grepl(rownames(mSet$dataSet$norm),pattern= "QC"),]
          mSet$dataSet$cls <<- mSet$dataSet$cls[which(!grepl(rownames(mSet$dataSet$norm),pattern= "QC")), drop = TRUE]
          mSet$dataSet$batches <<- mSet$dataSet$batches[!grepl(mSet$dataSet$batches$Sample,pattern= "QC"),]
          mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
        }
      }
      
      print(rownames(mSet$dataSet$norm))
      
      
      setProgress(.5)
      
      setProgress(.6)
      varNormPlots <- ggplotNormSummary(mSet)
      output$var1 <- renderPlot(varNormPlots$tl)
      output$var2 <- renderPlot(varNormPlots$bl)
      output$var3 <- renderPlot(varNormPlots$tr)
      output$var4 <- renderPlot(varNormPlots$br)
      
      setProgress(.7)
      sampNormPlots <-  ggplotSampleNormSummary(mSet)
      output$samp1 <- renderPlot(sampNormPlots$tl)
      output$samp2 <- renderPlot(sampNormPlots$bl)
      output$samp3 <- renderPlot(sampNormPlots$tr)
      output$samp4 <- renderPlot(sampNormPlots$br)
      setProgress(.8)
      
      mSet$dataSet$adducts <<- selected_adduct_list
      update.UI()  
      setProgress(.9)
      
    })
  })
  
  # MAIN EXECUTION OF ANALYSES
  
  observeEvent(input$tab_time, {
    if(input$nav_general != "analysis") return(NULL)
    # get excel table stuff.
    switch(input$tab_time,
           ipca = {
             # --------------------
             output$plot_ipca <- plotly::renderPlotly({
               if("ipca" %not in% names(mSet$analSet)){
                 iPCA.Anal(mSet, file.path(options$work_dir, "ipca_3d_0_.json"))
               }
               fac.lvls <- unique(mSet$analSet$ipca$score[[input$ipca_factor]])
               chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(fac.lvls)
               # ---------------
               df <- t(as.data.frame(mSet$analSet$ipca$score$xyz))
               x <- gsub(input$ipca_x,pattern = "\\(.*$", replacement = "")
               y <- gsub(input$ipca_y,pattern = "\\(.*$", replacement = "")
               z <- gsub(input$ipca_z,pattern = "\\(.*$", replacement = "")
               x.num <- as.numeric(gsub(x, pattern = "PC", replacement = ""))
               y.num <- as.numeric(gsub(y, pattern = "PC", replacement = ""))
               z.num <- as.numeric(gsub(z, pattern = "PC", replacement = ""))
               # ---------------
               plots <- plotly::plot_ly()
               for(class in fac.lvls){
                 row = which(mSet$analSet$ipca$score[[input$ipca_factor]] == class)
                 # ---------------------
                 xc=df[row, x.num]
                 yc=df[row, y.num]
                 zc=df[row, z.num]
                 # --- plot ellipse ---
                 o <- rgl::ellipse3d(cov(cbind(xc,yc,zc)), 
                                     centre=c(mean(xc), 
                                              mean(yc), 
                                              mean(zc)), 
                                     level = 0.95)
                 mesh <- c(list(x = o$vb[1, o$ib]/o$vb[4, o$ib], 
                                y = o$vb[2, o$ib]/o$vb[4, o$ib], 
                                z = o$vb[3, o$ib]/o$vb[4, o$ib]))
                 plots <- plots %>% add_trace(
                   x=mesh$x, 
                   y=mesh$y, 
                   z=mesh$z, 
                   type='mesh3d',
                   alphahull=0,
                   opacity=0.1,
                   hoverinfo="none"
                 )
               }
               adj_plot <<- plotly_build(plots)
               rgbcols <- toRGB(chosen.colors)
               c = 1
               for(i in seq_along(adj_plot$x$data)){
                 item = adj_plot$x$data[[i]]
                 if(item$type == "mesh3d"){
                   adj_plot$x$data[[i]]$color <- rgbcols[c]
                   c = c + 1
                 }
               }
               # --- render! ---
               ipca_plot <<- adj_plot %>%
                 add_trace(
                   x = df[,x.num], 
                   y = df[,y.num], 
                   z = df[,z.num], 
                   opacity=1,
                   type = "scatter3d",
                   color= mSet$analSet$ipca$score[[input$ipca_factor]], colors=chosen.colors
                 ) %>%  layout(scene = list(
                   aspectmode="cube",
                   xaxis = list(
                     titlefont = list(size = 20),
                     title = mSet$analSet$ipca$score$axis[x.num]),
                   yaxis = list(
                     titlefont = list(size = 20),
                     title = mSet$analSet$ipca$score$axis[y.num]),
                   zaxis = list(
                     titlefont = list(size = 20),
                     title = mSet$analSet$ipca$score$axis[z.num])))
               # --- return ---
               ipca_plot
             })
             # ------------------
             ipca.table <- as.data.table(data.table(gsub(mSet$analSet$ipca$score$axis, pattern = "\\(.*$| ", replacement = ""),
                                                    gsub(mSet$analSet$ipca$score$axis, pattern = "PC\\d*| |\\(|%\\)", replacement = "")),
                                         keep.rownames = T)
             colnames(ipca.table) <- c("Principal Component", "% variance")
             output$ipca_tab <-DT::renderDataTable({
               req(ipca.table)
               # -------------
               DT::datatable(ipca.table, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
           },
           meba = {
             if("MB" %not in% names(mSet$analSet)){
               performMB(10, dir=options$work_dir)
             }
             output$meba_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$MB$stats, 
                             selection = 'single',
                             colnames = c("Compound", "Hotelling/T2 score"),
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
           },
           manova = {
             NULL
           },
           asca = {
             if("asca" %not in% names(mSet$analSet)){
               Perform.ASCA(1, 1, 2, 2)
               CalculateImpVarCutoff(0.05, 0.9, dir=options$work_dir)
             }
             output$asca_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$asca$sig.list$Model.ab, 
                             selection = 'single',
                             colnames = c("Compound", "Leverage", "SPE"),
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
             
           })
  })
  
  # outliers
  # 
  # outliers <- names(which(abs(mSet$analSet$pca$x[,1]) > 50))
  # # remove
  # csv <- as.data.table(fread(csv_loc))
  # csv_no_out <- csv[!Sample %in% outliers]
  # fwrite(csv_no_out,file = "ayr_no_outliers.csv")
  
  observeEvent(input$tab_stat, {
    if(input$nav_general != "analysis") return(NULL)
    nvars <<- length(levels(mSet$dataSet$cls))
    # get excel table stuff.
    switch(input$tab_stat,
           pca = {
             if(!"pca" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- PCA.Anal(mSet)
               })
             }
             output$plot_pca <- plotly::renderPlotly({
               df <- mSet$analSet$pca$x
               x <- input$pca_x
               y <- input$pca_y
               z <- input$pca_z
               #x =1;y=2;z=3
               x.var <- round(mSet$analSet$pca$variance[x] * 100.00, digits=1)
               y.var <- round(mSet$analSet$pca$variance[y] * 100.00, digits=1)
               z.var <- round(mSet$analSet$pca$variance[z] * 100.00, digits=1)
               fac.lvls <- length(levels(mSet$dataSet$cls))
               
               chosen.colors <<- if(fac.lvls == length(color.vec())) color.vec() else rainbow(fac.lvls)
               # --- add ellipses ---
               classes <- mSet$dataSet$cls
               plots <- plotly::plot_ly(showlegend = F)
               
               for(class in levels(classes)){
                 row = which(classes == class)
                 # ---------------------
                 xc=mSet$analSet$pca$x[row, x]
                 yc=mSet$analSet$pca$x[row, y]
                 zc=mSet$analSet$pca$x[row, z]
                 # --- plot ellipse ---
                 o <- rgl::ellipse3d(cov(cbind(xc,yc,zc)), 
                                     centre=c(mean(xc), 
                                              mean(yc), 
                                              mean(zc)), 
                                     level = 0.95)
                 mesh <- c(list(x = o$vb[1, o$ib]/o$vb[4, o$ib], 
                                y = o$vb[2, o$ib]/o$vb[4, o$ib], 
                                z = o$vb[3, o$ib]/o$vb[4, o$ib]))
                 plots = plots %>% add_trace(
                   x=mesh$x, 
                   y=mesh$y, 
                   z=mesh$z, 
                   type='mesh3d',
                   alphahull=0,
                   opacity=0.1,
                   hoverinfo="none"
                 )
               }
               adj_plot <<- plotly_build(plots)
               rgbcols <- toRGB(chosen.colors)
               c = 1
               for(i in seq_along(adj_plot$x$data)){
                 item = adj_plot$x$data[[i]]
                 if(item$type == "mesh3d"){
                   adj_plot$x$data[[i]]$color <- rgbcols[c]
                   adj_plot$x$data[[i]]$visible <- TRUE
                   c = c + 1
                 }
               }
               # --- return ---
               pca_plot <<- adj_plot %>% add_trace(
                 hoverinfo = 'text',
                 text = rownames(df),
                 x = mSet$analSet$pca$x[,x], 
                 y = mSet$analSet$pca$x[,y], 
                 z = mSet$analSet$pca$x[,z], 
                 visible = rep(T, times=fac.lvls),
                 type = "scatter3d",
                 opacity=1,
                 color= mSet$dataSet$cls, colors=chosen.colors
               ) %>%  layout(scene = list(
                 aspectmode="cube",
                 xaxis = list(
                   titlefont = list(size = 20),
                   title = fn$paste("$x ($x.var %)")),
                 yaxis = list(
                   titlefont = list(size = 20),
                   title = fn$paste("$y ($y.var %)")),
                 zaxis = list(
                   titlefont = list(size = 20),
                   title = fn$paste("$z ($z.var %)")))) 
               # --- return ---
               pca_plot
             })
             # === LEGEND ===
             output$pca_legend <- plotly::renderPlotly({
               frame <- data.table(x = c(1), 
                                   y = fac.lvls)
               p <- ggplot(data=frame,
                           aes(x, 
                               y, 
                               color=factor(y),
                               fill=factor(y))) + 
                 geom_point(shape = 21, size = 5, stroke = 5) +
                 scale_colour_manual(values=chosen.colors) +
                 theme_void() + 
                 theme(legend.position="none")
               # --- return ---
               ggplotly(p, tooltip = NULL) %>% config(displayModeBar = F)
             })
             # ==============
             pca.table <- as.data.table(round(mSet$analSet$pca$variance * 100.00,
                                              digits = 2),
                                        keep.rownames = T)
             colnames(pca.table) <- c("Principal Component", "% variance")
             output$pca_tab <-DT::renderDataTable({
               req(pca.table)
               # -------------
               DT::datatable(pca.table, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })},
           plsdaa = {NULL},
           heatmap_mult = {
             req(input$color_ramp)
             req(mSet$analSet$aov)
             used.analysis <- "aov"
             used.values <- "p.value"
             sourceTable <- sort(mSet$analSet[[used.analysis]][[used.values]])[1:100]
             withProgress({
               output$heatmap <- plotly::renderPlotly({
                 x <- mSet$dataSet$norm[,names(sourceTable)]
                 final_matrix <- t(x)
                 translator <- data.table(Sample=rownames(mSet$dataSet$norm),Group=mSet$dataSet$prenorm.cls)
                 group_assignments <- translator[,"Group",on=colnames(final_matrix)]$Group
                 hm_matrix <<- heatmaply::heatmapr(final_matrix, 
                                                   Colv=T, 
                                                   Rowv=T,
                                                   col_side_colors = group_assignments,
                                                   k_row = NA)
                 heatmaply::heatmaply(hm_matrix,
                                      Colv=F, Rowv=F,
                                      branches_lwd = 0.3,
                                      margins = c(60,0,NA,50),
                                      colors = color.function(256),
                                      col_side_colors = function(n){
                                        if(n == length(color.vec())) color.vec() else rainbow(n)
                                        rainbow(n)
                                      },
                                      subplot_widths = c(.9,.1),
                                      subplot_heights =  c(.1,.05,.85),
                                      column_text_angle = 90)
               })
             })
           },
           heatmap_biv = {
             #req(input$color_ramp)
             req(mSet$analSet$tt)
             req(mSet$analSet$fc)
             req(input$heatmode)
             
             if(input$heatmode){
               used.analysis <- "tt"
               used.values <- "p.value"
             }else{
               used.analysis <- "fc"
               used.values <- "fc.all"
             }
             #sourceTable <- mSet$analSet[[used.analysis]]$inx.imp[mSet$analSet[[used.analysis]]$inx.imp == TRUE]
             sourceTable <- sort(mSet$analSet[[used.analysis]][[used.values]])[1:100]
             # --- render ---
             withProgress({
               output$heatmap <- plotly::renderPlotly({
                 x <- mSet$dataSet$norm[,names(sourceTable)]
                 final_matrix <- t(x)
                 translator <- data.table(Sample=rownames(mSet$dataSet$norm),Group=mSet$dataSet$cls)
                 group_assignments <- translator[,"Group",on=colnames(final_matrix)]$Group
                 hm_matrix <<- heatmaply::heatmapr(final_matrix, 
                                                   Colv=T, 
                                                   Rowv=T,
                                                   col_side_colors = group_assignments,
                                                   k_row = NA)
                 heatmaply::heatmaply(hm_matrix,
                                      Colv=F, Rowv=F,
                                      branches_lwd = 0.3,
                                      margins = c(60,0,NA,50),
                                      colors = color.function(256),
                                      col_side_colors = function(n){
                                        if(n == length(color.vec())) color.vec() else rainbow(n)
                                        rainbow(n)
                                      },
                                      subplot_widths = c(.9,.1),
                                      subplot_heights =  c(.1,.05,.85),
                                      column_text_angle = 90)
               })               
             })
           },
           tt = {
             if(!"tt" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- Ttests.Anal(mSet,
                                      nonpar = FALSE, 
                                      threshp = 0.05, 
                                      paired = FALSE,
                                      equal.var = TRUE
                                      #  multicorr = "BH"
                 )
               })
             }
             res <<- mSet$analSet$tt$sig.mat
             if(is.null(res)) res <<- data.table("No significant hits found")
             output$tt_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
             output$tt_overview_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotTT(color.function, 20)
             })
           },
           fc = {
             if(!"fc" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- FC.Anal.unpaired(mSet,
                                           2.0, 
                                           1)                 
               })
             }
             res <<- mSet$analSet$fc$sig.mat
             if(is.null(res)) res <<- data.table("No significant hits found")
             output$fc_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
             output$fc_overview_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotFC(color.function, 20)
             })
           },
           aov = {
             if(!"aov" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- ANOVA.Anal(mSet, 
                                     thresh=0.05,
                                     nonpar = F)
               })
             }
             res <<- mSet$analSet$aov$sig.mat
             if(is.null(res)) res <<- data.table("No significant hits found")
             output$aov_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
           },
           rf={
             if(!"rf" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- RF.Anal(mSet, 500,7,1)
               })
             }
             vip.score <<- as.data.table(mSet$analSet$rf$importance[, "MeanDecreaseAccuracy"],keep.rownames = T)
             colnames(vip.score) <<- c("rn", "accuracyDrop")
             output$rf_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(vip.score, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
           },
           volc = {
             if(!"volc" %in% names(mSet$analSet)){
               withProgress({
                 mSet <<- Volcano.Anal(mSet,FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
               })
             }
             output$volc_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(mSet$analSet$volc$sig.mat, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
             output$volc_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggPlotVolc(color.function, 20)
             })
           }
    )
  })
  
  observeEvent(input$do_plsda, {
    print("Doing plsda...")
    library(e1071)
    print(input$plsda_type)
    switch(input$plsda_type,
           normal={
             withProgress({
               mSet <<- PLSR.Anal(mSet,
                                  TRUE)
               setProgress(0.3)
               
               mSet <<- PLSDA.CV(mSet,compNum = 3)
               mSet <<- PLSDA.Permut(mSet,num = 100)
               
               output$plsda_cv_plot <- renderPlot({ggPlotClass(cf = cf)})
               output$plsda_perm_plot <- renderPlot({ggPlotPerm(cf = cf)})
               
               print("here")
               # - - - - - - - - - - - 
               setProgress(0.6)
               # - - overview table - -
               plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                                  / mSet$analSet$plsr$Xtotvar 
                                                  * 100.0,
                                                  digits = 2),
                                            keep.rownames = T)
               colnames(plsda.table) <- c("Principal Component", "% variance")
               plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
               setProgress(0.9)
               # --- coordinates ---
               coords <<- mSet$analSet$plsr$scores
               colnames(coords) <<- paste0("PC", 1:ncol(coords))
               # --- vip table ---
               colnames(mSet$analSet$plsda$vip.mat) <- paste0("PC", 1:ncol(mSet$analSet$plsda$vip.mat))
               compounds_pc <- as.data.table(mSet$analSet$plsda$vip.mat,keep.rownames = T)
               ordered_pc <- setorderv(compounds_pc, input$plsda_vip_cmp, -1)
               plsda_tab <<- cbind(ordered_pc[1:50, c("rn")], 
                                   Rank=c(1:50))
             })
           },
           sparse={
             mSet <<- SPLSR.Anal(mSet, 5, 10, "same")
             plsda.table <- data.table("Principal Component" = paste0("PC", 1:length(mSet$analSet$splsr$explained_variance$X)),
                                       "% variance" = round(100 * mSet$analSet$splsr$explained_variance$X, 2))
             coords <- data.frame(signif(mSet$analSet$splsr$variates$X, 5))
             compounds_pc <- as.data.table(abs(mSet$analSet$splsr$loadings$X),keep.rownames = T)
             colnames(compounds_pc)[2:ncol(compounds_pc)] <- paste0("PC", 1:(ncol(compounds_pc)-1))
             ordered_pc <- setorderv(compounds_pc, "PC1", -1)
             plsda_tab <<- cbind(ordered_pc[1:50, c("rn")],
                                 Rank=c(1:50))
             ## TODO ##
           },
           ortho={
             mSet <<- OPLSR.Anal(mSet, reg=TRUE)
             coords <- data.frame(signif(mSet$analSet$oplsda$coefficients, 5))
             coords <- data.frame(mSet$analSet$oplsda$scoreMN,
                                  mSet$analSet$oplsda$orthoScoreMN)
             plsda.table = data.table('Principal Component'=c(paste0("PC", 1:4)))
             compounds_pc <- as.data.table(abs(cbind(mSet$analSet$oplsda$loadingMN,
                                                     mSet$analSet$oplsda$orthoLoadingMN)),
                                           keep.rownames = T)
             colnames(compounds_pc)[2:ncol(compounds_pc)] <- paste0("PC", 1:(ncol(compounds_pc)-1))
             ordered_pc <- setorderv(compounds_pc, "PC1", -1)
             plsda_tab <<- cbind(ordered_pc[1:50, c("rn")],
                                 Rank=c(1:50))
             
             # mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
             # mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "png", 72, width=NA, "all");
             # mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
             # mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
             # mSet<-PlotOPLS.Permutation(mSet, "opls_perm_1_", "png", 72, 100, width=NA)
             # ## TODO ##
             
           })
    # -----------
    output$plot_plsda_3d <- plotly::renderPlotly({
      # x <- input$plsda_x
      # y <- input$plsda_y
      # z <- input$plsda_z
      x <- "PC1"
      y <- "PC2"
      z <- "PC3"
      x.var <- plsda.table[`Principal Component` == x, 
                           `% variance`]
      y.var <- plsda.table[`Principal Component` == y, 
                           `% variance`]
      z.var <- plsda.table[`Principal Component` == z, 
                           `% variance`]
      fac.lvls <- unique(mSet$dataSet$cls)
      chosen.colors <- if(fac.lvls == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
      # --- add ellipses ---
      classes <- mSet$dataSet$cls
      plots <- plotly::plot_ly(showlegend = F)
      
      for(class in levels(classes)){
        row = which(classes == class)
        # ---------------------
        xc=coords[row, x]
        yc=coords[row, y]
        zc=coords[row, z]
        # --- plot ellipse ---
        o <- rgl::ellipse3d(cov(cbind(xc,yc,zc)), 
                            centre=c(mean(xc), 
                                     mean(yc), 
                                     mean(zc)), 
                            level = 0.95)
        mesh <- c(list(x = o$vb[1, o$ib]/o$vb[4, o$ib], 
                       y = o$vb[2, o$ib]/o$vb[4, o$ib], 
                       z = o$vb[3, o$ib]/o$vb[4, o$ib]))
        plots = plots %>% add_trace(
          x=mesh$x, 
          y=mesh$y, 
          z=mesh$z, 
          type='mesh3d',
          alphahull=0,
          opacity=0.1,
          hoverinfo="none"
        )
      }
      adj_plot <<- plotly_build(plots)
      rgbcols <- toRGB(chosen.colors)
      c = 1
      for(i in seq_along(adj_plot$x$data)){
        item = adj_plot$x$data[[i]]
        if(item$type == "mesh3d"){
          adj_plot$x$data[[i]]$color <- rgbcols[c]
          adj_plot$x$data[[i]]$visible <- TRUE
          c = c + 1
        }
      }
      # ---------------
      adj_plot %>%
        add_trace(
          x = coords[,x], 
          y = coords[,y], 
          z = coords[,z], 
          type = "scatter3d",
          color = mSet$dataSet$cls, 
          colors=chosen.colors,
          opacity=1
        ) %>%  layout(scene = list(
          aspectmode="cube",
          xaxis = list(
            titlefont = list(size = 20),
            title = fn$paste("$x ($x.var %)")),
          yaxis = list(
            titlefont = list(size = 20),
            title = fn$paste("$y ($y.var %)")),
          zaxis = list(
            titlefont = list(size = 20),
            title = fn$paste("$z ($z.var %)")
          )
        ))
    })
    output$plsda_tab <- DT::renderDataTable({
      req(plsda.table)
      # -------------
      DT::datatable(plsda.table, 
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })
    # -------------
    output$plsda_vip_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(plsda_tab,
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })
  })
  
  observeEvent(input$do_rf, {
    
    curr <- as.data.table(mSet$dataSet$preproc)
    curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
    
    config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
    config <- config[!is.na(config$Sample),]
    keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))
    
    config <- cbind(config, Label=mSet$dataSet$cls)
    
    curr <- cbind(config, curr[keep_curr])
    
    predictor = "Label"
    
    inTrain <- caret::createDataPartition(y = config$Label,
                                          ## the outcome data are needed
                                          p = 0.6,#input$rf_train_perc/100,
                                          ## The percentage of data in the training set
                                          list = FALSE)
    trainY <- curr[inTrain, 
                   ..predictor][[1]]
    
    training <- curr[ inTrain,]
    testing <- curr[-inTrain, ]
    
    trainX <- apply(training[, -c(1:ncol(config)), with=FALSE], 2, as.numeric)
    testX <- apply(testing[, -c(1:ncol(config)), with=FALSE], 2, as.numeric)
    
    testY <- testing[,predictor, 
                     with=FALSE][[1]]
    
    model = randomForest::randomForest(trainX, trainY, ntree=500,importance=TRUE)
    result <- randomForest::rfcv(trainX, trainY, cv.fold=10, recursive=TRUE)    
    with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
    
    importance = as.data.table(model$importance, keep.rownames = T)
    rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
    rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
    rf_tab <<- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
    
    prediction <- stats::predict(model,
                                 testX, "prob")[,2]
    
    data <- data.frame(D = as.numeric(as.factor(testY))-1,
                       D.str = testY)
    data <- cbind(data, prediction)
    
    roc_coord <- data.frame(D = rep(data[, "D"], length(3)), M = data[, 3], name = rep(names(data)[3], each = nrow(data)), stringsAsFactors = FALSE)
    #roc_coord <- plotROC::melt_roc(data, "D", m = 3:ncol(data))
    
    ggplot2::ggplot(roc_coord, 
                    ggplot2::aes(d = D, m = M, color = name)) + 
      plotROC::geom_roc(labelsize=0,show.legend = TRUE) + 
      plotROC::style_roc() + 
      theme(axis.text=element_text(size=19),
            axis.title=element_text(size=19,face="bold"),
            legend.title=element_text(size=19),
            legend.text=element_text(size=19))
    
    
  })
  
  observeEvent(input$do_lasnet, {
    
    plot.many <- function(res.obj = models, which_alpha = 1){
      
      predictions <- if(length(res.obj) > 1) do.call("cbind", lapply(res.obj, function(x) x$prediction)) else data.frame(res.obj[[1]]$prediction)
      
      colnames(predictions) <- if(length(res.obj) > 1) sapply(res.obj, function(x) x$alpha) else res.obj[[1]]$alpha
      testY = res.obj[[1]]$labels
      
      if(length(unique(testY)) > 2){
        
        return("not supported yet")
        # https://stats.stackexchange.com/questions/112383/roc-for-more-than-2-outcome-categories
        #
        # for (type.id in length(unique(testY))) {
        #   type = as.factor(iris.train$Species == lvls[type.id])
        #   
        #   nbmodel = NaiveBayes(type ~ ., data=iris.train[, -5])
        #   nbprediction = predict(nbmodel, iris.test[,-5], type='raw')
        #   
        #   score = nbprediction$posterior[, 'TRUE']
        #   actual.class = iris.test$Species == lvls[type.id]
        #   
        #   pred = prediction(score, actual.class)
        #   nbperf = performance(pred, "tpr", "fpr")
        #   
        #   roc.x = unlist(nbperf@x.values)
        #   roc.y = unlist(nbperf@y.values)
        #   lines(roc.y ~ roc.x, col=type.id+1, lwd=2)
        #   
        #   nbauc = performance(pred, "auc")
        #   nbauc = unlist(slot(nbauc, "y.values"))
        #   aucs[type.id] = nbauc
        # }
      }else{
        data <- data.frame(D = as.numeric(as.factor(testY))-1,
                           D.str = testY)
        data <- cbind(data, predictions)
        if(length(res.obj) > 1){
          roc_coord <- plotROC::melt_roc(data, "D", m = 3:ncol(data))
        }else{
          roc_coord <- data.frame(D = rep(data[, "D"], length(3)),
                                  M = data[, 3], 
                                  name = rep(names(data)[3], each = nrow(data)), 
                                  stringsAsFactors = FALSE)
          print(head(roc_coord))
        }
      }
      
      names(roc_coord)[which(names(roc_coord) == "name")] <- "alpha"
      
      roc_coord <- roc_coord[roc_coord$alpha %in% which_alpha,]
      # plot
      ggplot2::ggplot(roc_coord, 
                      ggplot2::aes(d = D, m = M, color = alpha)) + 
        plotROC::geom_roc(labelsize=0,show.legend = TRUE) + 
        plotROC::style_roc() + 
        theme(axis.text=element_text(size=19),
              axis.title=element_text(size=19,face="bold"),
              legend.title=element_text(size=19),
              legend.text=element_text(size=19))
    }
    
    withProgress({
      
      # - - - use filtered data, but not normalized - - -
      
      curr <- as.data.table(mSet$dataSet$preproc)
      curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
      
      config <- mSet$dataSet$batches[match(rownames(mSet$dataSet$preproc),mSet$dataSet$batches$Sample),]
      config <- config[!is.na(config$Sample),]
      keep_curr <- match(mSet$dataSet$batches$Sample,rownames(mSet$dataSet$preproc))
      
      config <- cbind(config, Label=mSet$dataSet$cls)
      
      curr <- cbind(config, curr[keep_curr])
      
      setProgress(0.2)
      
      # - - - - - - - - - - - - - - - - - - - - - - - - - 
      
      input = list(lasnet_alpha = 1, lasnet_trainfrac = 60)
      
      alphas <- input$lasnet_alpha 
      
      models <- glmnet_all_alpha(curr = curr,
                                 nvars = ncol(config) + 1, 
                                 cl = 0,
                                 a = alphas,
                                 perc_train = input$lasnet_trainfrac/100)
      
      #setProgress(0.4)
      
      # - - store results in mset - - -
      
      output$lasnet_roc_plot <- renderPlot({ plot.many(models, 
                                                                             which_alpha = alphas) })
      
      setProgress(0.6)
      
      lasnet_tables <<- lapply(models, function(x){
        table = x$model$beta
        keep = which(table[,1] > 0)
        table = data.frame("beta" = table[keep,1],
                           "absbeta" = abs(table[keep,1]),
                           row.names = rownames(table)[keep])
        colnames(table) <- c("beta", "abs_beta")
        # - - -
        table
      })
      
      lasnet_stab_tables <<- lapply(models, function(x){
        table = data.frame("perc_chosen" = c(x$int_cv_feat))
        rownames(table) <- names(x$int_cv_feat)
        colnames(table) <- "perc_chosen"
        # - - -
        table
      })
      
      names(lasnet_tables) <<- alphas
      names(lasnet_stab_tables) <<- alphas
      
      setProgress(0.8)
      
      # - - plot ROC curves - - -
      
      lapply(1:length(lasnet_tables), function(i){
        output[[paste0("alpha_", names(lasnet_tables)[i], "_lasnet_tab")]] <- DT::renderDataTable({
          DT::datatable(data = as.data.frame(lasnet_tables[i]),
                        selection = 'single',
                        autoHideNavigation = T,
                        options = list(lengthMenu = c(5, 10, 15), 
                                       pageLength = 5))
        })  
        output[[paste0("alpha_stab_", names(lasnet_tables)[i], "_lasnet_tab")]] <- DT::renderDataTable({
          DT::datatable(data = as.data.frame(lasnet_stab_tables[i]),
                        selection = 'single',
                        autoHideNavigation = T,
                        options = list(lengthMenu = c(5, 10, 15), 
                                       pageLength = 5))
        }) 
      })
      
      
      setProgress(0.9)
      
      output$lasnet_table_ui <- renderUI({
        do.call(tabsetPanel, c(id='t',lapply(1:length(lasnet_tables), function(i) {
          tabPanel(
            title=paste0("alpha=",names(lasnet_tables)[i]), 
            tabsetPanel(id="t2", 
                        tabPanel(title="final model",div(DT::dataTableOutput(outputId = paste0("alpha_", names(lasnet_tables)[i], "_lasnet_tab")),style='font-size:80%')),
                        tabPanel(title="feature stats",div(DT::dataTableOutput(outputId = paste0("alpha_stab_", names(lasnet_tables)[i], "_lasnet_tab")),style='font-size:80%'))
            ))
        })))
      })
      
    })
    
  })
  
  observe({
    # ---------------------------------
    db_search_list <- lapply(db_list, 
                             FUN = function(db){
                               # -----------------------
                               if(!input[[paste0("search_", db)]]){
                                 c(file.path(options$db_dir, paste0(db, ".full.db")))
                               }
                               else{NA}
                             }
    )
    db_search_list <<- db_search_list[!is.na(db_search_list)]
  })  
  
  observe({
    # ---------------------------------
    add_search_list <- lapply(db_list, 
                              FUN = function(db){
                                if(is.null(input[[paste0("add_", db)]])){
                                  return(NA)
                                }else{
                                  # -----------------------
                                  if(!input[[paste0("add_", db)]]){
                                    c(file.path(options$db_dir, paste0(db, ".full.db")))
                                  }
                                  else{NA}
                                }
                              }        
    )
    add_search_list <<- add_search_list[!is.na(add_search_list)]
  })  
  
  observe({
    # ---------------------------------
    enrich_db_list <- lapply(db_list, 
                             FUN = function(db){
                               if(is.null(input[[paste0("enrich_", db)]])){
                                 return(NA)
                               }else{
                                 # -----------------------
                                 if(!input[[paste0("enrich_", db)]]){
                                   c(file.path(options$db_dir, paste0(db, ".full.db")))
                                 }
                                 else{NA}
                               }
                             }
    )
    enrich_db_list <<- enrich_db_list[!is.na(enrich_db_list)]
  })  
  
  observe({
    # --------------
    wanted.adducts.pos <- pos_adducts[input$pos_add_tab_rows_selected, "Name"]
    wanted.adducts.neg <- neg_adducts[input$neg_add_tab_rows_selected, "Name"]
    # ---------
    selected_adduct_list <<- rbind(wanted.adducts.neg, 
                                   wanted.adducts.pos)$Name
  })
  
  res.update.tables <<- c("tt", 
                          "fc", 
                          "aov",
                          "rf",
                          "asca", 
                          "meba",
                          "plsda_vip",
                          "enrich_pw",
                          paste0("alpha_", seq(0,1,0.2), "_lasnet"),
                          paste0("alpha_stab_", seq(0,1,0.2), "_lasnet"))
  
  observeEvent(input$enriched_rows_selected, {
    curr_row = input$enriched_rows_selected
    # LOAD MEMBERS IN GROUP
    curr_pw <- enrich_tab[curr_row,]
    pw_memb_here <- gsaRes$gsc[[curr_pw$Name]]
    pw_memb_tab <<- as.data.frame(pw_memb_here)
    rownames(pw_memb_tab) <- pw_memb_here
    pw_memb_tab[,1] <- c("Yes")
    colnames(pw_memb_tab) <- "Present"
    # get missing members
    # pw_memb_all <- gset_proc$gsc[[curr_pw$Name]]
    # pw_memb_all_tab <<- as.data.frame(pw_memb_all)
    # rownames(pw_memb_all_tab) <- pw_memb_all
    # pw_memb_all_tab[,1] <- c("No")
    # colnames(pw_memb_all_tab) <- "Present"
    # join
    enrich_overview_tab <<- pw_memb_tab #rbind(pw_members_tab, pw_memb_all_tab)
    #
    output$enrich_pw_tab <-DT::renderDataTable({
      DT::datatable(enrich_overview_tab,
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), 
                                   pageLength = 5))
    })  
  })
  
  lapply(unique(res.update.tables), FUN=function(table){
    observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
      curr_row = input[[paste0(table, "_tab_rows_selected")]]
      # do nothing if not clicked yet, or the clicked cell is not in the 1st column
      if (is.null(curr_row)) return()
      print("a click has been registered")
      if(grepl(pattern = "lasnet",x = table)){
        print("trigger")
        alpha <- gsub("alpha|stab|_|lasnet", "", table)
        
        if(!grepl(pattern = "stab",x = table)){
          table <- "lasnet_a"
        }else{
          table <- "lasnet_b"
        }
        # - - -
        print(alpha)
        print(table)
      }
      curr_cpd <<- data.table::as.data.table(switch(table,
                                                    tt = mSet$analSet$tt$sig.mat,
                                                    fc = mSet$analSet$fc$sig.mat,
                                                    asca = mSet$analSet$asca$sig.list$Model.ab,
                                                    aov = mSet$analSet$aov$sig.mat,
                                                    rf = vip.score,
                                                    enrich_pw = enrich_overview_tab,
                                                    meba = mSet$analSet$MB$stats,
                                                    plsda_vip = plsda_tab,
                                                    lasnet_a = lasnet_tables[[alpha]],
                                                    lasnet_b = lasnet_stab_tables[[alpha]])
                                             , keep.rownames = T)[curr_row, rn]
      
      if(grepl(pattern = "lasnet",x = table)){
        outplot_name = "lasnet_specific_plot"
      }else{
        outplot_name = paste0(table, "_specific_plot")
      }
      
      output$curr_cpd <- renderText(curr_cpd)
      
      output[[outplot_name]] <- plotly::renderPlotly({
        # --- ggplot ---
        if(table == 'meba'){
          ggplotMeba(curr_cpd, draw.average = T, cols = color.vec(),cf=color.function)
        }else{
          ggplotSummary(curr_cpd, cols = color.vec(), cf=color.function)
        }
      })
      # if(input$autosearch & length(db_search_list > 0)){
      #   match_table <<- multimatch(curr_cpd, db_search_list, mSet$dataSet$grouping)
      #   output$match_tab <-DT::renderDataTable({
      #     DT::datatable(if(mSet$dataSet$grouping == 'pathway') match_table else{match_table[,-c("Description","Structure")]},
      #                   selection = 'single',
      #                   autoHideNavigation = T,
      #                   options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      #   })  
      # }
    })
  })
  
  observeEvent(plotly::event_data("plotly_click"),{
    d <- plotly::event_data("plotly_click")
    req(d)
    # --------------------------------------------------------------
    switch(input$tab_stat,
           pca = {
             if(!"z" %in% names(d)){
               which_group = 1
               which_group <- d$curveNumber + 1
               traceLoc <- length(unique(mSet$dataSet$cls)) + 1
               scatter <- pca_plot$x$attrs[[traceLoc]]
               idx <- scatter$color == which_group
               if(pca_plot$x$data[[which_group]]$visible == "legendonly"){
                 pca_plot$x$data[[which_group]]$visible = TRUE
                 scatter$visible[idx] <- T
               }else{ # hide
                 pca_plot$x$data[[which_group]]$visible = "legendonly"
                 scatter$visible[idx] <- F
               }
               pca_plot$x$attrs[[traceLoc]] <- scatter
               pca_plot <<- pca_plot
               output$plot_pca <- plotly::renderPlotly({pca_plot})
             }
           },
           tt = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$tt$p.value)) return(NULL)
             curr_cpd <<- d$key
             output$tt_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           fc = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$fc$fc.log)) return(NULL)
             curr_cpd <<- d$key
             output$fc_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           aov = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$aov$p.value)) return(NULL)
             curr_cpd <<- d$key
             output$aov_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           rf = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(vip.score)) return(NULL)
             curr_cpd <<- d$key
             output$rf_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           heatmap_biv = {
             if(!exists("hm_matrix")) return(NULL)
             if(d$y > length(hm_matrix$matrix$rows)) return(NULL)
             curr_cpd <<- hm_matrix$matrix$rows[d$y]
             output$heatmap_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           heatmap_mult = {
             if(!exists("hm_matrix")) return(NULL)
             if(d$y > length(hm_matrix$matrix$rows)) return(NULL)
             curr_cpd <<- hm_matrix$matrix$rows[d$y]
             output$heatmap_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           volc = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% rownames(mSet$analSet$volcano$sig.mat)) return(NULL)
             curr_cpd <<- d$key
             output$volc_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           },
           lasnet = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$lasnet$cpds)) return(NULL)
             curr_cpd <<- d$key
             output$lasnet_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)
             })
           })
    # ----------------------------
    output$curr_cpd <- renderText(curr_cpd)
    # if(input$autosearch & length(db_search_list) > 0){
    #   match_table <<- multimatch(curr_cpd, db_search_list, mSet$dataSet$grouping)
    #   output$match_tab <-DT::renderDataTable({
    #     DT::datatable(if(mSet$dataSet$grouping == 'pathway') match_table else{match_table[,-c("Description", "Structure")]},
    #                   selection = 'single',
    #                   autoHideNavigation = T,
    #                   options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    #   })  
    # }
  })
  
  # --- find matches ---
  
  output$find_mol_icon <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('www/search.png'))
    # Return a list containing the filename and alt text
    list(src = filename,
         width=70,
         height=70)
  }, deleteFile = FALSE)
  
  observeEvent(input$search_cpd, {
    req(db_search_list)
    # ----------------
    if(length(db_search_list) > 0){
      match_table <<- multimatch(curr_cpd, db_search_list, mSet$dataSet$grouping)
      output$match_tab <-DT::renderDataTable({
        DT::datatable(match_table[,-c("description","structure")],
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })  
    }
  })
  
  observeEvent(input$match_tab_rows_selected,{
    curr_row = input$match_tab_rows_selected
    curr_row <<- input$match_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_def <<- match_table[curr_row,'description']
    output$curr_definition <- renderText(curr_def$description)
    curr_struct <<- match_table[curr_row,'structure'][[1]]
    output$curr_struct <- renderPlot({plot.mol(curr_struct,style = "cow")})
  })
  
  observeEvent(input$browse_db,{
    req(db_search_list)
    # -------------------
    cpd_list <- lapply(db_search_list, FUN=function(match.table){
      browse_db(match.table)
    })
    # ------------------
    browse_table <<- unique(as.data.table(rbindlist(cpd_list)))
    
    output$browse_tab <-DT::renderDataTable({
      DT::datatable(browse_table[,-c("Description", "Charge")],
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
    })
  })
  
  observeEvent(input$revsearch_cpd, {
    req(db_search_list)
    req(input$browse_tab_rows_selected)
    # -------------------
    search_cmd <- browse_table[curr_row,c('Formula', 'Charge')]
    # -------------------
    cpd_list <- lapply(db_search_list, FUN=function(match.table){
      get_mzs(search_cmd$Formula, search_cmd$Charge, match.table)})
    # ------------------
    hits_table <<- unique(as.data.table(rbindlist(cpd_list)))
    output$hits_tab <-DT::renderDataTable({
      DT::datatable(hits_table,
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
    })
  })
  
  observeEvent(input$browse_tab_rows_selected,{
    curr_row = input$browse_tab_rows_selected
    curr_cpd = browse_table[curr_row, Formula]
    curr_row <<- input$browse_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_def <<- browse_table[curr_row,'Description']
    output$browse_definition <- renderText(curr_def$Description)
    # --- search ---
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = color.vec(),cf=color.function)})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
  })
  
  observeEvent(input$hits_tab_rows_selected,{
    curr_row = input$hits_tab_rows_selected
    curr_row <<- input$hits_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_cpd <<- hits_table[curr_row, mzmed.pgrp]
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = color.vec(),cf=color.function)})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=color.function)})
  })
  
  observeEvent(input$go_enrich,{
    enrich_db_list <- paste0("./backend/db/", c("kegg", "wikipathways", "smpdb", "metacyc"), ".full.db")
    withProgress({
      all_pathways <- lapply(enrich_db_list, FUN=function(db){
        conn <- RSQLite::dbConnect(RSQLite::SQLite(), db) # change this to proper var later
        dbname <- gsub(basename(db), pattern = "\\.full\\.db", replacement = "")
        #RSQLite::dbGetQuery(conn, "SELECT distinct * FROM pathways limit 10;")
        gset <- RSQLite::dbGetQuery(conn,"SELECT DISTINCT c.baseformula AS cpd,
                                    p.name as name
                                    FROM pathways p
                                    JOIN base c
                                    ON c.pathway = p.identifier 
                                    ")
        RSQLite::dbDisconnect(conn)
        # --- returny ---
        gset$name <- paste(gset$name, " (", dbname, ")", sep="")
        gset
      })
      setProgress(0.33)
      # --- only anova top 100 for now ---
      used.analysis <- input$enrich_stats
      if(used.analysis == "rf"){
        ranking <- vip.score[accuracyDrop > 0,]
        ranking <- ranking[order(-accuracyDrop),]
        vec <- ranking$accuracyDrop
        names(vec) <- ranking$rn
        sigvals <- switch(input$enrich_vals,
                          sig=vec[1:500],
                          t50=vec[1:50],
                          t100=vec[1:100],
                          t200=vec[1:200],
                          t500=vec[1:500]
        )
      }else{
        stat.tab <- switch(input$enrich_stats,
                           tt="p.value",
                           aov="p.value",
                           fc="fc.all")
        sigvals <- switch(input$enrich_vals,
                          sig=mSet$analSet[[used.analysis]][[stat.tab]][mSet$analSet[[used.analysis]]$inx.imp],
                          t50=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:50],
                          t100=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:100],
                          t200=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:200],
                          t500=sort(mSet$analSet[[used.analysis]][[stat.tab]])[1:500]
        )        
      }
      
      # -----------------------
      gset <- rbindlist(all_pathways)
      gset_proc <<- piano::loadGSC(gset)
      setProgress(0.66)
      
      # - - - - - - 
      
      sigvals <- rownames(lasnet_tables[[1]])
      matches <- rbindlist(lapply(sigvals, function(mz){
        subtables <- lapply(enrich_db_list, function(db){
          matches = get_matches(mz, db, F, "mz")
        }) 
        tab <- rbindlist(subtables)
        tab$mz = c(mz)
        tab
      }))
      
      lasnet_vals <- lasnet_tables[[1]]
      lasnet_vals$mz <- rownames(lasnet_vals)
      matches <- merge(matches, lasnet_vals, by = "mz")
      sigvals <- unique(matches[, c("Baseformula", "beta")])
      
      # - - - - - -
      
      gsaRes <<- piano::runGSA(sigvals, 
                               gsc = gset_proc)
      enrich_tab <<- piano::GSAsummaryTable(gsaRes)[,1:3]
      # --- render ---
      output$enriched <- DT::renderDataTable({
        # -------------
        DT::datatable(enrich_tab, 
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
        
      })
      setProgress(1)
    })
  })
  
  # --- ON CLOSE ---
  session$onSessionEnded(function() {
    if(any(!is.na(session_cl))) parallel::stopCluster(session_cl)
    R.utils::gcDLLs() # flush dlls
    #save(mSet$dataSet, mSet$analSet, file = file.path(options$work_dir, paste0(options$proj_name, ".RData")))
  })
})


