# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

shinyServer(function(input, output, session) {
  
  # ================================= DEFAULTS ===================================
  
  source('./backend/scripts/joanna/shiny_general.R')

  #theme_update(plot.title = element_text(hjust = 0.5))
  
  shinyOptions(progress.style="old")
  
  parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
    "mape",
    "flattenlist"
  ))
  
  parallel::clusterEvalQ(session_cl, library(data.table))
  
  # - - - - - - - - - - - - -
  
  spinnyimg <- reactiveVal("www/electron.png")
  
  output$spinny <- renderText({spinnyimg()})
  
  output$taskbar_image <- renderImage({
    list(src = file.path(getwd(), 
                         "www", 
                         getOptions("user_options.txt")$taskbar_image), 
         width = 120,
         height = 120,
         style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;")
  }, deleteFile = FALSE)
  
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
                            sardine(fadeImageButton("add_wikidata", img.path = "wikidata.png")),
                            sardine(fadeImageButton("add_vmh", img.path = "vmh.png"))
                            
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
  
  
  # --- render text ---
  
  lapply(global$constants$default.text, FUN=function(default){
    output[[default$name]] = renderText(default$text)
  })
  
  stat.ui.bivar <- reactive({
    navbarPage(inverse=F,h2("Standard analysis"), id="tab_stat",
               # TODO: T-SNE
               tabPanel(h3("PCA"), value = "pca", #icon=icon("cube"),
                        fluidRow(column(12,align="center",plotly::plotlyOutput("plot_pca",height = "600px", width="600px"))),
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
                                                             plotOutput("pca_scree"))
                                        ))
                        )
               ),
               tabPanel(h3("PLSDA"), value = "plsda",
                        fluidRow(
                          selectInput("plsda_type", 
                                      label="PLSDA subtype:", 
                                      choices=list("Normal"="normal",
                                                   "Orthogonal"="ortho",
                                                   "Sparse"="sparse"), 
                                      selected=1),
                          actionButton("do_plsda",label="Go")
                        ),
                        hr(),
                        navbarPage(inverse=F,"",
                                   tabPanel("", icon=icon("globe"),
                                            plotly::plotlyOutput("plot_plsda_3d", 
                                                                 height="800px")
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
                                   ))
               ),
               tabPanel(h3("Heatmap"), value="heatmap_biv",
                        plotly::plotlyOutput("heatmap",width="100%",height="700px"),
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
               ,tabPanel(h3("ML"), value = "ml",
                         fluidRow(
                           column(width=4,
                                  selectInput("ml_method", 
                                              label = "Type:", 
                                              selected = "rf", 
                                              choices = list("Random Forest" = "rf",
                                                             "Lasso" = "ls",
                                                             "Group lasso" = "gls"),
                                              multiple = F)
                           ),
                           column(width=1, br(), tags$a(img(src="help.png"), href="https://regex101.com")),
                           column(width=2, textInput("ml_train_regex", label = "Regex for train:")),
                           column(width=2, textInput("ml_test_regex", label = "Regex for test:")),
                           column(width=2, textInput("ml_name", label="Name:"))
                         ),
                         fluidRow(
                           column(width=5,sliderInput("ml_train_perc", 
                                                      label = "Percentage in training", 
                                                      min = 1,
                                                      max = 100,
                                                      step = 1,
                                                      value = 60, 
                                                      post = "%"),
                                  selectInput("ml_folds", label="Fold CV",choices = c("5", "10", "20", "50", "LOOCV"),multiple = F)),
                           column(width=2, 
                                  #switchButton(inputId = "ml_saturation", label = "SatMode", value = FALSE, col = "BW", type = "OO"),
                                  actionButton("do_ml",label="Go",width = "50px"),style = "margin-top: 35px;", align="left"),
                           column(width=5,sliderInput("ml_attempts", 
                                                      label = "Attempts", 
                                                      min = 1,
                                                      max = 100,
                                                      step = 1,
                                                      value = 20, 
                                                      post = "x")
                           )),
                         
                         hr()
                         ,
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
                                             #,uiOutput("ml_table_ui")
                                    )
                         )
               ))
  })
  
  stat.ui.multivar <- reactive({
    navbarPage(inverse=F,h2("Standard analysis"), id="tab_stat",
               tabPanel("", value = "intro", icon=icon("comment-o"),
                        helpText("Info text here")
               ), # pca_legend
               # TODO: t-sne
               tabPanel(h3("PCA"), value = "pca", #icon=icon("cube"),
                        fluidRow(column(10, plotly::plotlyOutput("plot_pca",height = "600px")),
                                 column(2, br(),br(), br(),plotly::plotlyOutput("pca_legend",height = "400px"))
                        ),
                        fluidRow(column(3,
                                        selectInput("pca_x", label = "X axis:", choices = paste0("PC", 1:30),selected = "PC1",width="100%"),
                                        selectInput("pca_y", label = "Y axis:", choices = paste0("PC", 1:30),selected = "PC2",width="100%"),
                                        selectInput("pca_z", label = "Z axis:", choices = paste0("PC", 1:30),selected = "PC3",width="100%")
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
               # tabPanel(h3("RandomForest"), value="rf",
               #          fluidRow(plotly::plotlyOutput('rf_specific_plot',width="100%")),
               #          navbarPage(inverse=F,"",
               #                     tabPanel("", icon=icon("table"),
               #                              div(DT::dataTableOutput('rf_tab',width="100%"),style='font-size:80%'))
               #          )
               # ),
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
    selectInput('your.time', 'What end time do you want to pick?', choices = get_times(global$paths$patdb))
  }
  
  observeEvent(input$exp_type,{
    req(input$exp_type)
    # ---------------------
    update.UI()
  })
  
  # -----------------
  
  output$pos_add_tab <-DT::renderDataTable({
    # -------------
    DT::datatable(global$vectors$pos_adducts,
                  selection = list(mode = 'multiple', selected = c(1:3, nrow(global$vectors$pos_adducts)), target = 'row'),
                  options = list(pageLength = 5, dom = 'tp'), 
                  rownames = F)
  })
  
  output$neg_add_tab <-DT::renderDataTable({
    # -------------
    DT::datatable(global$vectors$neg_adducts,
                  selection = list(mode = 'multiple', selected = c(1, 2, 14, 15, nrow(global$vectors$neg_adducts)), target = 'row'),
                  options = list(pageLength = 5, dom = 'tp'), 
                  rownames = F)
  })
  
  observeEvent(input$sel_all_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', selected = c(1:nrow(global$vectors$pos_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$neg_adducts,
                    selection = list(mode = 'multiple', selected = c(1:nrow(global$vectors$neg_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  observeEvent(input$sel_no_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', selected = c(0), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$neg_adducts,
                    selection = list(mode = 'multiple', selected = c(0), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  
  observeEvent(input$sel_comm_adducts, {
    output$pos_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$pos_adducts,
                    selection = list(mode = 'multiple', selected = c(1:3, nrow(global$vectors$pos_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
    
    output$neg_add_tab <-DT::renderDataTable({
      # -------------
      DT::datatable(global$vectors$neg_adducts,
                    selection = list(mode = 'multiple', selected = c(1, 2, 14:15, nrow(global$vectors$neg_adducts)), target = 'row'),
                    options = list(pageLength = 5, dom = 'tp'), 
                    rownames = F)
    })
  })
  
  # =================================================
  
  lapply(global$constants$images, FUN=function(image){
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
    # - - - - - -
    pkg_tbl <- get.package.table() #TODO: sort by 'No' first!! (ascending?) - or translate to numeric factor first?
    #pkg_tbl <- rbind(pkg_tbl, data.table("lalalala", 'No', "1.1.1"))
    #pkg_tbl <- pkg_tbl[order(pkg_tbl$Installed,decreasing = T),]
    output$package_tab <- DT::renderDataTable({
    # - - - - - -
      DT::datatable(pkg_tbl,
                    selection = 'none',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(10, 20, 30), pageLength = 10), 
                    rownames = F)
    })
  }
  )
  
  observeEvent(input$update_packages, {
    pacman::p_load(char = global$constants$packages, update = T, character.only = T)
    # - - refresh package list - - 
    pkg_tbl <- get.package.table() 
    output$package_tab <- DT::renderDataTable({
      # - - - - - - - - - - -
      DT::datatable(pkg_tbl,
                    selection = 'none',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(10, 20, 30), pageLength = 10), 
                    rownames = F)
    })
  })

  observe({
    shinyFileChoose(input, 'outlist_pos', roots=global$paths$volumes, filetypes=c('csv'))
    shinyFileChoose(input, 'outlist_neg', roots=global$paths$volumes, filetypes=c('csv'))
    shinyFileChoose(input, 'excel', roots=global$paths$volumes, filetypes=c('xls', 'xlsm', 'xlsx'))
    shinyFileChoose(input, 'database', roots=global$paths$volumes, filetypes=c('sqlite3', 'db', 'sqlite'))
    shinyFileChoose(input, 'taskbar_image_path', roots=global$paths$volumes, filetypes=c('png', 'jpg', 'jpeg', 'bmp'))
  })
  
  observe({
    # - - - - 
    if(!is.list(input$taskbar_image_path)) return()
    print(input$taskbar_image_path)
    img_path <- parseFilePaths(global$paths$volumes, input$taskbar_image_path)$datapath
    new_path <- file.path(getwd(), "www", basename(img_path))
    
    print(img_path)
    print(new_path)
    
    if(img_path != new_path) file.copy(img_path, new_path, overwrite = T)
    # - - -
    output$taskbar_image <- renderImage({
      list(src = new_path, 
           width = 120,
           height = 120,
           style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;") # blocky bg
      #
    }, deleteFile = FALSE)
    # - - -
    setOption('user_options.txt', 'taskbar_image', basename(new_path))
  })
  
  observe({  
    shinyDirChoose(input, "get_db_dir", 
                   roots=global$paths$volumes, 
                   session = session)
     
    test <<- input$get_db_dir
    
    if(typeof(input$get_db_dir) != "list") return()

    print(input$get_db_dir)
    given_dir <- parseDirPath(global$paths$volumes, 
                              input$get_db_dir)
    if(is.null(given_dir)) return()
    # - - connect - -
    setOption("user_options.txt", "db_dir", given_dir)
    output$curr_db_dir <- renderText({getOptions('user_options.txt')$db_dir})
  })
  
  observe({
    shinyDirChoose(input, "get_work_dir",
                   roots = global$paths$volumes,
                   session = session)

    if(typeof(input$get_work_dir) != "list") return()
    
    print(input$get_work_dir)
    given_dir <- parseDirPath(global$paths$volumes,
                              input$get_work_dir)
    if(is.null(given_dir)) return()
    setOption("user_options.txt", "work_dir", given_dir)
    output$curr_exp_dir <- renderText({getOptions('user_options.txt')$work_dir})
  })
  
  observeEvent(input$set_proj_name, {
    proj_name <<- input$proj_name
    if(proj_name == "") return(NULL)
    global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(proj_name,".db", sep=""))
    # --- connect ---
    setOption("user_options.txt", "proj_name", proj_name)
    output$proj_name <<- renderText(proj_name)
  })
  
  observeEvent(input$set_ppm, {
    ppm <<- input$ppm
    output$ppm <<- renderText(ppm)
    # --- connect ---
    setOption("user_options.txt", "ppm", ppm)
  })
  
  observeEvent(input$color_ramp,{
    output$ramp_plot <- plotly::renderPlotly({
      # global$functions$color.function <<- switch(input$color_ramp,
      #                                            "rb"=rainbow,
      #                                            "y2b"=ygobb,
      #                                            "ml1"=matlab.like2,
      #                                            "ml2"=matlab.like,
      #                                            "m2g"=magenta2green,
      #                                            "c2y"=cyan2yellow,
      #                                            "b2y"=blue2yellow,
      #                                            "g2r"=green2red,
      #                                            "b2g"=blue2green,
      #                                            "b2r"=blue2red,
      #                                            "b2p"=cm.colors,
      #                                            "bgy"=topo.colors,
      #                                            "gyw"=terrain.colors,
      #                                            "ryw"=heat.colors,
      #                                            "bw"=blackwhite.colors)

      brew.cols <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", # - - sequential - -
                   "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                   "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", 
                   "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", # - - diverging - -
                   "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3" # - - qualitative - -
                   )
      
      brew.opts <- lapply(brew.cols, function(opt) colorRampPalette(RColorBrewer::brewer.pal(10, opt)))
      names(brew.opts) <- brew.cols
      
      base.opts <- list("rb"=rainbow,
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
      
      all.opts <- append(base.opts, brew.opts)
      
      global$functions$color.function <<- all.opts[[input$color_ramp]]
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        titlefont = list(size = 20)
      )
      
      # --- re-render --- 
      #TODO: SHOULD BE GENERAL PURPOSE FUNCTION w/ listeners
      plotly::plot_ly(z = volcano, 
                      colors = global$functions$color.function(100), 
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
    output$ggplot_theme_example <- renderPlot({
      p <- ggplot(mtcars) + geom_boxplot(aes(x = wt, y = mpg,
                                           colour = factor(gear)))
      p + plot.theme()
    })
  })
  
  observeEvent(input$change_css, {
    # - - - - - - - - - - - - - - -
    setOption("user_options.txt", "col1", input$bar.col.1)
    setOption("user_options.txt", "col2", input$bar.col.2)
    setOption("user_options.txt", "col3", input$bar.col.3)
    setOption("user_options.txt", "col4", input$bar.col.4)
    
    setOption("user_options.txt", "font1", input$font.1)
    setOption("user_options.txt", "font2", input$font.2)
    setOption("user_options.txt", "font3", input$font.3)
    setOption("user_options.txt", "font4", input$font.4)
    
    setOption("user_options.txt", "size1", input$size.1)
    setOption("user_options.txt", "size2", input$size.2)
    setOption("user_options.txt", "size3", input$size.3)
    setOption("user_options.txt", "size4", input$size.4)
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
    facs <- switch(input$exp_type,
           stat = {
             mSet$dataSet$cls
           },
           time_std = {
             if("facA" %not in% names(mSet$dataSet) & "facB" %not in% names(mSet$dataSet)) return(c("Blue", "Pink"))
             lbl.fac <- if(mSet$dataSet$facA.lbl == "Time") "facB" else "facA"
             mSet$dataSet[[lbl.fac]]
           },
           time_fin = {
             mSet$dataSet$cls
           },
           time_min = {           
             mSet$dataSet$cls
           })
    default.colours <- rainbow(length(facs))
    # -------------
    func <- unlist(lapply(seq_along(facs), function(i) {
      input[[paste("col", i, sep="_")]]
    }))
    # - - -
    global$functions$color.vec <<- func
    # - - - 
    func
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
      DF = rhandsontable::hot_to_r(input$adduct_tab)
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
  
  lapply(global$vectors$db_list, FUN=function(db){
    observeEvent(input[[paste0("check_", db)]],{
      db_folder_files <- list.files(getOptions("user_options.txt")$db_dir)
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
  
  lapply(global$vectors$db_list, FUN=function(db){
    observeEvent(input[[paste0("build_", db)]], {
      # ---------------------------
      library(RCurl)
      library(XML)
      library(SPARQL)
      # ---------------------------
      withProgress({
        
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
          "kegg.charge",
          "xmlParse",
          "getURL",
          "mape",
          "flattenlist"
        ))
        
        pkgs = c("data.table", "enviPat", "KEGGREST", "XML", "SPARQL", "RCurl")
        parallel::clusterCall(session_cl, function(pkgs) {
          for (req in pkgs) {
            require(req, character.only = TRUE)
          }
        }, pkgs = pkgs)
        build.base.db(db,
                      outfolder = getOptions("user_options.txt")$db_dir, 
                      cl = session_cl)
        shiny::setProgress(session=session, 0.5)
        build.extended.db(db, 
                          outfolder = getOptions("user_options.txt")$db_dir,
                          adduct.table = adducts, 
                          cl = session_cl, 
                          fetch.limit = 500) #TODO: figure out the optimal fetch limit...
      })
    })
  })
  
  observeEvent(input$score_iso, {

    req(input$iso_score_method)
    
    if(!data.table::is.data.table(global$tables$last_matches)) return(NULL)
    
    if("score" %in% colnames(global$tables$last_matches)){
      global$tables$last_matches <<- global$tables$last_matches[,-"score"]
    }
    
    # - - - 
    
    score_table <- score.isos(global$paths$patdb, method=input$iso_score_method, inshiny=T) 
    
    global$tables$last_matches <<- global$tables$last_matches[score_table, on = c("baseformula", "adduct")]
    
    output$match_tab <-DT::renderDataTable({
      DT::datatable(global$tables$last_matches[,-c("description","structure", "baseformula")],
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })  
  })
  
  # ================== DATA IMPORT ===========================
  
  observeEvent(input$create_db,{
    
    
    # --------------------
    
    global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name, ".db"))
    
    # --------------------
    
    withProgress({
      
      shiny::setProgress(session=session, value= .1,message = "Loading outlists into memory...")
      
      switch(input$new_proj,
             `From DB` = {
               
               req(input$database, input$excel)
               
               db_path <- parseFilePaths(global$paths$volumes, input$database)$datapath
               
               file.copy(db_path, global$paths$patdb, overwrite = T)
               
               shiny::setProgress(session=session, value= .30)
               
               exp_vars <- load.excel(parseFilePaths(global$paths$volumes, input$excel)$datapath, global$paths$patdb)
               
               shiny::setProgress(session=session, value= .60)
               
               #reindex.pat.db(global$paths$patdb)
               
             },
             `From CSV` = {
               
               req(input$outlist_neg, input$outlist_pos, input$excel)
               
               build.pat.db(global$paths$patdb,
                            ppm = ppm,
                            pospath = parseFilePaths(global$paths$volumes, input$outlist_pos)$datapath,
                            negpath = parseFilePaths(global$paths$volumes, input$outlist_neg)$datapath,
                            overwrite = T)
               
               shiny::setProgress(session=session, value= .95,message = "Adding excel sheets to database...")
               
               exp_vars <<- load.excel(parseFilePaths(global$paths$volumes, input$excel)$datapath, global$paths$patdb)}
             )
             
      })
  })
  
  observeEvent(input$import_db, {
    req(input$pat_db)
    global$paths$patdb <<- input$pat_db$datapath
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
    global$paths$csv_loc <<- input$pat_csv$datapath
    
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
    # ---------
    req(getOptions("user_options.txt")$proj_name)
    req(getOptions("user_options.txt")$work_dir)
    # ---------
    withProgress({
      shiny::setProgress(session=session, value= 1/4)
      tbl <- get.csv(global$paths$patdb,
                     time.series = if(input$exp_type == "time_std") T else F,
                     group_adducts = if(length(global$vectors$add_search_list) == 0) F else T,
                     groupfac = input$group_by,
                     which_dbs = global$vectors$add_search_list,
                     which_adducts = selected_adduct_list
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
      shiny::setProgress(session=session, value= 2/4)
      global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
      fwrite(tbl.adj, global$paths$csv_loc, sep="\t")
      # --- overview table ---
      
      as.numi <- as.numeric(colnames(tbl.adj)[1:100])
      
      exp.vars <- which(is.na(as.numi))
      
      #nvars <- length(tbl.adj[,1:(which(colnames(tbl.adj) == "$"))])
      
      shiny::setProgress(session=session, value= 3/4)
      
      output$csv_tab <-DT::renderDataTable({
        overview_tab <- if(input$exp_type == "time_std"){
          t(data.table(keep.rownames = F,
                       Identifiers = ncol(tbl.adj) - length(exp.vars),
                       Samples = nrow(tbl.adj),
                       Times = length(unique(tbl.adj$Time))
          ))
        } else{
          t(data.table(keep.rownames = F,
                       Identifiers = ncol(tbl.adj) - length(exp.vars),
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
  
  # ===================== METABOSTART ========================
  
  observeEvent(input$check_csv, {
    # ----------------------
    
    csv <- fread(global$paths$csv_loc, 
                 header = T)
    
    as.numi <- as.numeric(colnames(csv)[1:100])
    
    exp.vars <- which(is.na(as.numi))
    
    opts <<- colnames(csv[,..exp.vars])
    
    batch <<- which(sapply(exp.vars, function(x) length(unique(csv[,..x][[1]])) < nrow(csv)))
    
    numi <<- which(sapply(exp.vars, function(x) is.numeric(csv[,..x][[1]])))
    
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
    req(global$paths$csv_loc)
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
      shiny::setProgress(session=session, value= .1)
      # ------------------------------
      #Below is your R command history: 
      mSet <<- switch(input$exp_type,
                      stat = {
                      #  mSet<- #DEBUG
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
      
      # ------- load and re-save csv --------

      csv_orig <- fread(global$paths$csv_loc, 
                        data.table = TRUE,
                        header = T)
      
      # replace all 0's with NA
      
      csv_orig[,(1:ncol(csv_orig)) := lapply(.SD,function(x){ ifelse(x == 0, NA, x)})]
      
      # - - -
      
      csv_orig$Sample <- gsub(csv_orig$Sample, pattern=" ", replacement="")
      
      as.numi <- as.numeric(colnames(csv_orig)[1:100])
      
      exp.vars <- which(is.na(as.numi))
      
      # --- batches ---
      
      # input <- list(batch_var = "Batch",
      #               exp_var = "Group",
      #               perc_limit = .99,
      #               filt_type = "none",
      #               miss_type = "rf",
      #               norm_type = "SumNorm",
      #               trans_type = "LogNorm",
      #               scale_type = "AutoNorm",
      #               ref_var = "none",
      #               remove_outliers = TRUE
      #               )
    
      batches <- input$batch_var
      condition <- input$exp_var
      
      if(is.null(batches)) batches <- ""
      
      batch_corr <- if(length(batches) == 1 & batches[1] == "") FALSE else TRUE
      
      if("Batch" %in% batches){ batches = c(batches, "Injection") }
      
      first_part <<- csv_orig[,..exp.vars, with=FALSE]
      first_part[first_part == "" | is.null(first_part)] <- "Unknown"
      
      csv_temp <- cbind(first_part[,!duplicated(names(first_part)),with=FALSE],
                        "Label" = first_part[,..condition][[1]],
                        csv_orig[,-..exp.vars,with=FALSE])
      
      # --- remove the rest ---
      
      csv_subset <- csv_temp
      
      # --- remove outliers? ---
      if(input$remove_outliers){
        
        
        sums <- rowSums(csv_subset[, -exp.vars,with=FALSE],na.rm = TRUE)
        names(sums) <- csv_subset$Sample
        outliers = c(car::Boxplot(as.data.frame(sums)))
        
        csv_temp_no_out <- csv_subset[!(Sample %in% outliers),]
      } else{
        csv_temp_no_out <- csv_subset
      }
      
      # - - - remove peaks that are missing in all - - -
      
      csv_temp_no_out <- csv_temp_no_out[,which(unlist(lapply(csv_temp_no_out, function(x)!all(is.na(x))))),with=F]
      
      # # - - - low signal samples - - -
      
      complete.perc <- rowMeans(!is.na(csv_temp_no_out))
      keep_samps <- csv_temp_no_out$Sample[which(complete.perc > .2)]
      
      csv_temp_no_out <- csv_temp_no_out[Sample %in% keep_samps,]
      
      covar_table <- first_part[Sample %in% keep_samps,]
      
      batchview = if(condition == "Batch") TRUE else FALSE
      
      if(any(grepl("QC", csv_temp_no_out$Sample))){
        samps <- which(!grepl(csv_temp_no_out$Sample, pattern = "QC"))
        batchnum <- unique(csv_temp_no_out[samps, "Batch"][[1]])
        keep_samps_post_qc <- covar_table[which(covar_table$Batch %in% batchnum),"Sample"][[1]]
        covar_table <- covar_table[which(covar_table$Batch %in% batchnum),]
        csv_temp_no_out <- csv_temp_no_out[which(csv_temp_no_out$Sample %in% keep_samps_post_qc),-"Batch"]
        #csv_temp_filt <- csv_temp_no_out#[, -"Batch"]
      }
      
      # remove all except sample and time in saved csv
      exp_var_names <- colnames(csv_temp_no_out)[exp.vars]
      keep_cols <- c("Sample", "Label")# "Group", "Time")
      remove <- which(!(exp_var_names %in% keep_cols))
      
      csv_loc_no_out <- gsub(pattern = "\\.csv", replacement = "_no_out.csv", x = global$paths$csv_loc)
      
      if(file.exists(csv_loc_no_out)) file.remove(csv_loc_no_out)
      
      fwrite(csv_temp_no_out[,-remove,with=F], file = csv_loc_no_out)
      
      rownames(covar_table) <- covar_table$Sample
      
      # -------------------------------------
      
      mSet <<- Read.TextData(mSet, 
                            filePath = csv_loc_no_out, 
                            "rowu")
      
      file.remove(csv_loc_no_out)
      
      mSet$dataSet$covars <<- covar_table
      
      # - - - sanity check - - -
      
      mSet <<- SanityCheckData(mSet)
      
      mSet <<- RemoveMissingPercent(mSet, 
                                    percent = input$perc_limit/100)
      
      if(input$miss_type != "none"){
        if(input$miss_type == "pmm"){
          require(mice)
          #base <- mSet$dataSet$orig
          base <- mSet$dataSet$preproc
          imp <- mice::mice(base, printFlag = TRUE)
          #mSet$dataSet$norm <<- imp
        }else if(input$miss_type == "rf"){
          samples <- rownames(mSet$dataSet$preproc)
          
          w.missing <- mSet$dataSet$preproc#[,1:50]
          w.missing <- apply(w.missing, 2, as.numeric)

          library(doParallel)

          registerDoParallel(session_cl)
          
          auto.mtry <- floor(sqrt(ncol(mSet$dataSet$preproc)))
          
          mtry <- ifelse(auto.mtry > 100, 100, auto.mtry)
          
          imp <- missForest::missForest(w.missing, 
                                        parallelize = "variables",
                                        verbose = T,
                                        ntree = 10,
                                        mtry = mtry)
          
          mSet$dataSet$procr <<- imp$ximp
          rownames(mSet$dataSet$procr) <<- rownames(mSet$dataSet$preproc)
          # - - - - - - - - - - - - 
        }else{
          mSet <<- ImputeVar(mSet,
                             method = # "knn"
                              input$miss_type
                              )
        }
      }
      # ------------------------
      
      if(input$filt_type != "none"){
        
        mSet <<- FilterVariable(mSet,
                                filter = input$filt_type,
                                qcFilter = "F",
                                rsd = 25)
      }
      # ------------------------------------
      
      shiny::setProgress(session=session, value= .2)
      
      if(input$norm_type == "SpecNorm"){
        norm.vec <<- mSet$dataSet$covars[match(mSet$dataSet$covars$Sample,
                                               rownames(mSet$dataSet$preproc)),][[input$samp_var]]
        norm.vec <<- scale(x = norm.vec, center = 1)
        
      }else{
        norm.vec <<- rep(1, length(mSet$dataSet$cls))
      }
      
      mSet <<- Normalization(mSet,
                             rowNorm = input$norm_type,
                             transNorm = input$trans_type,
                             scaleNorm = input$scale_type,
                             ref = input$ref_var)
      
      shiny::setProgress(session=session, value= .4)
      
      smps <- rownames(mSet$dataSet$norm)
      
      qc_rows <- which(grepl(pattern = "QC", x = smps))
      
      # - - - - - - - - - - - - - - - - - - - - - - -
      
      has.qc <- length(qc_rows) > 0
        
      if(batch_corr){
        
        if("Batch" %in% input$batch_var & has.qc){
          
          smpnames = smps[qc_rows]
          
          batch.idx = as.numeric(as.factor(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$Sample),"Batch"][[1]]))
          seq.idx = as.numeric(mSet$dataSet$covars[match(smps, mSet$dataSet$covars$Sample),"Injection"][[1]])
          
          # --- QC CORRECTION ---
          
          corr_cols <- pbapply::pblapply(1:ncol(mSet$dataSet$norm), function(i){
            #
            vec = mSet$dataSet$norm[,i]
            #
            corr_vec = BatchCorrMetabolomics::doBC(Xvec = as.numeric(vec), 
                                                   ref.idx = as.numeric(qc_rows), 
                                                   batch.idx = batch.idx,
                                                   seq.idx = seq.idx)
            # ---------
            corr_vec
          })
          
          qc_corr_matrix <- as.data.frame(do.call(cbind, corr_cols))
          
          colnames(qc_corr_matrix) <- colnames(mSet$dataSet$norm)
          rownames(qc_corr_matrix) <- rownames(mSet$dataSet$norm)
          
          mSet$dataSet$norm <<- as.data.frame(qc_corr_matrix)
          
        }
        
        if(!batchview){
          mSet$dataSet$norm <<- mSet$dataSet$norm[!grepl(rownames(mSet$dataSet$norm),pattern= "QC"),]
          mSet$dataSet$cls <<- mSet$dataSet$cls[which(!grepl(rownames(mSet$dataSet$norm),pattern= "QC")), drop = TRUE]
          #mSet$dataSet$covars <<- mSet$dataSet$covars[!grepl(mSet$dataSet$covars$Sample,pattern= "QC"),]
          mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
        }
        
        left_batch_vars <- grep(input$batch_var, 
                                pattern =  ifelse(has.qc, "Batch|Injection|Sample", "Injection|Sample"),
                                value = T,
                                invert = T)
        
        if(length(left_batch_vars) > 2){ 
          NULL
        } else if(length(left_batch_vars) == 0){
          NULL
        } else{
          
          smp <- rownames(mSet$dataSet$norm)
          exp_lbl <- mSet$dataSet$cls
          
          csv <- as.data.table(cbind(Sample = smp, 
                                     Label = mSet$dataSet$cls,
                                     mSet$dataSet$norm))
          
          csv_edata <-t(csv[,!c(1,2)])
          colnames(csv_edata) <- csv$Sample
          
          if(length(left_batch_vars) == 1){
            csv_pheno <- data.frame(sample = 1:nrow(csv),
                                    batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$Sample),left_batch_vars[1], with=FALSE][[1]],
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
                                    batch1 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$Sample), left_batch_vars[1], with=FALSE][[1]],
                                    batch2 = mSet$dataSet$covars[match(smp, mSet$dataSet$covars$Sample), left_batch_vars[2], with=FALSE][[1]],
                                    outcome = as.factor(exp_lbl))
            batch_normalized = t(limma::removeBatchEffect(x = csv_edata,
                                                          #design = mod.pheno,
                                                          batch = csv_pheno$batch1,
                                                          batch2 = csv_pheno$batch2))
            rownames(batch_normalized) <- rownames(mSet$dataSet$norm)
          }
          
          mSet$dataSet$norm <<- as.data.frame(batch_normalized)
        }
      } else{
        if(!batchview){
          mSet$dataSet$norm <<- mSet$dataSet$norm[!grepl(rownames(mSet$dataSet$norm),pattern= "QC"),]
          mSet$dataSet$cls <<- mSet$dataSet$cls[which(!grepl(rownames(mSet$dataSet$norm),pattern= "QC")), drop = TRUE]
          mSet$dataSet$covars <<- mSet$dataSet$covars[!grepl(mSet$dataSet$covars$Sample,pattern= "QC"),]
          mSet$dataSet$cls.num <<- length(levels(mSet$dataSet$cls))
        }
      }
      
      # REORDER COVARS

      mSet$dataSet$covars <<- mSet$dataSet$covars[match(rownames(mSet$dataSet$norm), 
                                                        mSet$dataSet$covars$Sample), 
                                                  -"#"]
      
      shiny::setProgress(session=session, value= .5)
      
      shiny::setProgress(session=session, value= .6)
      varNormPlots <- ggplotNormSummary(mSet)
      output$var1 <- renderPlot(varNormPlots$tl)
      output$var2 <- renderPlot(varNormPlots$bl)
      output$var3 <- renderPlot(varNormPlots$tr)
      output$var4 <- renderPlot(varNormPlots$br)
      
      shiny::setProgress(session=session, value= .7)
      sampNormPlots <-  ggplotSampleNormSummary(mSet)
      output$samp1 <- renderPlot(sampNormPlots$tl)
      output$samp2 <- renderPlot(sampNormPlots$bl)
      output$samp3 <- renderPlot(sampNormPlots$tr)
      output$samp4 <- renderPlot(sampNormPlots$br)
      shiny::setProgress(session=session, value= .8)
      
      mSet$dataSet$adducts <<- selected_adduct_list
      update.UI()  
      shiny::setProgress(session=session, value= .9)
      
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
                 iPCA.Anal(mSet, file.path(getOptions("user_options.txt")$work_dir, "ipca_3d_0_.json"))
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
               performMB(10, dir=getOptions("user_options.txt")$work_dir)
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
               CalculateImpVarCutoff(0.05, 0.9, dir=getOptions("user_options.txt")$work_dir)
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
  # csv <- as.data.table(fread(global$paths$csv_loc))
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
             # CAN use sparse PCA? linear combinations of just a few orig dimensions.
             output$pca_scree <- renderPlot({
               df <- data.table(
                 pc = 1:length(names(mSet$analSet$pca$variance)),
                 var = mSet$analSet$pca$variance)
               p <- ggplot2::ggplot(data=df[1:20,]) + ggplot2::geom_line(mapping = aes(x=pc, y=var, colour=var), cex=3) + 
                 plot.theme(base_size = 10) +
                 ggplot2::scale_colour_gradientn(colours = global$functions$color.function(256)) +
                 theme(axis.text=element_text(size=10),
                       axis.title=element_text(size=19,face="bold"),
                       #legend.title=element_text(size=15),
                       #legend.text=element_text(size=12),
                       legend.position="none")
               # - - - - - 
               p
               #ggplotly(p)
             })
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
                   title = gsubfn::fn$paste("$x ($x.var %)")),
                 yaxis = list(
                   titlefont = list(size = 20),
                   title = gsubfn::fn$paste("$y ($y.var %)")),
                 zaxis = list(
                   titlefont = list(size = 20),
                   title = gsubfn::fn$paste("$z ($z.var %)")))) 
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
                 translator <- data.table(Sample=rownames(mSet$dataSet$norm),Group=mSet$dataSet$cls)
                 group_assignments <- translator[,"Group",on=colnames(final_matrix)]$Group
                 hm_matrix <<- heatmaply::heatmapr(final_matrix, 
                                                   Colv=T, 
                                                   Rowv=T,
                                                   col_side_colors = group_assignments,
                                                   #dist_method = ..., # use dbscan?
                                                   #hclust_method = ..., ???
                                                   k_row = NA)
                 heatmaply::heatmaply(hm_matrix,
                                      Colv=F, Rowv=F,
                                      branches_lwd = 0.3,
                                      margins = c(60,0,NA,50),
                                      colors = global$functions$color.function(256),
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
             
             sourceTable <- mSet$analSet[[used.analysis]]$inx.imp[mSet$analSet[[used.analysis]]$inx.imp == TRUE]
             if(length(sourceTable) == 0){
               ordered <- order(mSet$analSet[[used.analysis]][[used.values]], decreasing = F)
               sourceTable <- mSet$analSet[[used.analysis]][[used.values]][ordered][1:100]
               #sourceTable <- mSet$analSet[[used.analysis]]$inx.imp[mSet$analSet[[used.analysis]]$inx.imp == TRUE]
             }
             
             #sourceTable <- sort(mSet$analSet[[used.analysis]][[used.values]])[1:100]
             # --- render ---
             withProgress({
               output$heatmap <- plotly::renderPlotly({
                 x <- mSet$dataSet$norm[,names(sourceTable)]
                 final_matrix <- t(x)
                 translator <- data.table(Sample=rownames(mSet$dataSet$norm),Group=mSet$dataSet$cls)
                 group_assignments <- translator[,"Group",on=colnames(final_matrix)]$Group
                 
                 color.mapper <- {
                  classes <- levels(mSet$dataSet$cls)
                  cols <- sapply(1:length(classes), function(i) color.vec()[i])
                  names(cols) <- classes
                  # - - -
                  cols
                 }
                 
                 assignment.df <- as.data.frame(group_assignments)
                 colnames(assignment.df) <- "Group"
                 
                 heatmaply::heatmaply(final_matrix,
                                      Colv=T, 
                                      Rowv=T,
                                      branches_lwd = 0.3,
                                      margins = c(60, 0, NA, 50),
                                      colors = global$functions$color.function(256),
                                      col_side_colors = assignment.df,
                                      col_side_palette = color.mapper,
                                      subplot_widths = c(.9,.1),
                                      subplot_heights =  c(.1,.05,.85),
                                      column_text_angle = 90,
                                      xlab = "Sample",
                                      ylab = "m/z",
                                      showticklabels = c(T,F)
                                      #label_names = c("m/z", "sample", "intensity") #breaks side colours
                                      )
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
               ggPlotTT(global$functions$color.function, 20)
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
               ggPlotFC(global$functions$color.function, 20)
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
           # rf={
           #   if(!"rf" %in% names(mSet$analSet)){
           #     withProgress({
           #       mSet <<- RF.Anal(mSet, 500,7,1)
           #     })
           #   }
           #   vip.score <<- as.data.table(mSet$analSet$rf$importance[, "MeanDecreaseAccuracy"],keep.rownames = T)
           #   colnames(vip.score) <<- c("rn", "accuracyDrop")
           #   output$rf_tab <-DT::renderDataTable({
           #     # -------------
           #     DT::datatable(vip.score, 
           #                   selection = 'single',
           #                   autoHideNavigation = T,
           #                   options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           #   })
           # },
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
               ggPlotVolc(global$functions$color.function, 20)
             })
           }
    )
  })
  
  observeEvent(input$do_plsda, {
    
    library(e1071)
    
    switch(input$plsda_type,
           normal={
             withProgress({
               
               require(caret)

               mSet <<- PLSR.Anal(mSet,
                                  TRUE)
               shiny::setProgress(session=session, value= 0.3)
               
               mSet <<- PLSDA.CV(mSet,compNum = 5)
               mSet <<- PLSDA.Permut(mSet,num = 300, type = "accu")
               #swine_mset <<- PLSDA.Permut(swine_mset,num = 300, type = "accu")
               
               output$plsda_cv_plot <- renderPlot({ggPlotClass(cf = global$functions$color.function)})
               output$plsda_perm_plot <- renderPlot({ggPlotPerm(cf = global$functions$color.function)})
               
               
               # - - - - - - - - - - - 
               shiny::setProgress(session=session, value= 0.6)
               # - - overview table - -
               plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                                  / mSet$analSet$plsr$Xtotvar 
                                                  * 100.0,
                                                  digits = 2),
                                            keep.rownames = T)
               colnames(plsda.table) <- c("Principal Component", "% variance")
               plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
               shiny::setProgress(session=session, value= 0.9)
               # --- coordinates ---
               coords <- mSet$analSet$plsr$scores
               colnames(coords) <- paste0("PC", 1:ncol(coords))
               # --- vip table ---
               colnames(mSet$analSet$plsda$vip.mat) <<- paste0("PC", 1:ncol(mSet$analSet$plsda$vip.mat))
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
          adj_plot$x$data[[i]]$hoverinfo <- "none"
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
          opacity=1,
          hoverinfo = 'text',
          text = rownames(mSet$dataSet$norm)
        ) %>%  layout(scene = list(
          aspectmode="cube",
          xaxis = list(
            titlefont = list(size = 20),
            title = gsubfn::fn$paste("$x ($x.var %)")),
          yaxis = list(
            titlefont = list(size = 20),
            title = gsubfn::fn$paste("$y ($y.var %)")),
          zaxis = list(
            titlefont = list(size = 20),
            title = gsubfn::fn$paste("$z ($z.var %)")
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
  
  
  observeEvent(input$do_ml, {
    
    # do_ml, ml_attempts, ml_train_perc
    
    
    # -- NOTE : CHECK IF AUC INCREASES W/ INCREASING TRAINING SET SIZE --- !!!!!!!!!!!!!!!!
    # 'saturation point'
    #switchButton(inputId = "ml_saturation", label = "SatMode", value = FALSE, col = "BW", type = "OO"),
    
    withProgress({
      
      setProgress(value = 0)
      # prepare matric
      
      curr <- as.data.table(mSet$dataSet$preproc)
      
      curr[,(1:ncol(curr)) := lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
      
      config <- mSet$dataSet$covars[match(mSet$dataSet$covars$Sample,rownames(mSet$dataSet$preproc)),]
      config <- config[!is.na(config$Sample),]
      config <- cbind(config, Label=mSet$dataSet$cls)
      
      # bye bye NAs
      
      config <- config[,!apply(is.na(config), 2, any), with=FALSE]
      
      
      # - - - - - -
      keep_curr <- match(mSet$dataSet$covars$Sample,rownames(mSet$dataSet$preproc))
      
      curr <- cbind(config, curr[keep_curr])
      
      curr <- curr[which(!grepl(curr$Sample,
                                pattern = "QC"))]
      
      # remove cols with all NA
      
      curr <- curr[,colSums(is.na(curr))<nrow(curr),with=FALSE]
      
      configCols <- which(!(colnames(curr) %in% colnames(mSet$dataSet$norm)))
      mzCols <- which(colnames(curr) %in% colnames(mSet$dataSet$norm))
      
      curr[,(configCols):= lapply(.SD, function(x) as.factor(x)), .SDcols = configCols]
      curr[,(mzCols):= lapply(.SD, function(x) as.numeric(x)), .SDcols = mzCols]
      
      goes = as.numeric(input$ml_attempts)
      
      repeats <- pbapply::pblapply(1:goes, function(i){
        
        # train / test
        
        ml_train_regex <<- input$ml_train_regex
        ml_test_regex <<- input$ml_test_regex
        
        ml_train_perc <- input$ml_train_perc/100
        
        if(ml_train_regex == "" & ml_test_regex == ""){ # BOTH ARE NOT DEFINED
          test_idx = caret::createDataPartition(y = curr$Label, p = ml_train_perc, list = FALSE)
          train_idx = setdiff(1:nrow(curr), test_idx)
          inTrain = train_idx
          inTest = test_idx
        }else if(ml_train_regex != ""){ #ONLY TRAIN IS DEFINED
          train_idx = grep(config$Sample, pattern = ml_train_regex)
          test_idx = setdiff(1:nrow(curr), train_idx)
          reTrain <- caret::createDataPartition(y = config[train_idx, Label], p = ml_train_perc)
          inTrain <- train_idx[reTrain$Resample1]
          inTest = test_idx
        }else{ # ONLY TEST IS DEFINED
          test_idx = grep(config$Sample, pattern = ml_test_regex)
          train_idx = setdiff(1:nrow(curr), test_idx)
          reTrain <- caret::createDataPartition(y = config[train_idx, Label], p = ml_train_perc)
          inTrain <- train_idx[reTrain$Resample1]
          
          inTest <- test_idx
          #reTest <- caret::createDataPartition(y = config[test_idx, Label], p = ml_train_perc)
          #inTest <- test_idx[reTest$Resample1]
        }
        
        # - - - re-split - - -
        #reTrain <- caret::createDataPartition(y = config[train_idx, Label], p = ml_train_perc)
        #inTrain <- train_idx[reTrain$Resample1]
        #reTest <- caret::createDataPartition(y = config[test_idx, Label], p = ml_train_perc)
        #inTest <- test_idx[reTest$Resample1]
        
        # - - divide - -
        
        predictor = "Label"
        
        trainY <- curr[inTrain, 
                       ..predictor][[1]]
        testY <- curr[inTest,
                      ..predictor][[1]]
        
        group.cols <- grep(colnames(curr), pattern = "^Group",value = T)
        
        remove.cols <- c("Sample", "Label", group.cols, "Stool_condition")
        remove.idx <- which(colnames(curr) %in% remove.cols)
        training <- curr[inTrain,-remove.idx, with=FALSE]
        testing <- curr[inTest,-remove.idx, with=FALSE]
        
        predIdx <- which(colnames(curr) %in% colnames(config))
        
        #trainX <- apply(training, 2, as.numeric)
        #testX <- apply(testing, 2, as.numeric)
        
        training <- data.matrix(gdata::drop.levels(training))
        testing <- data.matrix(gdata::drop.levels(testing))
        
        setProgress(value = i/goes)
        
        switch(input$ml_method,
               rf = {
                 model = randomForest::randomForest(x = training, 
                                                    y = trainY, 
                                                    ntree = 500,
                                                    importance=TRUE)
                 
                 prediction <- stats::predict(model, 
                                              testing, 
                                              "prob")[,2]
                 
                 importance = as.data.table(model$importance, keep.rownames = T)
                 rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
                 rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
                 rf_tab <- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
                 list(type="rf",
                      feats = as.data.table(rf_tab, keep.rownames = T), 
                      model = model,
                      prediction = prediction,
                      labels = testY)
               }, 
               ls = {

                 nfold = switch(input$ml_folds, 
                                "5" = 5,
                                "10" = 10,
                                "20" = 20,
                                "50" = 50,
                                "LOOCV" = length(trainY))
                 family = "binomial"
                 
                 cv1 <- glmnet::cv.glmnet(training, trainY, family = family, type.measure = "auc", alpha = 1, keep = TRUE, nfolds=nfold)
                 cv2 <- data.frame(cvm = cv1$cvm[cv1$lambda == cv1[["lambda.min"]]], lambda = cv1[["lambda.min"]], alpha = 1)
                 
                 model <- glmnet::glmnet(as.matrix(training), trainY, family = family, lambda = cv2$lambda, alpha = cv2$alpha)
                 
                 prediction <- stats::predict(model,
                                              type = "response", 
                                              newx = testing, 
                                              s = "lambda.min")#[,2] # add if necessary
                 list(type = "ls",
                      model = model,
                      prediction = prediction, 
                      labels = testY)
               }, 
               gls = {
                 NULL
               })
      })
      
      if(!"ml" %in% names(mSet$analSet)){
        
        mSet$analSet$ml <<- list(ls=list(), rf=list())
      }
      
      #input <- list(ml_attempts = 50, ml_method = "ls", ml_name = "test")
      
      xvals <<- list(type = {unique(lapply(repeats, function(x) x$type))},
                     models = {lapply(repeats, function(x) x$model)},
                     predictions = {lapply(repeats, function(x) x$prediction)},
                     labels = {lapply(repeats, function(x) x$labels)})
      
      mSet$analSet$ml[[input$ml_method]][[input$ml_name]] <<- list("roc" = xvals,
                                                                  "bar" = repeats)
      
      output$ml_roc <- plotly::renderPlotly({plotly::ggplotly(ggPlotROC(xvals, input$ml_attempts, global$functions$color.function))})
      output$ml_bar <- plotly::renderPlotly({plotly::ggplotly(ggPlotBar(repeats, input$ml_attempts, global$functions$color.function, input$ml_top_x, input$ml_name))})
      
    })
  })
  
  observe({
    # ---------------------------------
    db_search_list <- lapply(global$vectors$db_list, 
                             FUN = function(db){
                               # -----------------------
                               if(!input[[paste0("search_", db)]]){
                                 c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db")))
                               }
                               else{NA}
                             }
    )
    global$vectors$db_search_list <<- db_search_list[!is.na(db_search_list)]
  })  
  
  observe({
    # ---------------------------------
    add_search_list <- lapply(global$vectors$db_list, 
                              FUN = function(db){
                                if(is.null(input[[paste0("add_", db)]])){
                                  return(NA)
                                }else{
                                  # -----------------------
                                  if(!input[[paste0("add_", db)]]){
                                    c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db")))
                                  }
                                  else{NA}
                                }
                              }        
    )
    global$vectors$add_search_list <<- add_search_list[!is.na(add_search_list)]
  })  
  
  observe({
    # ---------------------------------
    enrich_db_list <- lapply(global$vectors$db_list, 
                             FUN = function(db){
                               if(is.null(input[[paste0("enrich_", db)]])){
                                 return(NA)
                               }else{
                                 # -----------------------
                                 if(!input[[paste0("enrich_", db)]]){
                                   c(file.path(getOptions("user_options.txt")$db_dir, paste0(db, ".full.db")))
                                 }
                                 else{NA}
                               }
                             }
    )
    global$vectors$enrich_db_list <<- enrich_db_list[!is.na(enrich_db_list)]
  })  
  
  observe({
    # --------------
    wanted.adducts.pos <- global$vectors$pos_adducts[input$pos_add_tab_rows_selected, "Name"]
    wanted.adducts.neg <- global$vectors$neg_adducts[input$neg_add_tab_rows_selected, "Name"]
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
                          "ml",
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
      
      if(grepl(pattern = "lasnet",x = table)){
        
        alpha <- gsub("alpha|stab|_|lasnet", "", table)
        if(!grepl(pattern = "stab",x = table)){
          table <- "lasnet_a"
        }else{
          table <- "lasnet_b"
        }
        # - - -
        
        
      }
      curr_cpd <<- data.table::as.data.table(switch(table,
                                                    tt = mSet$analSet$tt$sig.mat,
                                                    fc = mSet$analSet$fc$sig.mat,
                                                    ml = ml_tab,
                                                    asca = mSet$analSet$asca$sig.list$Model.ab,
                                                    aov = mSet$analSet$aov$sig.mat,
                                                    rf = vip.score,
                                                    enrich_pw = enrich_overview_tab,
                                                    meba = mSet$analSet$MB$stats,
                                                    plsda_vip = plsda_tab)
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
          ggplotMeba(curr_cpd, draw.average = T, cols = color.vec(),cf=global$functions$color.function)
        }else{
          ggplotSummary(curr_cpd, cols = color.vec(), cf=global$functions$color.function)
        }
      })
    })
  })
  
  observeEvent(plotly::event_data("plotly_click"),{
    d <- plotly::event_data("plotly_click")
    req(d)
    # --------------------------------------------------------------
    switch(input$tab_stat,
           ml = {
             
             switch(input$ml_results, roc = {
               attempt = d$curveNumber - 1
               if(attempt > 1){
                 ml_type <- xvals$type[[1]]
                 model <- xvals$models[[attempt]]
                 output$ml_tab <- switch(ml_type,
                                         rf = {
                                           importance = as.data.table(model$importance, keep.rownames = T)
                                           rf_tab <- importance[which(MeanDecreaseAccuracy > 0), c("rn", "MeanDecreaseAccuracy")]
                                           rf_tab <- rf_tab[order(MeanDecreaseAccuracy, decreasing = T)]
                                           # - - - return - - -
                                           ml_tab <<- data.frame(MDA = rf_tab$MeanDecreaseAccuracy, row.names = rf_tab$rn) 
                                           DT::renderDataTable({
                                             DT::datatable(rf_tab,
                                                           selection = 'single',
                                                           autoHideNavigation = T,
                                                           options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                                           })
                                         },
                                         ls = {
                                           tab = model$beta
                                           keep = which(tab[,1] > 0)
                                           tab_new = data.frame("beta" = tab[keep,1],
                                                                "absbeta" = abs(tab[keep,1]),
                                                                row.names = rownames(tab)[keep])
                                           colnames(tab_new) <- c("beta", "abs_beta")
                                           ml_tab <<- tab_new[order(tab_new[,1],decreasing = TRUE),]
                                           DT::renderDataTable({
                                             DT::datatable(ml_tab,
                                                           selection = 'single',
                                                           autoHideNavigation = T,
                                                           options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
                                           })
                                         })
                 
               }
             }, bar = {
               #
               #d = list(x=1)
               curr_cpd <<- as.character(global$tables$ml_bar_tab[d$x,"mz"][[1]])
               
               output$ml_specific_plot <- plotly::renderPlotly({
                 # --- ggplot ---
                 ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
               })
             })
           },
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
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           fc = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$fc$fc.log)) return(NULL)
             curr_cpd <<- d$key
             output$fc_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           aov = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$aov$p.value)) return(NULL)
             curr_cpd <<- d$key
             output$aov_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           rf = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(vip.score)) return(NULL)
             curr_cpd <<- d$key
             output$rf_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           heatmap_biv = {
             if(!exists("hm_matrix")) return(NULL)
             if(d$y > length(hm_matrix$matrix$rows)) return(NULL)
             curr_cpd <<- hm_matrix$matrix$rows[d$y]
             output$heatmap_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           heatmap_mult = {
             if(!exists("hm_matrix")) return(NULL)
             if(d$y > length(hm_matrix$matrix$rows)) return(NULL)
             curr_cpd <<- hm_matrix$matrix$rows[d$y]
             output$heatmap_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           volc = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% rownames(mSet$analSet$volcano$sig.mat)) return(NULL)
             curr_cpd <<- d$key
             output$volc_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           },
           lasnet = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$lasnet$cpds)) return(NULL)
             curr_cpd <<- d$key
             output$lasnet_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)
             })
           })
    # ----------------------------
    output$curr_cpd <- renderText(curr_cpd)
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
    req(global$vectors$db_search_list)
    # ----------------
    if(length(global$vectors$db_search_list) > 0){
      global$tables$last_matches <<- unique(multimatch(curr_cpd, global$vectors$db_search_list))
      output$match_tab <- DT::renderDataTable({
        DT::datatable(global$tables$last_matches[,-c("description","structure", "baseformula")],
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
    curr_def <<- global$tables$last_matches[curr_row,'description']
    output$curr_definition <- renderText(curr_def$description)
    curr_struct <<- global$tables$last_matches[curr_row,'structure'][[1]]
    output$curr_struct <- renderPlot({plot.mol(curr_struct,style = "cow")})
    curr_formula <<- global$tables$last_matches[curr_row,'baseformula'][[1]]
    output$curr_formula <- renderText({curr_formula})
  })
  
  observeEvent(input$browse_db,{
    #req(db_search_list)
    # -------------------
    cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
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
    #req(db_search_list)
    req(input$browse_tab_rows_selected)
    # -------------------
    search_cmd <- browse_table[curr_row,c('Formula', 'Charge')]
    # -------------------
    cpd_list <- lapply(global$vectors$db_search_list, FUN=function(match.table){
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
    #TODO: this should be a function and not re-written
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = color.vec(),cf=global$functions$color.function)})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
  })
  
  observeEvent(input$hits_tab_rows_selected,{
    curr_row = input$hits_tab_rows_selected
    curr_row <<- input$hits_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_cpd <<- hits_table[curr_row, mzmed.pgrp]
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = color.vec(),cf=global$functions$color.function)})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec(),cf=global$functions$color.function)})
  })
  
  observeEvent(input$go_enrich,{
    #TODO: fix, currently broken (build in m/z > cpd search?)
    enrich_db_list <- paste0("./backend/db/", c("kegg", 
                                                "wikipathways", 
                                                "smpdb",
                                                "metacyc"), ".full.db")
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
      shiny::setProgress(session=session, value= 0.33)
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
      shiny::setProgress(session=session, value= 0.66)
      
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
      shiny::setProgress(session=session, value= 1)
    })
  })
  
  
  # === META ANALYSIS ===
  
  # --- combine tab ---
  
  observeEvent(input$venn_members, {
    
    if("ls" %in% input$venn_members | "rf" %in% input$venn_members | "plsda" %in% input$venn_members){
      uiElements <- lapply(intersect(c("rf", "ls", "plsda"), input$venn_members), function(name){
        options <- switch(name,
                          tt = NULL,
                          fc = NULL,
                          ls = names(mSet$analSet$ml$ls),
                          rf = names(mSet$analSet$ml$rf),
                          volc = NULL,
                          plsda = c("PC1", "PC2", "PC3"))
        
        shorthand = switch(name,
                           ls = "lasso",
                           rf = "random forest",
                           plsda = "PLSDA")
        selectInput(inputId = paste0(name,"_choice"), 
                    label = gsubfn::fn$paste("Which $shorthand group(s):"), 
                    multiple=TRUE,
                    choices = options)
      })
      
      output$venn_ml_ui <- renderUI({uiElements})
    }else{
      output$venn_ml_ui <- renderUI({br()})
    }
  })
  
  observeEvent(input$build_venn, {
    
    # - - - cats - - -
    
    top = input$venn_tophits
    
    categories <- input$venn_members
    
    #top = 200
    #categories = c("tt", "fc")
    
    tables <- lapply(categories, function(name){
      
      tbls <- switch(name,
                     tt = list(rownames(mSet$analSet$tt$sig.mat[order(mSet$analSet$tt$sig.mat[,2], decreasing = F),])),
                     fc = list(rownames(mSet$analSet$fc$sig.mat[order(abs(mSet$analSet$fc$sig.mat[,2]), decreasing = F),])),
                     ls = {
                       tbls_ls <- lapply(input$ls_choice, function(name){
                         mSet$analSet$ml$ls[[name]][order(mSet$analSet$ml$ls[[name]]$count, decreasing = T),]$mz})
                       names(tbls_ls) <- input$ls_choice
                       # - - - 
                       tbls_ls
                     },
                     rf = {
                       tbls_rf <- lapply(input$rf_choice, function(name){
                         mSet$analSet$ml$rf[[name]][order(mSet$analSet$ml$rf[[name]]$mda, decreasing = T),]$mz})
                       names(tbls_rf) <- input$rf_choice
                       # - - -
                       tbls_rf
                     },
                     plsda = {
                       tbls_plsda <- lapply(input$plsda_choice, function(name){
                         compounds_pc <- as.data.table(mSet$analSet$plsda$vip.mat,keep.rownames = T)
                         colnames(compounds_pc) <- c("rn", paste0("PC", 1:3))
                         ordered_pc <- setorderv(compounds_pc, name, -1)
                         ordered_pc[, c("rn")][[1]]               
                       })
                       names(tbls_plsda) <- paste0("pls_", input$plsda_choice)
                       # - - -
                       tbls_plsda
                     },
                     volc = list(rownames(mSet$analSet$volcano$sig.mat)))
      
      # - - - - - - - 
      # FIX NAMES!!!
      
      tbls_top <- lapply(tbls, function(tbl){
        if(length(tbl) < top){
          tbl
        }else{
          tbl[1:top]
        }
      })
      
      # - - return - -
      tbls_top
    })
    
    names(tables) <- categories
    # - - unlist - -
    
    flattened <<- flattenlist(tables)
    names(flattened) <<- gsub(x = names(flattened), pattern = "(.*\\.)(.*$)", replacement = "\\2")
    
    circles = length(flattened)
    
    # - - - - - - --
    
    # TODO
    # Hypergeometric testing?
    # SO : 'Calculate venn diagram hypergeometric p value using R'
    
    venn.plot <- VennDiagram::venn.diagram(x = flattened,
                                           filename = NULL)
    
    # - - - - - - - - -
    
    items <- strsplit(as.character(venn.plot), split = ",")[[1]]
    
    circ_values <<- data.frame(
      id = 1:length(grep(items, pattern="polygon"))
      #,value = c(3, 3.1, 3.1, 3.2, 3.15, 3.5)
    )
    
    txt_values <- data.frame(
      id = grep(items, pattern="text"),
      value = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$label))
    )
    
    txt_values$value <- gsub(x = txt_values$value, pattern = "(.*\\.)(.*$)", replacement = "\\2")
    categories <- c(categories, input$rf_choice, input$ls_choice, input$plsda_choice)
    
    x_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$x))
    y_c = unlist(lapply(grep(items, pattern="polygon"), function(i) venn.plot[[i]]$y))
    
    x_t = unlist(lapply(grep(items, pattern="text"), function(i) venn.plot[[i]]$x))
    y_t = unlist(lapply(grep(items, pattern="text"), function(i)venn.plot[[i]]$y))
    
    positions_c <- data.frame(
      id = rep(circ_values$id, each = length(x_c)/length(circ_values$id)),
      x = x_c,
      y = y_c
    )
    
    positions_t <- data.frame(
      id = rep(txt_values$id, each = length(x_t)/length(txt_values$id)),
      x = x_t,
      y = y_t
    )
    
    datapoly <- merge(circ_values, positions_c, by=c("id"))
    datatxt <- merge(txt_values, positions_t, by=c("id"))
    
    numbers <- datatxt[!(datatxt$value %in% names(flattened)),]
    headers <- datatxt[(datatxt$value %in% names(flattened)),]
    
    if(circles == 2){
      occur <- table(numbers$y)
      newy <- names(occur[occur == max(occur)])
      # - - -
      numbers$y <- as.numeric(c(newy))
    }
    
    p <- ggplot(datapoly, 
                aes(x = x, 
                    y = y)) + geom_polygon(colour="black", alpha=0.5, aes(fill=id, group=id)) +
      geom_text(mapping = aes(x=x, y=y, label=value), data = numbers, size = 5) +
      geom_text(mapping = aes(x=x, y=y, label=value), data = headers, fontface="bold", size = 7) +
      theme_void() +
      theme(legend.position="none") + 
      scale_fill_gradientn(colours = global$functions$color.function(circles)) +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
    
    
    output$venn_plot <- plotly::renderPlotly({
      
      ggplotly(p, tooltip = "label")
      #%>% layout(plot_bgcolor='transparent')
      #%>% layout(paper_bgcolor='transparent')
    
      })
    updateSelectizeInput(session, "intersect_venn", choices = names(flattened))
  })
  
  observeEvent(input$intersect_venn, {

    if(length(input$intersect_venn) > 1){
      venn_overlap <<- Reduce("intersect", lapply(input$intersect_venn, function(x){
        flattened[[x]]
      })
      )
    }
    output$venn_tab <- DT::renderDataTable({
      # -------------
      DT::datatable(data.table(mz = venn_overlap), 
                    selection = 'single',
                    autoHideNavigation = T,
                    options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })
  })
  
  observeEvent(input$venn_tab_rows_selected, {
    curr_cpd <<- venn_overlap[input$venn_tab_rows_selected]
    output$curr_cpd <- renderText(curr_cpd)
  })
  
  # - - ON CLOSE - -
  
  observe({
    if (input$nav_general == "stop"){
      shinyalert::shinyalert(title = "Question", 
                             text = "Do you want to close metaboShiny?", 
                             type = "warning",
                             #imageUrl = "www/question.png", 
                             showCancelButton = T, 
                             cancelButtonText = "No",
                             showConfirmButton = T, 
                             confirmButtonText = "Yes",
                             callbackR = function(x) {if(x == TRUE){
                               print("closing metaboShiny ~ヾ(＾∇＾)")
                               if(any(!is.na(session_cl))) parallel::stopCluster(session_cl)
                               R.utils::gcDLLs() # flush dlls
                               stopApp() 
                             }else{
                               NULL
                             }})
    }
  })
})
