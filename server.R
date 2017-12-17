shinyServer(function(input, output, session) {
  
  # ================================= DEFAULTS ===================================
  
  tbl <<- NA
  tables <- list()
  db_search_list <- c()
  color.function <- rainbow
  patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
  shinyOptions(progress.style="old")
  ppm <<- 2
  nvars <<- 2
  
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
                            sardine(fadeImageButton("add_chebi", img.path = "chebilogo.png")),br(),
                            sardine(fadeImageButton("add_smpdb", img.path = "smpdb_logo_adj.png")),
                            sardine(fadeImageButton("add_wikipathways", img.path = "wikipathways.png")),
                            sardine(fadeImageButton("add_kegg", img.path = "kegglogo.gif", value = T))
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
                "wikipathways",
                "smpdb")
  
  images <<- list(list(name = 'cute_package', path = 'www/new-product.png', dimensions = c(80, 80)),
                  list(name = 'umc_logo_int', path = 'www/umcinternal.png', dimensions = c(120, 120)),
                  list(name = 'umc_logo_noise', path = 'www/umcnoise.png', dimensions = c(120, 120)),
                  list(name = 'hmdb_logo', path = 'www/hmdblogo.png', dimensions = c(150, 100)),
                  list(name = 'chebi_logo', path = 'www/chebilogo.png', dimensions = c(140, 140)),
                  list(name = 'wikipath_logo', path = 'www/wikipathways.png', dimensions = c(130, 150)),
                  list(name = 'kegg_logo', path = 'www/kegglogo.gif', dimensions = c(200, 150)),
                  list(name = 'pubchem_logo', path = 'www/pubchemlogo.png', dimensions = c(145, 90)),
                  list(name = 'smpdb_logo', path = 'www/smpdb_logo_adj.png', dimensions = c(200, 160)),
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
    navbarPage("Standard analysis", id="tab_stat",
               tabPanel("", value = "intro", icon=icon("comment-o"),
                        helpText("Info text here")
               ),
               tabPanel("PCA", value = "pca", #icon=icon("cube"),
                        plotly::plotlyOutput("plot_pca",height = "600px"),
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
               tabPanel("PLSDA", value = "plsda",
                        navbarPage("",
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
                        
               ), tabPanel("T-test", value="tt", 
                           fluidRow(plotly::plotlyOutput('tt_specific_plot',width="100%")),
                           navbarPage("",
                                      tabPanel("", icon=icon("table"),
                                               div(DT::dataTableOutput('tt_tab',width="100%"),style='font-size:80%'))
                                      ,tabPanel("", icon=icon("area-chart"),
                                                plotly::plotlyOutput('tt_overview_plot',height="300px")
                                      )
                           )),
               tabPanel("Fold-change", value="fc",
                        fluidRow(plotly::plotlyOutput('fc_specific_plot',width="100%")),
                        navbarPage("",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('fc_tab',width="100%"),style='font-size:80%'))
                                   ,tabPanel("", icon=icon("area-chart"),
                                             plotly::plotlyOutput('fc_overview_plot',height="300px")
                                   )
                        )),
               # =================================================================================
               tabPanel("RandomForest", value="rf",
                        fluidRow(plotly::plotlyOutput('rf_specific_plot',width="100%")),
                        navbarPage("",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('rf_tab',width="100%"),style='font-size:80%'))
                        )
               ),
               tabPanel("Heatmap", value="heatmap_biv",
                        plotly::plotlyOutput("heatmap",width="110%",height="700px"),
                        br(),
                        fluidRow(column(align="center",
                                        width=12,switchButton(inputId = "heatmode",
                                                              label = "Use data from:", 
                                                              value = TRUE, col = "BW", type = "TTFC"))
                        )
               ),
               tabPanel("Volcano", value="volc",
                        fluidRow(plotly::plotlyOutput('volc_plot',width="100%",height="600px"))),
               tabPanel("Enrichment", value="enrich_biv",
                        sidebarLayout(position="left",
                                      sidebarPanel = sidebarPanel(align="center",
                                                                  fluidRow(
                                                                    tags$b("Pathway DBs"),br(),
                                                                    sardine(fadeImageButton("enrich_smpdb", 
                                                                                            img.path = "smpdb_logo_adj.png", 
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
                                                            fluidRow(div(DT::dataTableOutput('enrich_pw_tab'),style='font-size:80%'))
                                      ))
               )
    )
  })
  
  stat.ui.multivar <- reactive({
    navbarPage("Standard analysis", id="tab_stat",
               tabPanel("", value = "intro", icon=icon("comment-o"),
                        helpText("Info text here")
               ),
               tabPanel("PCA", value = "pca", #icon=icon("cube"),
                        plotly::plotlyOutput("plot_pca",height = "600px"),
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
               tabPanel("ANOVA", value="aov",
                        fluidRow(plotly::plotlyOutput('aov_specific_plot',width="100%")),
                        navbarPage("",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('aov_tab',width="100%"),style='font-size:80%'))
                                   ,tabPanel("", icon=icon("area-chart"),
                                             plotly::plotlyOutput('aov_overview_plot',height="300px")
                                   )
                        )),
               # =================================================================================
               tabPanel("RandomForest", value="rf",
                        fluidRow(plotly::plotlyOutput('rf_specific_plot',width="100%")),
                        navbarPage("",
                                   tabPanel("", icon=icon("table"),
                                            div(DT::dataTableOutput('rf_tab',width="100%"),style='font-size:80%'))
                        )
               ),
               tabPanel("Heatmap", value="heatmap_mult",
                        plotly::plotlyOutput("heatmap",width="110%",height="700px")
               ),
               tabPanel("Enrichment", value="enrich_multi",
                        sidebarLayout(position="left",
                                      sidebarPanel = sidebarPanel(align="center",
                                                                  fluidRow(
                                                                    tags$b("Pathway DBs"),br(),
                                                                    sardine(fadeImageButton("enrich_smpdb", 
                                                                                            img.path = "smpdb_logo_adj.png", 
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
                                                            fluidRow(div(DT::dataTableOutput('enrich_pw_tab'),style='font-size:80%'))
                                      ))
               )
    )
    
  })
  
  
  time.ui <- reactive({
    navbarPage("Time Series", id="tab_time",
               tabPanel("iPCA", value = "ipca", 
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
               tabPanel("MEBA", value="meba", 
                        fluidRow(plotly::plotlyOutput('meba_specific_plot'),height="600px"),
                        fluidRow(div(DT::dataTableOutput('meba_tab', width="100%"),style='font-size:80%'))
               ),
               # =================================================================================
               tabPanel("ASCA", value="asca",
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
    print(given_dir)
    if(is.null(given_dir)) return()
    options$work_dir <<- given_dir
    output$exp_dir <- renderText(options$work_dir)
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
        showgrid = FALSE
      )
      # --- render ---
      plotly::plot_ly(z = volcano, 
                      colors = color.function(100), 
                      type = "heatmap",
                      showscale=FALSE)  %>%
        layout(xaxis = ax, yaxis = ax)
    })
  })
  
  
  color.pickers <- reactive({
    req(mSet$dataSet)
    # -------------
    switch(input$exp_type,
           stat = {
             facs <- levels(mSet$dataSet$filt.cls)
           },
           time_std = {
             lbl.fac <- if(mSet$dataSet$facA.lbl == "Time") "facB" else "facA"
             facs <- levels(mSet$dataSet[[lbl.fac]])
           },
           time_fin = {
             facs <- levels(mSet$dataSet$filt.cls)
           },
           time_min = {           
             facs <- levels(mSet$dataSet$filt.cls)
           })
    default.colours <- rainbow(length(facs))
    # -------------
    lapply(seq_along(facs), function(i) {
      colourInput(inputId = paste("col", i, sep="_"), 
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
             facs <- mSet$dataSet$filt.cls
           },
           time_std = {
             if("facA" %not in% names(mSet$dataSet) & "facB" %not in% names(mSet$dataSet)) return(c("Blue", "Pink"))
             lbl.fac <- if(mSet$dataSet$facA.lbl == "Time") "facB" else "facA"
             facs <- mSet$dataSet[[lbl.fac]]
           },
           time_fin = {
             facs <- mSet$dataSet$filt.cls
           },
           time_min = {           
             facs <- mSet$dataSet$filt.cls
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
        DF = wkz.adduct.confirmed
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
      withProgress({
        parallel::clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
          "isotopes",
          "subform.joanna", 
          "mergeform.joanna",
          "multiform.joanna",
          "check.ded.joanna",
          "data.table",
          "rbindlist",
          "isopattern",
          "keggFind",
          "keggGet",
          "kegg.charge",
          "regexpr",
          "regmatches"
        ))
        #setProgress(message = "Working...")
        build.base.db(db, outfolder=options$db_dir)
        setProgress(0.5)
        build.extended.db(db, 
                          outfolder=options$db_dir,
                          adduct.table = wkz.adduct.confirmed, 
                          cl=session_cl, 
                          fetch.limit=100)
      })
    })
  })
  
  # ================== DATA IMPORT ===========================
  
  observeEvent(input$create_db,{
    # --------------------
    patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
    # --------------------
    withProgress({
      setProgress(.25,message = "Loading outlists into memory...")
      req(input$outlist_neg, input$outlist_pos, input$excel)
      
      outlist_pos <- as.data.table(fread(file = input$outlist_pos$datapath, sep="\t", header=T))
      outlist_neg <- as.data.table(fread(file = input$outlist_neg$datapath, sep="\t", header=T))
      
      setProgress(.50,message = "Creating experiment database file...")
      
      build.pat.db(patdb,
                   ppm = ppm,
                   poslist = outlist_pos,
                   neglist = outlist_neg,
                   overwrite = T)
      
      setProgress(.75,message = "Adding excel sheets to database...")
      
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
  
  # ==================== CREATE CSV =======================
  
  observeEvent(input$create_csv, {
    req(options$proj_name)
    req(options$work_dir)
    # ---------
    withProgress({
      setProgress(1/4)
      # create csv
      print(if(input$broadvars) "individual_data" else "setup")
      tbl <- get.csv(patdb,
                     time.series = if(input$exp_type == "time_std") T else F,
                     exp.condition = input$exp_var,
                     group_adducts = if(length(add_search_list) == 0) F else T,
                     group_by = input$group_by,
                     which_dbs = add_search_list,
                     which_adducts = selected_adduct_list,
                     var_table = if(input$broadvars) "individual_data" else "setup"
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
      export_tbl <<- tbl.adj
      # save csv
      setProgress(2/4)
      csv_loc <<- file.path(options$work_dir, paste0(options$proj_name,".csv"))
      fwrite(tbl.adj, csv_loc, sep="\t")
      # --- overview table ---
      setProgress(3/4)
      output$csv_tab <-DT::renderDataTable({
        overview_tab <- if(input$exp_type == "time_std"){
          t(data.table(keep.rownames = F,
                       Identifiers = ncol(tbl.adj) - 3,
                       Samples = nrow(tbl.adj),
                       Times = length(unique(tbl.adj$Time)),
                       Groups = length(unique(tbl.adj$Label))
          ))
        } else{
          t(data.table(keep.rownames = F,
                       Identifiers = ncol(tbl.adj) - 3,
                       Samples = nrow(tbl.adj),
                       Groups = length(unique(tbl.adj$Label))
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
  
  observeEvent(input$check_excel, {
    # get excel table stuff.
    updateSelectInput(session, "exp_var",
                      choices = get_exp_vars(if(input$broadvars) "individual_data" else "setup")
    )
  })
  
  
  ref.selector <- reactive({
    # -------------
    if(input$norm_type == "ProbNorm"){
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
    # get excel table stuff.
    updateSelectInput(session, "ref_var",
                      choices = get_ref_vars(fac = "Label") # please add options for different times later, not difficult
    )
  })
  
  output$ref_select <- renderUI({ref.selector()})
  
  observeEvent(input$initialize,{
    req(input$filt_type)
    req(input$norm_type)
    req(input$trans_type)
    req(input$scale_type)
    # match
    withProgress({
      # get curr values from: input$ exp_type, filt_type, norm_type, scale_type, trans_type (later)
      setProgress(.1)
      # ------------------------------
      #Below is your R command history: 
      switch(input$exp_type,
             stat = {
               mSet <<- InitDataObjects("pktable", 
                                        "stat", 
                                        FALSE)
               mSet <<- Read.TextData(mSet, 
                                      filePath = csv_loc,
                                      "rowu")
               mSet
             },
             time_std = {
               mSet <<- InitDataObjects("pktable", 
                                        "ts", 
                                        FALSE)
               mSet <<- Read.TextData(mSet, 
                                      filePath = csv_loc,
                                      "rowu")
               mSet
             },
             time_fin = {
               mSet <<- InitDataObjects("pktable", 
                                        "stat", 
                                        FALSE)
               mSet <<- Read.TextData(mSet, 
                                      filePath = csv_loc,
                                      "rowu")
               mSet
             },
             time_min = {   
               mSet <<- InitDataObjects("pktable", 
                                        "stat", 
                                        FALSE)
               mSet <<- Read.TextData(mSet, 
                                      filePath = csv_loc,
                                      "rowu")
               mSet
             }
      )
      
      mSet <<- SanityCheckData(mSet)
      mSet <<- RemoveMissingPercent(mSet, percent = 0.5)
      mSet <<- ImputeVar(mSet,
                         method = "min")
      mSet <<- ReplaceMin(mSet)
      mSet <<- FilterVariable(mSet,
                              filter = input$filt_type,
                              qcFilter = "F", 
                              rsd = 25)
      # ------------------------------------
      setProgress(.2)
      mSet <<- Normalization(mSet,
                             rowNorm = input$norm_type,
                             transNorm = input$trans_type,
                             scaleNorm = input$scale_type,
                             ref = input$ref_var)
      mSet$dataSet$grouping <<- input$group_by
      # -------------------------------------
      if(mSet$dataSet$grouping == "pathway"){mSet$dataSet$norm <- abs(mSet$dataSet$norm)}
      setProgress(.3)
      varNormPlots <- ggplotNormSummary(mSet)
      output$var1 <- renderPlot(varNormPlots$tl)
      output$var2 <- renderPlot(varNormPlots$bl)
      output$var3 <- renderPlot(varNormPlots$tr)
      output$var4 <- renderPlot(varNormPlots$br)
      
      setProgress(.4)
      sampNormPlots <-  ggplotSampleNormSummary(mSet)
      output$samp1 <- renderPlot(sampNormPlots$tl)
      output$samp2 <- renderPlot(sampNormPlots$bl)
      output$samp3 <- renderPlot(sampNormPlots$tr)
      output$samp4 <- renderPlot(sampNormPlots$br)
      
      mSet$dataSet$adducts <<- selected_adduct_list
      update.UI()  
    })
  })
  
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
               print(fac.lvls)
               chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
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
                 print(row)
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
                   opacity=0.1
                 )
               }
               adj_plot <<- plotly_build(plots)
               rgbcols <- toRGB(chosen.colors)
               c = 1
               for(i in seq_along(adj_plot$x$data)){
                 print(i)
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
                   xaxis = list(
                     title = mSet$analSet$ipca$score$axis[x.num]),
                   yaxis = list(
                     title = mSet$analSet$ipca$score$axis[y.num]),
                   zaxis = list(
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
  
  observeEvent(input$tab_stat, {
    if(input$nav_general != "analysis") return(NULL)
    nvars <<- length(levels(mSet$dataSet$cls))
    # get excel table stuff.
    switch(input$tab_stat,
           pca = {
             mSet <<- PCA.Anal(mSet)
             output$plot_pca <- plotly::renderPlotly({
               df <- mSet$analSet$pca$x
               x <- input$pca_x
               y <- input$pca_y
               z <- input$pca_z
               #x =1;y=2;z=3
               x.var <- round(mSet$analSet$pca$variance[x] * 100.00, digits=1)
               y.var <- round(mSet$analSet$pca$variance[y] * 100.00, digits=1)
               z.var <- round(mSet$analSet$pca$variance[z] * 100.00, digits=1)
               fac.lvls <- unique(mSet$dataSet$filt.cls)
               chosen.colors <<- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
               # --- add ellipses ---
               classes <- mSet$dataSet$cls
               plots <- plotly::plot_ly()
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
                   opacity=0.1
                 )
               }
               adj_plot <<- plotly_build(plots)
               rgbcols <- toRGB(chosen.colors)
               c = 1
               for(i in seq_along(adj_plot$x$data)){
                 print(i)
                 item = adj_plot$x$data[[i]]
                 if(item$type == "mesh3d"){
                   adj_plot$x$data[[i]]$color <- rgbcols[c]
                   c = c + 1
                 }
               }
               # --- return ---
               pca_plot <- adj_plot %>% add_trace(
                 hoverinfo = 'text',
                 text = rownames(df),
                 x = mSet$analSet$pca$x[,1], 
                 y = mSet$analSet$pca$x[,2], 
                 z = mSet$analSet$pca$x[,3], 
                 type = "scatter3d",
                 opacity=1,
                 color= mSet$dataSet$filt.cls, colors=chosen.colors
               ) %>%  layout(scene = list(
                 xaxis = list(
                   title = fn$paste("$x ($x.var %)")),
                 yaxis = list(
                   title = fn$paste("$y ($y.var %)")),
                 zaxis = list(
                   title = fn$paste("$z ($z.var %)"))))     
               # --- return ---
               pca_plot
             })
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
           plsda = {
             mSet <<- PLSR.Anal(mSet,
                                F)
             mSet <<- PLSDA.CV(mSet,
                               "L",
                               5, 
                               "Q2")
             plsda.table <- as.data.table(round(mSet$analSet$plsr$Xvar 
                                                / mSet$analSet$plsr$Xtotvar 
                                                * 100.0,
                                                digits = 2),
                                          keep.rownames = T)
             colnames(plsda.table) <- c("Principal Component", "% variance")
             plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
             # -----------
             output$plot_plsda <- plotly::renderPlotly({
               x <- input$plsda_x
               y <- input$plsda_y
               z <- input$plsda_z
               x.var <- plsda.table[`Principal Component` == x, 
                                    `% variance`]
               y.var <- plsda.table[`Principal Component` == y, 
                                    `% variance`]
               z.var <- plsda.table[`Principal Component` == z, 
                                    `% variance`]
               fac.lvls <- unique(mSet$dataSet$filt.cls)
               colnames(mSet$analSet$plsr$Yscores) <- paste0("PC", 1:ncol(mSet$analSet$plsr$Yscores))
               chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
               # ---------------
               plotly::plot_ly(hoverinfo = 'text',
                               text = rownames(mSet$dataSet$norm)) %>%
                 add_trace(
                   x = mSet$analSet$plsr$Yscores[,x], 
                   y = mSet$analSet$plsr$Yscores[,y], 
                   z = mSet$analSet$plsr$Yscores[,z], 
                   type = "scatter3d",
                   color = mSet$dataSet$filt.cls, colors=chosen.colors
                 ) %>%  layout(scene = list(
                   xaxis = list(
                     title = fn$paste("$x ($x.var %)")),
                   yaxis = list(
                     title = fn$paste("$y ($y.var %)")),
                   zaxis = list(
                     title = fn$paste("$z ($z.var %)")
                   )
                 ))
             })
             output$plsda_tab <-DT::renderDataTable({
               req(plsda.table)
               # -------------
               DT::datatable(plsda.table, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
             # -------------
             output$plsda_vip_tab <-DT::renderDataTable({
               colnames(mSet$analSet$plsda$vip.mat) <- paste0("PC", 1:ncol(mSet$analSet$plsda$vip.mat))
               compounds_pc <- as.data.table(mSet$analSet$plsda$vip.mat,keep.rownames = T)
               ordered_pc <- setorderv(compounds_pc, input$plsda_vip_cmp, -1)
               plsda_tab <<- cbind(ordered_pc[1:50, c("rn")], 
                                   Rank=c(1:50))
               # -------------
               DT::datatable(plsda_tab,
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             })
           },
           heatmap_mult = {
             req(input$color_ramp)
             req(mSet$analSet$aov)
             used.analysis <- "aov"
             used.values <- "p.value"
             sourceTable <- sort(mSet$analSet[[used.analysis]][[used.values]])[1:100]
             
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
                                    col_side_palette = function(n){
                                      if(n == length(color.vec())) color.vec() else rainbow(n)
                                      rainbow(n)
                                    },
                                    subplot_widths = c(.9,.1),
                                    subplot_heights =  c(.1,.05,.85),
                                    column_text_angle = 90)
             })
           },
           heatmap_biv = {
             req(input$color_ramp)
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
                                    col_side_palette = function(n){
                                      if(n == length(color.vec())) color.vec() else rainbow(n)
                                      rainbow(n)
                                    },
                                    subplot_widths = c(.9,.1),
                                    subplot_heights =  c(.1,.05,.85),
                                    column_text_angle = 90)
             })
           },
           tt = {
             mSet <<- Ttests.Anal(mSet,
                                  nonpar = F, 
                                  threshp = 0.1, 
                                  paired = F,
                                  equal.var = T)
             res <- mSet$analSet$tt$sig.mat
             if(is.null(res)) res <- data.table("No significant hits found")
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
             mSet <<- FC.Anal.unpaired(mSet,
                                       2.0, 
                                       1)
             res <- mSet$analSet$fc$sig.mat
             if(is.null(res)) res <- data.table("No significant hits found")
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
             mSet <<- ANOVA.Anal(mSet, 
                                 thresh=0.05,
                                 nonpar = F)
             res <- mSet$analSet$aov$sig.mat
             if(is.null(res)) res <- data.table("No significant hits found")
             output$aov_tab <-DT::renderDataTable({
               # -------------
               DT::datatable(res, 
                             selection = 'single',
                             autoHideNavigation = T,
                             options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
               
             })
           },
           rf={
             mSet <<- RF.Anal(mSet, 500,7,1)
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
             mSet <<- Volcano.Anal(mSet,FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
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
    print(enrich_db_list)
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
                          "enrich_pw")
  
  observeEvent(input$enriched_rows_selected, {
    curr_row = input$enriched_rows_selected
    print(curr_row)
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
  
  lapply(res.update.tables, FUN=function(table){
    observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
      curr_row = input[[paste0(table, "_tab_rows_selected")]]
      # do nothing if not clicked yet, or the clicked cell is not in the 1st column
      if (is.null(curr_row)) return()
      curr_cpd <<- data.table(switch(table,
                                     tt = mSet$analSet$tt$sig.mat,
                                     fc = mSet$analSet$fc$sig.mat,
                                     asca = mSet$analSet$asca$sig.list$Model.ab,
                                     aov = mSet$analSet$aov$sig.mat,
                                     rf = vip.score,
                                     enrich_pw = enrich_overview_tab,
                                     meba = mSet$analSet$MB$stats,
                                     plsda_vip = plsda_tab)
                              , keep.rownames = T)[curr_row, rn]
      output$curr_cpd <- renderText(curr_cpd)
      output[[paste0(table, "_specific_plot")]] <- plotly::renderPlotly({
        # --- ggplot ---
        if(table == 'meba'){
          ggplotMeba(curr_cpd, draw.average = T, cols = color.vec())
        }else{
          ggplotSummary(curr_cpd, cols = color.vec())
        }
      })
      if(input$autosearch & length(db_search_list > 0)){
        match_list <- lapply(db_search_list, FUN=function(match.table){
          get_matches(curr_cpd, match.table, searchid=mSet$dataSet$grouping)
        })
        match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
        output$match_tab <-DT::renderDataTable({
          DT::datatable(if(mSet$dataSet$grouping == 'pathway') match_table else{match_table[,-"Description"]},
                        selection = 'single',
                        autoHideNavigation = T,
                        options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
        })  
      }
    })
  })
  
  observe({
    d <- plotly::event_data("plotly_click")
    req(d)
    print(d)
    # -----------------------------
    switch(input$tab_stat,
           tt = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$tt$p.value)) return(NULL)
             curr_cpd <<- d$key
             output$tt_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec())
             })
           },
           fc = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$fc$fc.log)) return(NULL)
             curr_cpd <<- d$key
             output$fc_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec())
             })
           },
           aov = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(mSet$analSet$aov$p.value)) return(NULL)
             curr_cpd <<- d$key
             output$aov_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec())
             })
           },
           rf = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% names(vip.score)) return(NULL)
             curr_cpd <<- d$key
             output$rf_specific_plot <- plotly::renderPlotly({
               # --- ggplot ---
               ggplotSummary(curr_cpd, cols = color.vec())
             })
           },
           heatmap = {
             if(!exists("hm_matrix")) return(NULL)
             if(d$y > length(hm_matrix$matrix$rows)) return(NULL)
             curr_cpd <<- hm_matrix$matrix$rows[d$y]
           },
           volc = {
             if('key' %not in% colnames(d)) return(NULL)
             if(d$key %not in% rownames(mSet$analSet$volcano$sig.mat)) return(NULL)
             curr_cpd <<- d$key
           })
    # ----------------------------
    output$curr_cpd <- renderText(curr_cpd)
    if(input$autosearch & length(db_search_list) > 0){
      match_list <- lapply(db_search_list, FUN=function(match.table){
        get_matches(curr_cpd, match.table, searchid=input$group_by)
      })
      match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
      output$match_tab <-DT::renderDataTable({
        DT::datatable(if(mSet$dataSet$grouping == 'pathway') match_table else{match_table[,-"Description"]},
                      selection = 'single',
                      autoHideNavigation = T,
                      options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })  
    }
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
      match_list <- lapply(db_search_list, FUN=function(match.table){
        get_matches(curr_cpd, match.table, searchid=input$group_by)
      })
      print(match_list)
      match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
      output$match_tab <-DT::renderDataTable({
        DT::datatable(if(mSet$dataSet$grouping == 'pathway') match_table else{match_table[,-"Description"]},
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
    curr_def <<- if(mSet$dataSet$grouping == 'pathway') "Unavailable" else(match_table[curr_row,'Description'])
    output$curr_definition <- renderText(curr_def$Description)
  })
  
  observeEvent(input$browse_db,{
    req(input$checkGroup)
    # -------------------
    cpd_list <- lapply(input$checkGroup, FUN=function(match.table){
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
    req(input$checkGroup)
    req(input$browse_tab_rows_selected)
    # -------------------
    search_cmd <- browse_table[curr_row,c('Formula', 'Charge')]
    # -------------------
    cpd_list <- lapply(input$checkGroup, FUN=function(match.table){
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
    curr_row <<- input$browse_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_def <<- browse_table[curr_row,'Description']
    output$browse_definition <- renderText(curr_def$Description)
  })
  
  observeEvent(input$hits_tab_rows_selected,{
    curr_row = input$hits_tab_rows_selected
    curr_row <<- input$hits_tab_rows_selected
    if (is.null(curr_row)) return()
    # -----------------------------
    curr_cpd <<- hits_table[curr_row, mzmed.pgrp]
    output$meba_specific_plot <- plotly::renderPlotly({ggplotMeba(curr_cpd, draw.average=T, cols = color.vec())})
    output$asca_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
    output$fc_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
    output$tt_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
    output$aov_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
    output$plsda_specific_plot <- plotly::renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
  })
  
  observeEvent(input$go_enrich,{
    withProgress({
      all_pathways <- lapply(enrich_db_list, FUN=function(db){
        conn <- dbConnect(RSQLite::SQLite(), db) # change this to proper var later
        dbname <- gsub(basename(db), pattern = "\\.full\\.db", replacement = "")
        dbGetQuery(conn, "SELECT distinct * FROM pathways limit 10;")
        gset <- dbGetQuery(conn,"SELECT DISTINCT c.baseformula AS cpd,
                           p.name as name
                           FROM pathways p
                           JOIN base c
                           ON c.pathway = p.identifier 
                           ")
        dbDisconnect(conn)
        # --- returny ---
        gset$name <- paste(gset$name, " (", dbname, ")", sep="")
        gset
      })
      head(all_pathways)
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


