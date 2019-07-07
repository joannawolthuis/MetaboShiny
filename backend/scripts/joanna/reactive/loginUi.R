output$currUI <- renderUI({
  if(logged$status == "notlogged"){
    fluidRow(align="center",
             br(),br(),br(),br(),br(),
             imageOutput("login_header",inline = T),
             textInput("username", "username:"),
             passwordInput("password", "password:"),
             tags$style(type="text/css", "#string { height: 50px; width: 100%; text-align:center;
                          font-size: 30px; display: block;}"),
             shinyWidgets::circleButton("login", icon = icon("arrow-right")),
             br(),br(),
             div(style="width:300px;",verbatimTextOutput("login_status", placeholder = FALSE)))
    
  }else if(logged$status == "logged"){
    
    print("rendering...")
    # read settings
    
    opts <- getOptions(lcl$paths$opt.loc)
    
    logged$text <<- "loaded options!"
    logged$text <<- "starting MetaboShiny..."
    
    
    # generate CSS for the interface based on user settings for colours, fonts etc.
    bar.css <<- nav.bar.css(opts$col1, opts$col2, opts$col3, opts$col4)
    font.css <<- app.font.css(opts$font1, opts$font2, opts$font3, opts$font4,
                              opts$size1, opts$size2, opts$size3, opts$size4, online=online)
    
    # google fonts
    
    # name = the name of the font in Google's Library (https://fonts.google.com)
    # family = how you want to refer to the font inside R
    # regular.wt = the weight of font used in axis, labels etc.
    # bolt.wt = the weight of font used in the title
    
    # === GOOGLE FONT SUPPORT FOR GGPLOT2 ===
    
    online = internetWorks()

    # Download a webfont
    if(online){
      lapply(c(opts[grepl(pattern = "font", names(opts))]), function(font){
        try({
          sysfonts::font_add_google(name = font,
                                    family = font,
                                    regular.wt = 400,
                                    bold.wt = 700)
        })
      })
    }
    
    # Perhaps the only tricky bit is remembering to run the following function to enable webfonts
    showtext::showtext_auto(enable = T)
    
    # ======================================
    
    # set taskbar image as set in options
    taskbar_image <- opts$task_img
    
    # parse color opts
    lcl$aes$mycols <<- get.col.map(lcl$paths$opt.loc) # colours for discrete sets, like group A vs group B etc.
    lcl$aes$theme <<- opts$gtheme # gradient function for heatmaps, volcano plot etc.
    lcl$aes$spectrum <<- opts$gspec # gradient function for heatmaps, volcano plot etc.

    # load existing file
    bgcol <<- opts$col1
    
    # - - load custom dbs - -
    
    # load in custom databases
    has.customs <- dir.exists(file.path(lcl$paths$db_dir, "custom"))
    
    if(has.customs){
      
      customs = list.files(path = file.path(lcl$paths$db_dir, "custom"),
                           pattern = "\\.RData")
      
      dbnames = unique(tools::file_path_sans_ext(customs))
      
      for(db in dbnames){
        # add name to global
        dblist <- gbl$vectors$db_list
        dblist <- dblist[-which(dblist == "custom")]
        if(!(db %in% dblist)){
          dblist <- c(dblist, db, "custom")
          gbl$vectors$db_list <- dblist
        }
        metadata.path <- file.path(lcl$paths$db_dir, "custom", paste0(db, ".RData"))
        load(metadata.path)
        
        # add description to global
        gbl$constants$db.build.info[[db]] <- meta.dbpage
        
        # add image to global
        maxi = length(gbl$constants$images)
        gbl$constants$images[[maxi + 1]] <- meta.img
      }
    }
    
    # init stuff that depends on opts file
    
    lcl$proj_name <- opts$proj_name
    lcl$paths$patdb <- file.path(opts$work_dir, paste0(opts$proj_name, ".db"))
    lcl$paths$csv_loc <- file.path(opts$work_dir, paste0(opts$proj_name, ".csv"))
    lcl$texts <- list(
      list(name='curr_exp_dir',text=lcl$paths$work_dir),
      list(name='curr_db_dir',text=lcl$paths$db_dir),
      list(name='ppm',text=opts$ppm),
      list(name='proj_name',text=opts$proj_name)
    )
    
    a = c("fish_no_out.csv", "fish.csv", "lungcancer_no_out.csv", "lungcancer.csv")
    #gsub(a,pattern = "(_no_out\\.csv)|(\\.csv)", replacement="")
    #dput(list.files(opts$work_dir,pattern = "\\.csv"))
    lcl$vectors$project_names <- unique(gsub(list.files(opts$work_dir,pattern = "\\.csv"),pattern = "(_no_out\\.csv)|(\\.csv)", replacement=""))
    print(lcl$vectors$proj_names)
    
    updateSelectizeInput(session,
                         "proj_name",
                         choices = lcl$vectors$proj_names,
                         selected = opts$proj_name)
    # create default text objects in UI
    lapply(lcl$texts, FUN=function(default){
      output[[default$name]] = renderText(default$text)
    })
    
    lcl$aes$font <- list(family = opts$font4,
                          ax.num.size = 11,
                          ax.txt.size = 15,
                          ann.size = 20,
                          title.size = 25)
    
    # other default stuff that needs opts
    
    library(showtext)
    
    online = internetWorks()
    
    # import google fonts
    for(font in unlist(opts[grep(names(opts), pattern = "font")])){
      if(font %in% sysfonts::font.families()){
        NULL
      }else{
        if(online) sysfonts::font_add_google(font,db_cache = T)
      }
    }
    
    # other stuff
    
    # create color pickers based on amount of colours allowed in global
    output$colorPickers <- renderUI({
      lapply(c(1:gbl$constants$max.cols), function(i) {
        colourpicker::colourInput(inputId = paste("col", i, sep="_"),
                                  label = paste("Choose colour", i),
                                  value = lcl$aes$mycols[i],
                                  allowTransparent = F)
      })
    })
    
    # create color1, color2 etc variables to use in plotting functions
    # and update when colours picked change
    observe({
      values <- unlist(lapply(c(1:gbl$constants$max.cols), function(i) {
        input[[paste("col", i, sep="_")]]
      }))
      
      if(!any(is.null(values))){
        if(lcl$paths$opt.loc != ""){
          set.col.map(optionfile = lcl$paths$opt.loc, values)
          lcl$aes$mycols <- values
        }
      }
    })
    
    updateSelectInput(session, "ggplot_theme", selected = opts$gtheme)
    updateSelectInput(session, "color_ramp", selected = opts$gspec)
    
    opts <<- opts
    # logged in!
    
    
    require(shinyjs)
    titlejs=paste0("document.title ='-`* MetaboShiny *`-'")
    runjs(titlejs)
    
    # - - - - -
    
    tagList(
      tags$title('MetaboShiny'),
      tags$style(type="text/css", bar.css),
      navbarPage(inverse=TRUE,#tags$head(tags$script(src="sparkle.js")),
                 title=div(h1("MetaboShiny"), class="outlined", tags$style(type="text/css", font.css)), # make it use the sparkle.js for unnecessary sparkle effects ;)
                 id="nav_general",
                 # this tab shows the available databases, if they are installed, and buttons to install them. generated as output$db_build_ui in 'server'
                 tabPanel("database", icon = icon("database",class = "outlined"), value="database",
                          uiOutput("db_build_ui") #%>% shinycssloaders::withSpinner() # see server, is autogenerated now
                 ),
                 tabPanel("data import", icon = icon("upload", class = "outlined"),
                          fluidRow(column(12, align="center", 
                                          textInput("proj_name_new", label = "STEP 1: What is your project name?", value = lcl$proj_name))),
                          hr(),
                          fluidRow(column(3, align="center",
                                          imageOutput("merge_icon",inline = T),
                                          radioButtons("importmode", label = "", 
                                                       choices = list("Peaks are in a .db file"="db", "Peaks are in two .csv files (pos/neg mode)"="csv"),
                                                       selected = "db"),
                                          tags$b("STEP 2: Click buttons to load data."),
                                          conditionalPanel(condition = "input.importmode == 'db'",
                                                           shinyFilesButton('database', 'Database', 'Select .db file', FALSE),
                                                           shinyFilesButton('metadata', 'Metadata', 'Select metadata in csv/xls(x)', FALSE)
                                                           ),
                                          conditionalPanel(condition = "input.importmode == 'csv'",
                                                           shinyFilesButton('outlist_pos', '+ peaks', 'Select .csv for - mode peaks', FALSE),
                                                           shinyFilesButton('metadata', 'Metadata', 'Select metadata in csv/xls(x)', FALSE),
                                                           shinyFilesButton('outlist_neg', '- peaks', 'Select .csv for + mode peaks', FALSE)
                                                           )
                                          )
                                   ,column(2, align="center", #ok
                                           tags$b("STEP 3: Merge data and metadata"),br(),br(),
                                             shinyWidgets::circleButton("create_db", icon = icon("long-arrow-alt-right", class = "fa-2x"), size = "lg"))
                                   ,column(2, align="center", # issue lol
                                          imageOutput("db_icon")
                                          )
                                   ,column(2, align="center",
                                           tags$b("STEP 4: Convert to input-ready format"),
                                           br(),br(),
                                          shinyWidgets::circleButton("create_csv", icon = icon("long-arrow-alt-right", class = "fa-2x"), size = "lg"))
                                   ,column(3, align="center",
                                          imageOutput("laptop_icon", inline=T),br(),br(),
                                          div(DT::dataTableOutput('csv_tab'),style='font-size:80%')
                                          )
                                   ),
                          fluidRow(column(3, align="center",
                                          tags$i("Input files chosen?"),br(),br(),
                                          imageOutput("proj_merge_check")
                                          ),
                                   column(2, align="center",
                                          tags$i("Database present?"),br(),br(),
                                          imageOutput("proj_db_check"),offset = 2),
                                   column(3, align="center",
                                          tags$i("Final table present?"),br(),br(),
                                          imageOutput("proj_csv_check", inline=T),br(),br(),
                                          tags$b("STEP 5: If "), icon("check-circle"), tags$b(" continue to normalization"),
                                          offset = 2))
                          ),
                 # this tab shows the options for creating a new project, either from 2 csv files and an excel, or from a .db file and an excel.
                 # tabPanel("data import", icon = icon("upload",class = "outlined"), value="link",
                 #          fluidRow(column(12, align="center",
                 #                          h1("Create project"))),
                 #          fluidRow(column(6, 
                 #                          tabsetPanel(id="new_proj",selected = "From DB",
                 #                                      # first tab is the 'new' method that uses a .db file. generally used for very large projects by powerusers
                 #                                      # may want to hide later..
                 #                                      tabPanel(id="db", title="From DB",
                 #                                               br(),br(),
                 #                                               fluidRow(
                 #                                                 column(6, align="center",
                 #                                                        #imageOutput("db_icon", inline = T),
                 #                                                        br(),br()
                 #                                                 ),
                 #                                                 column(6, align="center",
                 #                                                        br(),br(),
                 #                                                        imageOutput("excel_icon",inline = T)
                 #                                                 )
                 #                                               ),
                 #                                               fluidRow(
                 #                                                 column(6, align="center",
                 #                                                        shinyFilesButton('database', 'DATABASE', 'Please select a database file', FALSE)
                 #                                                 ),
                 #                                                 column(6, align="center",
                 #                                                        shinyFilesButton('metadata', 'METADATA', 'Please select a metadata file', FALSE)
                 #                                                 )
                 #                                               )
                 #                                      ),
                 #                                      # second tab is the 'old' method that uses 2 csv files (with columns: "mzmed", "Sample1", "Sample2" etc)
                 #                                      tabPanel(id="csv", title="From CSV",
                 #                                               br(),br(),
                 #                                               fluidRow(column(4,  align="center",
                 #                                                               imageOutput("pos_icon",inline = T),
                 #                                                               br(),br(),
                 #                                                               shinyFilesButton('outlist_pos', '+ mode peaks', 'Please select a csv file', FALSE)
                 #                                               ),
                 #                                               column(4,  align="center",
                 #                                                      imageOutput("excel_icon_2",inline = T),
                 #                                                      br(),br(),
                 #                                                      shinyFilesButton('metadata', 'METADATA', 'Please select a metadata file', FALSE),
                 #                                                      hr(),
                 #                                                      sliderInput("ppm",label = "m/z accuracy",
                 #                                                                  min = 1, max = 50,
                 #                                                                  value = 2, post = " ppm")
                 #                                               ),
                 #                                               column(4,  align="center",
                 #                                                      imageOutput("neg_icon",inline = T),
                 #                                                      br(),br(),
                 #                                                      shinyFilesButton('outlist_neg', '- mode peaks', 'Please select a csv file', FALSE)
                 #                                               )
                 #                                               )
                 #                                      )
                 #                          ),
                 #                          hr(),
                 #                          fluidRow( align="center",
                 #                                    column(9,
                 #                                           textInput("proj_name_new", label = "Project name:", value = "my_metshi"),
                 #                                           shinyWidgets::circleButton("create_db", "Go", icon = icon("arrow-right"),size = "lg")
                 #                                    ))),
                 #                   column(6, align="center",
                 #                          fluidRow(align="center", h3("Generate CSV from database")),
                 #                          br(),
                 #                          shinyWidgets::circleButton("create_csv", icon=icon("arrow-right"), size = "lg"),
                 #                          hr()
                 #                          ))
                 # ),
                 # this tab is used to perform normalization of your data. settings are processed as input$filt_type etc. in 'server'.
                 tabPanel("normalize",  icon = icon("shower",class = "outlined"), value="filter",
                          fluidRow(column(3, aligh="center",
                                          selectInput('samp_var', 'Which variable represents sample amount/concentration?', choices = c("")), #TODO: only show this when normalize by sample specific factor (specnorm) is selected
                                          selectizeInput('batch_var', 'What are your batch variables?', choices = c("batch"), multiple=TRUE, options = list(maxItems = 2L)),
                                          actionButton("check_csv", "Get options", icon=icon("refresh")),
                                          hr(),
                                          shinyWidgets::sliderTextInput("perc_limit","Max. missing feature percent:",
                                                                        choices=c(0, 0.0001, 0.001, 0.01, 0.1, seq(1, 100, 1)),
                                                                        selected=1, grid = T),
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
                                          selectInput('scale_type', 'How will you scale your data?', choices = list("Autoscale/Z-transform"="AutoNorm",
                                                                                                                    "Mean-center"="MeanCenter",
                                                                                                                    "Pareto Scaling"="ParetoNorm",
                                                                                                                    "Range scaling"="RangeNorm",
                                                                                                                    "None"="NULL")),
                                          selectInput('miss_type', 'How to deal with missing values?', choices = list("Half feature minimum"="colmin",
                                                                                                                      "Half sample minimum"="rowmin",
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
                                                       value = FALSE, col = "BW", type = "YN"),
                                          actionButton("initialize", "Go", icon=icon("hand-o-right")),
                                          hr(),
                                          imageOutput("dataset_icon",inline = T),
                                          fileInput("pat_dataset", "Import dataset",
                                                    multiple = F,
                                                    accept = c(".RData")),
                                          actionButton("import_dataset", "Import", icon = icon("hand-peace-o")),
                                          imageOutput("dataset_upload_check",inline = T)
                          ), column(9,
                                    # show the summary plots post-normalization
                                    navbarPage(inverse=F,h3("explore"),
                                               tabPanel("m/z values",# icon=icon("braille"),
                                                        fluidRow(column(6,plotOutput("var1",height='300px')),
                                                                 column(6,plotOutput("var3", height='300px'))
                                                        ),
                                                        fluidRow(column(6,plotOutput("var2", height='500px')),
                                                                 column(6,plotOutput("var4", height='500px')))
                                               ),
                                               tabPanel("samples",# icon=icon("tint"),
                                                        fluidRow(column(6,plotOutput("samp1",height='300px')),
                                                                 column(6,plotOutput("samp3", height='300px'))
                                                        ),
                                                        fluidRow(column(6,plotOutput("samp2", height='500px')),
                                                                 column(6,plotOutput("samp4", height='500px')))
                                               )
                                    )
                          )
                          )),
                 # this tab is the main analysis tab. all tabs for all analyses are listed here, but the visibility is changed depending on the current experiment
                 tabPanel("analyse",  icon = icon("bar-chart",class = "outlined"), value = "analysis",
                          sidebarLayout(position="right",
                                        mainPanel = mainPanel(width = 8,
                                                              tabsetPanel(id="statistics",selected = "pca",
                                                                          #navbarPage(inverse=F, "", id="statistics", selected = "pca", collapsible = T,
                                                                          # TODO: T-SNE
                                                                          # this tab shows general information, mostly a message with 'please give me some data' :-)
                                                                          tabPanel(icon("star"), value = "inf",
                                                                                   fluidRow(column(width=12, align="center",
                                                                                                   br(),br(),br(),br(),
                                                                                                   #hr(),
                                                                                                   #icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg"),
                                                                                                   br(),br(),
                                                                                                   h2("Please select a variable of interest in the sidebar!"),br(),
                                                                                                   icon("exchange", "fa-4x"),
                                                                                                   br(),br(),br()
                                                                                                   #hr()
                                                                                                   #icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg"), icon("arrow-right","fa-lg")
                                                                                   ))),
                                                                          tabPanel("dimension reduction", value = "dimred",  icon=icon("cube"),
                                                                                   navbarPage(inverse=T, icon("cube"), id = "dimred",
                                                                                              # loading this tab performs PCA. summary and loading tables, alongside a 2d/3d PCA plot, are available here.
                                                                                              tabPanel("pca", value = "pca", #icon=icon("cube"),
                                                                                                       fluidRow(align="center",column(12,plotly::plotlyOutput("plot_pca",height = "600px", width="600px"))),#%>% shinycssloaders::withSpinner())),
                                                                                                       fluidRow(align="center",column(12,
                                                                                                                                      switchButton("pca_2d3d", label = "", col = "BW", type = "2d3d", value=T))),
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
                                                                                                                                            plotOutput("pca_scree", width = "100%", height="250px")
                                                                                                                                   ),
                                                                                                                                   tabPanel(title="Loadings",
                                                                                                                                            div(DT::dataTableOutput('pca_load_tab',width="100%"),style='font-size:80%'))
                                                                                                                       ))
                                                                                                       )
                                                                                              ),
                                                                                              # TODO: enable the sparse and orthogonal PLS-DA options in metaboanalystR
                                                                                              # this tab is used to perform pls-da. it triggers on 'go' button as it is a time costly analysis.
                                                                                              tabPanel("pls-da", value = "plsda",
                                                                                                       fluidRow(align="center",column(12,plotly::plotlyOutput("plot_plsda",height = "500px", width="500px"))),
                                                                                                       fluidRow(align="center",column(12,
                                                                                                                                      switchButton("plsda_2d3d", label = "", col = "BW", type = "2d3d"))),
                                                                                                       hr(),
                                                                                                       fluidRow(column(3,
                                                                                                                       div(style="display:inline-block",
                                                                                                                           selectInput("plsda_type",
                                                                                                                                       label="Type:",
                                                                                                                                       choices=list("Normal"="normal")
                                                                                                                                       #,
                                                                                                                                       #             "Orthogonal"="ortho",
                                                                                                                                       #             "Sparse"="sparse")
                                                                                                                                       ,width = '100px',
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
                                                                                              tabPanel("t-sne", value = "tsne",
                                                                                                       helpText("working on it")
                                                                                              )
                                                                                   )),
                                                                          tabPanel("per m/z", value = "permz", icon=icon("fingerprint"),
                                                                                   navbarPage(inverse=T, icon("fingerprint"), id = "permz",
                                                                                              tabPanel("t-test", value="tt",
                                                                                                       fluidRow(plotly::plotlyOutput('tt_specific_plot',width="100%")),
                                                                                                       fluidRow(align="center",
                                                                                                                sardine(switchButton("tt_nonpar", "Non-parametric?", col="BW", type="YN", value = T)),
                                                                                                                #sardine(uiOutput("tt_parbutton")),
                                                                                                                sardine(switchButton("tt_eqvar", "Equal variance?", col="BW", type="YN", value = T))
                                                                                                       ),
                                                                                                       navbarPage(inverse=F,"",
                                                                                                                  tabPanel("", icon=icon("table"),
                                                                                                                           div(DT::dataTableOutput('tt_tab',width="100%"),style='font-size:80%'))
                                                                                                                  ,tabPanel("", icon=icon("area-chart"),
                                                                                                                            plotly::plotlyOutput('tt_overview_plot',height="300px") %>% shinycssloaders::withSpinner()
                                                                                                                  )
                                                                                                       )),
                                                                                              tabPanel("anova", value="aov",
                                                                                                       fluidRow(plotly::plotlyOutput('aov_specific_plot',width="100%")),
                                                                                                       navbarPage(inverse=F,"",
                                                                                                                  tabPanel("", icon=icon("table"),
                                                                                                                           div(DT::dataTableOutput('aov_tab',width="100%"),style='font-size:80%'))
                                                                                                                  ,tabPanel("", icon=icon("area-chart"),
                                                                                                                            plotly::plotlyOutput('aov_overview_plot',height="300px") %>% shinycssloaders::withSpinner()
                                                                                                                  )
                                                                                                       )),
                                                                                              tabPanel("fold-change", value="fc",
                                                                                                       fluidRow(plotly::plotlyOutput('fc_specific_plot',width="100%")),
                                                                                                       navbarPage(inverse=F,"",
                                                                                                                  tabPanel("", icon=icon("table"),
                                                                                                                           div(DT::dataTableOutput('fc_tab',width="100%"),style='font-size:80%'))
                                                                                                                  ,tabPanel("", icon=icon("area-chart"),
                                                                                                                            plotly::plotlyOutput('fc_overview_plot',height="300px") %>% shinycssloaders::withSpinner()
                                                                                                                  ))
                                                                                              ),
                                                                                              tabPanel("meba", value="meba",
                                                                                                       fluidRow(plotly::plotlyOutput('meba_specific_plot',height="600px")),
                                                                                                       fluidRow(div(DT::dataTableOutput('meba_tab', width="100%"),style='font-size:80%'))
                                                                                              ),
                                                                                              tabPanel("asca", value="asca",
                                                                                                       fluidRow(plotly::plotlyOutput('asca_specific_plot', height="600px")),
                                                                                                       fluidRow(div(DT::dataTableOutput('asca_tab',width="100%"),style='font-size:80%'))
                                                                                              )
                                                                                   )
                                                                          ),
                                                                          tabPanel("overview analyses", value = "overview", icon=icon("globe"),
                                                                                   navbarPage(inverse=T, icon("globe"), id = "overview",
                                                                                              tabPanel("volcano plot", value="volc",
                                                                                                       fluidRow(plotly::plotlyOutput('volc_plot',width="100%",height="600px") %>% shinycssloaders::withSpinner()),
                                                                                                       fluidRow(div(DT::dataTableOutput('volc_tab',width="100%"),style='font-size:80%'))
                                                                                              ),
                                                                                              tabPanel("heatmap", value="heatmap",
                                                                                                       plotly::plotlyOutput("heatmap",width="100%",height="700px") %>% shinycssloaders::withSpinner(),
                                                                                                       br(),
                                                                                                       fluidRow(column(align="center",
                                                                                                                       width=12,
                                                                                                                       sliderInput("heatmap_topn", "Use top ... from table:", value=100, min = 10, max = 200))
                                                                                                       ),
                                                                                                       fluidRow(column(align="center",
                                                                                                                       width=12,
                                                                                                                       uiOutput("heatbutton"),
                                                                                                                       switchButton("heatsign", label = "Only significant hits?", col = "GB", type = "YN"),
                                                                                                                       switchButton("heatlimits", label = "Color based on -all- metabolites?", col = "GB", type = "YN")
                                                                                                       ))
                                                                                              ),
                                                                                              # this tab is used to find overlapping features of interest between analyses
                                                                                              # TODO: enable this with multiple saved mSets in mSet$storage
                                                                                              tabPanel(title="venn", value="venn", #icon=icon("comments"),
                                                                                                       sidebarLayout(position = "left",
                                                                                                                     sidebarPanel = sidebarPanel(
                                                                                                                       fluidRow(div(DT::dataTableOutput('venn_unselected'),style='font-size:80%'), align="center"),
                                                                                                                       fluidRow(shinyWidgets::circleButton("venn_add", icon=icon("arrow-down"), size="sm"),
                                                                                                                                shinyWidgets::circleButton("venn_remove", icon=icon("arrow-up"), size="sm"),
                                                                                                                                align="center"),
                                                                                                                       fluidRow(div(DT::dataTableOutput('venn_selected'),style='font-size:80%'),align="center"),
                                                                                                                       hr(),
                                                                                                                       fluidRow(
                                                                                                                         sliderInput("venn_tophits", label = "Only include top:", min = 1, max = 200, post = " hits", value=20)
                                                                                                                         ,align="center"),
                                                                                                                       fluidRow(
                                                                                                                         shinyWidgets::circleButton("venn_build", icon=icon("hand-pointer-o"),size="default")
                                                                                                                         ,align="center")
                                                                                                                     ),
                                                                                                                     mainPanel = mainPanel(
                                                                                                                       hr(),
                                                                                                                       plotOutput("venn_plot",inline = F),
                                                                                                                       # find the overlapping compounds between the groups you want to compare (user select)
                                                                                                                       # TODO: enable this with clicking the numbers/areas
                                                                                                                       fluidRow(selectInput("intersect_venn", label = "Show hits from (only):", selected = 1,choices = "",multiple = T),
                                                                                                                                align="center"),
                                                                                                                       fluidRow(uiOutput("venn_pval"), align="center"),
                                                                                                                       br(),
                                                                                                                       fluidRow(div(DT::dataTableOutput('venn_tab'),style='font-size:80%'),
                                                                                                                                align="center")
                                                                                                                     ))
                                                                                              )
                                                                                   )
                                                                          ),
                                                                          # this tab enables machine learning
                                                                          tabPanel("machine learning", value = "ml", icon=icon("signature"),
                                                                                   br(),
                                                                                   navbarPage(inverse=F, icon("signature"), id = "ml",
                                                                                              tabPanel("initialize", value="init",
                                                                                                       fluidRow(
                                                                                                         column(width=3,align="center",
                                                                                                                selectInput("ml_perf_metr", label=h2("Performance metric"),
                                                                                                                            choices = c("boot", "boot632", "optimism_boot",
                                                                                                                                        "boot_all", "cv", "repeatedcv",
                                                                                                                                        "LOOCV", "LGOCV", "none", "oob",
                                                                                                                                        "timeslice", "addaptive_cv", "adaptive_boot",
                                                                                                                                        "adaptive_LGOCV"),
                                                                                                                            multiple = F, selected = "repeatedcv"),
                                                                                                                sliderInput("ml_train_perc",
                                                                                                                            label = h2("Percentage in training"),
                                                                                                                            min = 1,
                                                                                                                            max = 100,
                                                                                                                            step = 1,
                                                                                                                            value = 60,
                                                                                                                            post = "%"),
                                                                                                                selectInput("ml_folds", label=h2("Fold CV"),choices = c("5",
                                                                                                                                                                        "10",
                                                                                                                                                                        "20",
                                                                                                                                                                        "50",
                                                                                                                                                                        "LOOCV"),
                                                                                                                            multiple = F),
                                                                                                                sliderInput("ml_attempts",
                                                                                                                            label = "Attempts",
                                                                                                                            min = 1,
                                                                                                                            max = 100,
                                                                                                                            step = 1,
                                                                                                                            value = 20,
                                                                                                                            post = "x")
                                                                                                         ),
                                                                                                         column(width=6,align="center",
                                                                                                                selectInput("ml_method",
                                                                                                                            label = h2("Used algorithm"),
                                                                                                                            selected = "glmnet",
                                                                                                                            choices = {
                                                                                                                              lst = as.list(gbl$constants$ml.models)
                                                                                                                              # names(lst) <- sapply(gbl$constants$ml.models, function(mdl) caret.mdls[[mdl]]$label)
                                                                                                                              lst
                                                                                                                            },
                                                                                                                            multiple = F),
                                                                                                                div(uiOutput("ml_params"), style = "font-size:60%"),
                                                                                                                selectizeInput("ml_preproc", label = h2("Data reprocessing"),
                                                                                                                               choices = c("center", "scale"),
                                                                                                                               selected = c("center", "scale"), multiple=T),
                                                                                                                shinyWidgets::circleButton("do_ml",
                                                                                                                                           icon = h3(paste("Go"),
                                                                                                                                                     icon("hand-pointer-o", "fa-lg")),
                                                                                                                                           status = "default",
                                                                                                                                           size = "lg")
                                                                                                         ),
                                                                                                         column(width=3,align="center",
                                                                                                                fluidRow(textOutput("ml_train_ss"),
                                                                                                                         actionButton("ml_train_ss", label = "train on:", icon = icon("arrow-up"))),
                                                                                                                fluidRow(textOutput("ml_test_ss"),
                                                                                                                         actionButton("ml_test_ss", label = "test on:", icon = icon("arrow-up"))),
                                                                                                                br(),
                                                                                                                textInput("ml_name", label=h3("Name:"), value = "all"))
                                                                                                       )
                                                                                              ),
                                                                                              tabPanel("results", value="res", icon=icon("poll"),
                                                                                                       br(),
                                                                                                       div(selectInput("show_which_ml", label = "Plot which model?", choices = c())),
                                                                                                       navbarPage(title=icon("poll"),id="ml_results",inverse=F,
                                                                                                                  tabPanel(title = "roc",value = "roc",icon=icon("area-chart"),
                                                                                                                           plotlyOutput("ml_roc",height = "600px"),
                                                                                                                           div(DT::dataTableOutput("ml_tab",width="100%"),style='font-size:80%')),
                                                                                                                  tabPanel("importance",value= "bar",icon=icon("star"),
                                                                                                                           fluidRow(plotlyOutput("ml_bar", width = "100%", height="600px")),
                                                                                                                           fluidRow(
                                                                                                                             column(12, sliderInput("ml_top_x",
                                                                                                                                                    label = "Show top:",
                                                                                                                                                    min = 10,
                                                                                                                                                    max = 200,
                                                                                                                                                    step=10,
                                                                                                                                                    value=20), align="center")
                                                                                                                           )
                                                                                                                  )
                                                                                                       )
                                                                                              )
                                                                                   )
                                                                          )
                                                              )
                                        ),
                                        # this is the sidebar that shows in the analysis tab. contains a lot of settings on the current variable of interest, plot themes and colours, and venn diagrams.
                                        sidebarPanel =
                                          sidebarPanel(align="center",width = 4,
                                                       tabsetPanel(id = "anal_sidebar", selected="switch/subset",#type = "pills",
                                                                   tabPanel(title="search", icon=icon("search"),
                                                                            br(),
                                                                            bsCollapse(bsCollapsePanel(title=h2("Settings"), style="info",
                                                                                                       tabsetPanel(id="tab_iden_1", selected = "start",
                                                                                                                   # forward searching
                                                                                                                   tabPanel(title=icon("database"),value="start",
                                                                                                                            uiOutput("db_search_select"),
                                                                                                                            div(id = "curly-brace", div(id = "left", class = "brace"),
                                                                                                                                div(id = "right", class = "brace")),
                                                                                                                            br(),br(),
                                                                                                                            shinyWidgets::circleButton("select_db_all",
                                                                                                                                                       icon = icon("shopping-cart"),
                                                                                                                                                       size = "default") # icon("fingerprint"), size = "sm")
                                                                                                                   ), # clicky buttons for database selection; this is generated in 'server'
                                                                                                                   tabPanel(title=icon("chart-bar"),
                                                                                                                            plotly::plotlyOutput("curr_plot", height="300px", width="100%") %>% shinycssloaders::withSpinner()
                                                                                                                   ),
                                                                                                                   tabPanel(title=icon("magic"),
                                                                                                                            h2("MagicBall settings"),
                                                                                                                            fluidRow(align="center",switchButton(inputId = "magicball_pubchem_cids",
                                                                                                                                                                 label = "Check PubChem for predicted formulas?",
                                                                                                                                                                 col = "BW", type = "YN", value = F)
                                                                                                                            ),
                                                                                                                            fluidRow(align="center",switchButton(inputId = "magicball_pubchem_details",
                                                                                                                                                                 label = "Get detailed PubChem matches? (SLOW!)",
                                                                                                                                                                 col = "BW", type = "YN", value = F)
                                                                                                                            ),
                                                                                                                            fluidRow(align="center", helpText("Considered adducts:")),
                                                                                                                            fluidRow(div(DT::dataTableOutput('magicball_add_tab'),style='font-size:100%'),
                                                                                                                                     align="center")
                                                                                                                   ),
                                                                                                                   tabPanel(title=icon("star-half-alt"),
                                                                                                                            selectInput("iso_score_method",
                                                                                                                                        "Which method used to score compounds of same weight?",
                                                                                                                                        selected="mscore",
                                                                                                                                        choices=list("M-score"="mscore"
                                                                                                                                                     #"Chi-square"="chisq",
                                                                                                                                                     #"Mean absolute percentage error"="mape",
                                                                                                                                                     #"SIRIUS"="sirius",
                                                                                                                                                     #"Network-based"="network"
                                                                                                                                        )),
                                                                                                                            sliderInput("int_prec", label = "Intensity imprecision", min = 1, max = 100, value = 2, post = "%"),
                                                                                                                            shinyWidgets::circleButton("score_iso", icon = icon("award"), size = "sm") # icon("fingerprint"), size = "sm")
                                                                                                                   )
                                                                                                       ))),
                                                                            tabsetPanel(id="tab_iden_2",
                                                                                        # forward searching
                                                                                        tabPanel(title="mz > molecule",
                                                                                                 hr(),
                                                                                                 fluidRow(
                                                                                                   tags$button(
                                                                                                     id = "search_mz",
                                                                                                     class = "btn btn-default action-button",
                                                                                                     img(src = "detective.png",
                                                                                                         height = "50px")
                                                                                                   ),
                                                                                                   div(
                                                                                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                                                                                 style="position:relative;
                                                                                                         top:10px;")),
                                                                                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                                                                                 style="position:relative;
                                                                                                         top:25px;")),
                                                                                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                                                                                 style="position:relative;
                                                                                                         top:10px;")),
                                                                                                     sardine(h2(textOutput("curr_mz"),style="padding:10px;")),
                                                                                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                                                                                 style="position:relative;
                                                                                                         top:10px;")),
                                                                                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                                                                                 style="position:relative;
                                                                                                         top:25px;")),
                                                                                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                                                                                 style="position:relative;
                                                                                                         top:10px;")),
                                                                                                     style="background-color:white;
                                                                                             height:55px;
                                                                                             width:115%;
                                                                                             position:relative;
                                                                                             right:30px;
                                                                                             border-top: 1px solid #DFDCDC;
                                                                                             border-bottom: 1px solid #DFDCDC;")
                                                                                                 ),
                                                                                                 bsCollapse(bsCollapsePanel(title=h2("Compound info"), style="warning",
                                                                                                                            tabsetPanel(id="tab_iden_3",
                                                                                                                                        tabPanel(title=icon("atlas"),
                                                                                                                                                 wellPanel(id = "def",style = "overflow-y:scroll; max-height: 200px",
                                                                                                                                                           textOutput("curr_definition"))
                                                                                                                                        ),
                                                                                                                                        tabPanel(title=icon("atom"),
                                                                                                                                                 textOutput("curr_formula"),
                                                                                                                                                 plotOutput("curr_struct", height="300px")
                                                                                                                                        )
                                                                                                                            ))),
                                                                                                 bsCollapse(bsCollapsePanel(title=h2("Search results"), style="error",
                                                                                                                            tabsetPanel(id="tab_iden_4",selected = "start",
                                                                                                                                        # forward searching
                                                                                                                                        tabPanel(title=icon("table"), value="start",
                                                                                                                                                 div(DT::dataTableOutput('match_tab', width="100%"),style='font-size:80%'),
                                                                                                                                                 hr(),
                                                                                                                                                 fluidRow(
                                                                                                                                                   switchButton(inputId = "auto_copy",
                                                                                                                                                                label = "Auto-copy name to clipboard??",
                                                                                                                                                                value = TRUE, col = "GB", type = "YN"),
                                                                                                                                                   align="center"),
                                                                                                                                                 helpText("Undo filtering"),
                                                                                                                                                 fluidRow(
                                                                                                                                                   shinyWidgets::circleButton("undo_match_filt", icon = icon("undo-alt"), size = "sm") # icon("fingerprint"), size = "sm")
                                                                                                                                                 )
                                                                                                                                        ),
                                                                                                                                        tabPanel(title=icon("database"), value="pie_db",
                                                                                                                                                 fluidRow(align = "center",
                                                                                                                                                          plotly::plotlyOutput("match_pie_db") %>% shinycssloaders::withSpinner()
                                                                                                                                                 )
                                                                                                                                        ),
                                                                                                                                        tabPanel(title=icon("plus"), value = "pie_add",
                                                                                                                                                 fluidRow(align = "center",
                                                                                                                                                          plotly::plotlyOutput("match_pie_add") %>% shinycssloaders::withSpinner()
                                                                                                                                                 )
                                                                                                                                        ),
                                                                                                                                        tabPanel(title=icon("cloud"), value = "word_cloud",
                                                                                                                                                 fluidRow(align = "center",
                                                                                                                                                          conditionalPanel("input.wc_cloudbar == true",
                                                                                                                                                                           tagList(
                                                                                                                                                                             wordcloud2::wordcloud2Output("wordcloud_desc"),
                                                                                                                                                                             tags$script(HTML(
                                                                                                                                                                               "$(document).on('click', '#canvas', function() {",
                                                                                                                                                                               'word = document.getElementById("wcSpan").innerHTML;',
                                                                                                                                                                               "Shiny.onInputChange('selected_word_desc', word);",
                                                                                                                                                                               "});"
                                                                                                                                                                             ))
                                                                                                                                                                           )),
                                                                                                                                                          conditionalPanel("input.wc_cloudbar == false",
                                                                                                                                                                           plotlyOutput("wordbar_desc")
                                                                                                                                                          ),
                                                                                                                                                          sliderInput("wc_topn", "Top words shown:", min = 1, max = 100,step = 1,value = 30, width="60%"),
                                                                                                                                                          switchButton("wc_cloudbar", label = "", col = "BW", type = "CLBR",value = T)
                                                                                                                                                 )
                                                                                                                                        ),
                                                                                                                                        tabPanel(title=icon("searchengin"),
                                                                                                                                                 textInput('pm_query', "Search for:"),
                                                                                                                                                 sliderInput('pm_year', "Paper publication range:",
                                                                                                                                                             min = 1900, max = as.numeric(format(Sys.Date(), '%Y')),
                                                                                                                                                             value = c(2000,as.numeric(format(Sys.Date(), '%Y'))),
                                                                                                                                                             step = 1,sep = ""
                                                                                                                                                 ),
                                                                                                                                                 sliderInput("pm_max",
                                                                                                                                                             "Stop after ... papers:",
                                                                                                                                                             min = 1,
                                                                                                                                                             max = 1000,
                                                                                                                                                             value = 500),
                                                                                                                                                 shinyWidgets::circleButton("search_pubmed", icon = icon("search"), size = "sm"),
                                                                                                                                                 tabsetPanel(selected = 1,
                                                                                                                                                             tabPanel(title = icon("cloud"),
                                                                                                                                                                      conditionalPanel("input.wc_cloudbar_pm == true",
                                                                                                                                                                                       tagList(
                                                                                                                                                                                         wordcloud2::wordcloud2Output("wordcloud_desc_pm"),
                                                                                                                                                                                         tags$script(HTML(
                                                                                                                                                                                           "$(document).on('click', '#canvas', function() {",
                                                                                                                                                                                           'word = document.getElementById("wcSpan").innerHTML;',
                                                                                                                                                                                           "Shiny.onInputChange('selected_word_desc_pm', word);",
                                                                                                                                                                                           "});"
                                                                                                                                                                                         ))
                                                                                                                                                                                       )),
                                                                                                                                                                      conditionalPanel("input.wc_cloudbar_pm == false",
                                                                                                                                                                                       plotlyOutput("wordbar_desc_pm")
                                                                                                                                                                      ),
                                                                                                                                                                      sliderInput("wc_topn_pm", "Top words shown:", min = 1, 
                                                                                                                                                                                  max = 100, step = 1,value = 30, width = "60%"),
                                                                                                                                                                      switchButton("wc_cloudbar_pm", label = "", col = "BW", type = "CLBR",value = F)
                                                                                                                                                             ),
                                                                                                                                                             tabPanel(title = icon("table"),
                                                                                                                                                                      div(DT::dataTableOutput('pm_tab', width="100%"),
                                                                                                                                                                          style='font-size:80%')
                                                                                                                                                             )
                                                                                                                                                 )
                                                                                                                                                 
                                                                                                                                        )
                                                                                                                            )))
                                                                                                 
                                                                                        ),
                                                                                        # reverse searching
                                                                                        tabPanel(title="molecule > mz",
                                                                                                 br(),
                                                                                                 actionButton("browse_db", "Browse compounds", icon=icon("eye")),
                                                                                                 hr(),
                                                                                                 tabsetPanel(
                                                                                                   tabPanel(NULL, icon = icon("database"),
                                                                                                            wellPanel(id = "def",style = "overflow-y:scroll; max-height: 200px",
                                                                                                                      textOutput("browse_definition")),
                                                                                                            div(DT::dataTableOutput('browse_tab'),style='font-size:80%'),
                                                                                                            hr(),
                                                                                                            actionButton("revsearch_mz", "Find hits", icon=icon("search"))
                                                                                                   ),
                                                                                                   tabPanel(NULL, icon = icon("search-location"),
                                                                                                            div(DT::dataTableOutput('hits_tab'),style='font-size:80%')
                                                                                                   )
                                                                                                 )
                                                                                                 
                                                                                        ))
                                                                   ),
                                                                   tabPanel("switch/subset", icon=icon("exchange")
                                                                            ,h2("Current experiment:")
                                                                            ,div(
                                                                              sardine(h2(textOutput("curr_name"),style="padding:10px;")),
                                                                              style="background-color:white;
                                                                      height:55px;
                                                                      width:115%;
                                                                      position:relative;
                                                                      right:30px;
                                                                      border-top: 1px solid #DFDCDC;
                                                                      border-bottom: 1px solid #DFDCDC;")
                                                                            ,hr()
                                                                            ,h2("Change variable of interest")
                                                                            ,selectInput("stats_var", label="Do statistics on:", choices = c("label"))
                                                                            ,shinyWidgets::circleButton("change_cls", icon = icon("hand-pointer-o"), size = "sm")
                                                                            ,fluidRow(column(12, align="center", uiOutput("timebutton")))
                                                                            ,hr()
                                                                            ,h2("Subset data")
                                                                            ,selectInput("subset_var", label="Subset data based on:", choices = c("label"))
                                                                            ,selectizeInput("subset_group", label="Group(s) in subset:", choices = c(), multiple=TRUE)
                                                                            ,shinyWidgets::circleButton("change_subset", icon = icon("hand-pointer-o"), size = "sm")
                                                                            ,shinyWidgets::circleButton("reset_subset", icon = icon("undo"), size = "sm")
                                                                            
                                                                   ),
                                                                   # this tab is used to select user plot theme and user colours (discrete and continuous)
                                                                   tabPanel("plot aesthetics", icon=icon("paint-brush"),
                                                                            h2("Summary plot style"),br(),
                                                                            selectizeInput("ggplot_sum_style", multiple=T, label = "Style(s)", choices = list("Box"="box",
                                                                                                                                                              "Violin"="violin",
                                                                                                                                                              "Beeswarm"="beeswarm",
                                                                                                                                                              "Scatterplot"="scatter"),
                                                                                           selected = c("violin")
                                                                            ),
                                                                            selectInput("ggplot_sum_stats", label = "Stats shown", choices = list("median", "mean", "none")),
                                                                            h2("Shape")
                                                                            ,selectInput("shape_var", label="Marker shape based on:", choices = c("label"))
                                                                            ,h2("Color")
                                                                            ,selectInput("col_var", label="Marker color based on:", choices = c("label"))
                                                                            ,h2("Hover text")
                                                                            ,selectInput("txt_var", label="Marker hover text based on:", choices = c("label")),
                                                                            h2("Plot theme"),
                                                                            selectInput("ggplot_theme", label = "Theme", choices = list("Grid, white bg"="bw",
                                                                                                                                        "No grid, white bg"="classic",
                                                                                                                                        "Grid, gray bg"="gray",
                                                                                                                                        "Minimal"="min",
                                                                                                                                        "Grid, black bg"="dark",
                                                                                                                                        "Grid, white bg, gray axes"="light",
                                                                                                                                        "Line drawing"="line"),
                                                                                        selected = opts$gtheme),
                                                                            fluidRow(plotOutput("ggplot_theme_example",inline = F, width="100%")),
                                                                            h2("Continuous data"),
                                                                            # the below options need to match with the corresponding function storage in 'global'. if you want to add more it'll go here!
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
                                                                                                                                             "Pastel2", "Set1", "Set2", "Set3"),selected = opts$gspec
                                                                            ),
                                                                            # preview plot
                                                                            fluidRow(plotly::plotlyOutput("ramp_plot",inline = T, width="100%") %>% shinycssloaders::withSpinner()),
                                                                            h2("Discrete data"),
                                                                            uiOutput("colorPickers") # colour pickers generated in server.R. default settings taken from user_options.txt.
                                                                   ))
                                          )
                          )),
                 
                 # report tab
                 tabPanel("report",
                          icon = icon("file-invoice", class = "outlined"),
                          value="reportTab",
                          fluidRow(
                            column(width=12, align="center",
                                   h2("Report"),
                                   br(),
                                   helpText("Report contents:"),
                                   div(DT::dataTableOutput('report_unselected',  width="100%"))
                            )#close column
                          )#close fluidrow
                 ),#close tabpanel
                 
                 # this tab is used to change general settings.
                 tabPanel("settings",  icon = icon("cog",class = "outlined"), value="options",
                          navbarPage(inverse=TRUE,"Settings", id="tab_settings",
                                     tabPanel("Mode", icon=icon("box-open"),
                                              switchButton(inputId = "db_only", label = "Run in database-only mode?", 
                                                           value = switch(opts$mode, dbonly=T, complete=F), 
                                                           col = "BW", type = "YN")
                                     ),
                                     tabPanel("Project", icon=icon("gift"),
                                              #textInput(inputId="proj_name", label="Project name", value = ''),
                                              selectizeInput(inputId="proj_name",
                                                             label="Project name",
                                                             choices=lcl$vectors$project_names, # existing projects in user folder (generated in 'global')
                                                             selected = opts$proj_name,
                                                             options=list(create = TRUE)), # let users add new names
                                              actionButton("set_proj_name", label="Apply"),
                                              helpText("This name will be used in all save files."),
                                              textOutput("proj_name")
                                     ),
                                     # change list of adducts used, or add your own
                                     # TODO: fix, i think this is currently non-functional
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
                                     # change toolbar colour, text font and size
                                     tabPanel("Aesthetic", icon=icon("child"),
                                              h3("Change app settings"),
                                              hr(),
                                              h2("Navigation bar colours"),
                                              colourpicker::colourInput(inputId = "bar.col.1",
                                                                        label = paste("Active background"),
                                                                        value = opts$col1,
                                                                        allowTransparent = FALSE),
                                              colourpicker::colourInput(inputId = "bar.col.2",
                                                                        label = paste("Inactive background"),
                                                                        value = opts$col2,
                                                                        allowTransparent = FALSE),
                                              colourpicker::colourInput(inputId = "bar.col.3",
                                                                        label = paste("Active tab"),
                                                                        value = opts$col3,
                                                                        allowTransparent = FALSE),
                                              colourpicker::colourInput(inputId = "bar.col.4",
                                                                        label = paste("Inactive tab"),
                                                                        value = opts$col4,
                                                                        allowTransparent = FALSE),
                                              br(),
                                              h2("Fonts (Google fonts)"),
                                              textInput(inputId="font.1", label="h1", value = opts$font1),
                                              textInput(inputId="font.2", label="h2", value = opts$font2),
                                              textInput(inputId="font.3", label="h3", value = opts$font3),
                                              textInput(inputId="font.4", label="body", value = opts$font4),
                                              br(), # TODO: font size modifier slider
                                              h2("Font size"),
                                              sliderInput("size.1", label="h1", value=as.numeric(opts$size1),min = 5, max=50),
                                              sliderInput("size.2", label="h2", value=as.numeric(opts$size2),min = 5, max=50),
                                              sliderInput("size.3", label="h3", value=as.numeric(opts$size3),min = 5, max=50),
                                              sliderInput("size.4", label="body", value=as.numeric(opts$size4),min = 5, max=50),
                                              br(),
                                              h3("Taskbar image"),
                                              div(imageOutput("taskbar_image",inline = T)),
                                              shinyFilesButton('taskbar_image_path',
                                                               'Select image',
                                                               'Please select an image file',
                                                               FALSE),
                                              hr(),
                                              actionButton("change_css", "Save settings (restart to apply)") # need to reload CSS to enable new settings
                                     )
                          )
                 ),
                 # prompt user on opening the quit tab.
                 # TODO: add 'save project?' dialog
                 #tabPanel(title = "", value="stop", icon = icon("times-circle",class = "outlined")),
                 div(class="spinnylocation1",
                     div(class="plus", img(class="imagetop", src=opts$taskbar_image, width="100px", height="100px")),
                     div(class="minus", img(class="imagebottom", src=opts$taskbar_image, width="100px", height="100px"))
                 ),
                 div(class="line")
                 ,footer=fluidRow(
                   br(),br(),br(),
                   div(
                     #actionButton("show_window", label="", icon = icon("map-marked")),
                     actionButton("load_mset", label="load", icon = icon("folder-open"),style=gsubfn::fn$paste("background-color:$bgcol; border-color:$bgcol;")),
                     actionButton("save_mset", label="save", icon = icon("save"),style=gsubfn::fn$paste("background-color:$bgcol; border-color:$bgcol;")),
                     actionButton("debug", label="debug", icon = icon("bug"),style=gsubfn::fn$paste("background-color:$bgcol; border-color:$bgcol;"))
                     , style=gsubfn::fn$paste("position:fixed;bottom:0;width:100%;height:40px;z-index:1005;background-color:$bgcol;border-style:solid; border-color:black;border-width:1px;")),
                   align="center")
      )
    )
  }else if (logged$status == "setfolder"){
    # prompt file location
    
    fluidRow(align="center",
             br(),br(),br(),br(),br(),
             imageOutput("login_header",inline = T),
             helpText("Welcome to MetaboShiny! Please pick a folder to save your files and databases in with the below button."),
             shinyDirButton("get_work_dir", "Choose a folder",
                            title = "Browse",
                            buttonType = "default", class = NULL),
             verbatimTextOutput("curr_exp_dir_start"),
             hr(),
             shinyWidgets::circleButton("confirm_work_dir", icon = icon("check"))
    )
    
  }else if(logged$status == "db_only"){
    print("db only mode")
    
    print("rendering...")
    # read settings
    
    opts <- getOptions(lcl$paths$opt.loc)
    
    logged$text <<- "loaded options!"
    logged$text <<- "starting MetaboShiny..."
    
    online = internetWorks()
    
    # generate CSS for the interface based on user settings for colours, fonts etc.
    bar.css <<- nav.bar.css(opts$col1, opts$col2, opts$col3, opts$col4)
    font.css <<- app.font.css(opts$font1, opts$font2, opts$font3, opts$font4,
                              opts$size1, opts$size2, opts$size3, opts$size4, online=online)
    
    # google fonts
    
    # name = the name of the font in Google's Library (https://fonts.google.com)
    # family = how you want to refer to the font inside R
    # regular.wt = the weight of font used in axis, labels etc.
    # bolt.wt = the weight of font used in the title
    
    # === GOOGLE FONT SUPPORT FOR GGPLOT2 ===
    
    # Download a webfont
    if(online){
      lapply(c(opts[grepl(pattern = "font", names(opts))]), function(font){
        try({
          sysfonts::font_add_google(name = font,
                                    family = font,
                                    regular.wt = 400,
                                    bold.wt = 700)
        })
      })
    }
    
    # Perhaps the only tricky bit is remembering to run the following function to enable webfonts
    showtext::showtext_auto(enable = T)
    
    # ======================================
    
    # set taskbar image as set in options
    taskbar_image <- opts$task_img
    
    tagList(
      tags$title('MetaboShiny'),
      tags$style(type="text/css", bar.css),
      navbarPage(inverse=TRUE,#tags$head(tags$script(src="sparkle.js")),
                 title=div(h1("MetaboShiny"), class="outlined", tags$style(type="text/css", font.css)), # make it use the sparkle.js for unnecessary sparkle effects ;)
                 id="nav_general",
                 # this tab shows the available databases, if they are installed, and buttons to install them. generated as output$db_build_ui in 'server'
                 tabPanel("database", icon = icon("database",class = "outlined"), value="database",
                          uiOutput("db_build_ui") #%>% shinycssloaders::withSpinner() # see server, is autogenerated now
                 ),
                 tabPanel("search", icon = icon("search",class = "outlined"), value="search",
                          fluidRow(align="center",
                                   div(
                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                 style="position:relative;
                                                                                                         top:10px;")),
                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                 style="position:relative;
                                                                                                         top:25px;")),
                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                 style="position:relative;
                                                                                                         top:10px;")),
                                     sardine(div(textInput("search_mz",value = "m/z here",label = NULL),style="font-size:40px;")),
                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                 style="position:relative;
                                                                                                         top:10px;")),
                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                 style="position:relative;
                                                                                                         top:25px;")),
                                     sardine(div(icon("paw","fa-xs fa-rotate-90"),
                                                 style="position:relative;
                                                                                                         top:10px;"))),
                                   shinyWidgets::circleButton("search_mz", icon = icon("search"), size = "lg"),
                                   sidebarLayout(position="left",
                                                 sidebarPanel = sidebarPanel(align="center",
                                                                             tabsetPanel(
                                                                               tabPanel(title = "settings", icon = icon("cog"),
                                                                                        shinyWidgets::knobInput(
                                                                                          inputId = "ppm",
                                                                                          label = "ppm accuracy:",
                                                                                          value = 5,
                                                                                          min = 1,
                                                                                          displayPrevious = TRUE, 
                                                                                          lineCap = "round",
                                                                                          fgColor = opts$col1,
                                                                                          inputColor = opts$col1
                                                                                        )),
                                                                               tabPanel(title = "databases", icon=icon("help"),
                                                                                        uiOutput("db_search_select"),
                                                                                        div(id = "curly-brace", 
                                                                                            div(id = "left", class = "brace"),
                                                                                            div(id = "right", class = "brace")),
                                                                                        shinyWidgets::circleButton("select_db_all",
                                                                                                                   icon = icon("shopping-cart"),
                                                                                                                   size = "default"))                
                                                 )),
                                                 # this tab shows a summary of the created/loaded in csv file.
                                                 mainPanel = mainPanel(align="center",
                                                                       tabsetPanel(id="tab_iden_4",selected = "start",
                                                                                   # forward searching
                                                                                   tabPanel(title=icon("table"), value="start",
                                                                                            fluidRow(column(6, wellPanel(id = "def",style = "overflow-y:scroll; max-height: 200px",
                                                                                                                         textOutput("curr_definition"))),
                                                                                                     column(6, textOutput("curr_formula"),
                                                                                                            plotOutput("curr_struct", height="300px"))),
                                                                                            hr(),
                                                                                            div(DT::dataTableOutput('match_tab', width="100%"),style='font-size:80%'),
                                                                                            helpText("Undo filtering"),
                                                                                            fluidRow(
                                                                                              shinyWidgets::circleButton("undo_match_filt", icon = icon("undo-alt"), size = "sm") # icon("fingerprint"), size = "sm")
                                                                                            )
                                                                                   ),
                                                                                   tabPanel(title=icon("plus"), value = "pie_add",
                                                                                            fluidRow(align = "center",
                                                                                                     plotly::plotlyOutput("match_pie_add"),
                                                                                                     plotly::plotlyOutput("match_pie_db"),
                                                                                                     conditionalPanel("input.wc_cloudbar == true",
                                                                                                                      tagList(
                                                                                                                        wordcloud2::wordcloud2Output("wordcloud_desc"),
                                                                                                                        tags$script(HTML(
                                                                                                                          "$(document).on('click', '#canvas', function() {",
                                                                                                                          'word = document.getElementById("wcSpan").innerHTML;',
                                                                                                                          "Shiny.onInputChange('selected_word_desc', word);",
                                                                                                                          "});"
                                                                                                                        ))
                                                                                                                      )),
                                                                                                     conditionalPanel("input.wc_cloudbar == false",
                                                                                                                      plotlyOutput("wordbar_desc")
                                                                                                     ),
                                                                                                     sliderInput("wc_topn", "Top words shown:", min = 1, max = 100,step = 1,value = 30, width="60%"),
                                                                                                     switchButton("wc_cloudbar", label = "", col = "BW", type = "CLBR",value = T)
                                                                                            )),
                                                                                   tabPanel(title=icon("searchengin"),
                                                                                            textInput('pm_query', "Search for:"),
                                                                                            sliderInput('pm_year', "Paper publication range:",
                                                                                                        min = 1900, max = as.numeric(format(Sys.Date(), '%Y')),
                                                                                                        value = c(2000,as.numeric(format(Sys.Date(), '%Y'))),
                                                                                                        step = 1,sep = ""
                                                                                            ),
                                                                                            sliderInput("pm_max",
                                                                                                        "Stop after ... papers:",
                                                                                                        min = 1,
                                                                                                        max = 1000,
                                                                                                        value = 500),
                                                                                            shinyWidgets::circleButton("search_pubmed", icon = icon("search"), size = "sm"),
                                                                                            tabsetPanel(selected = 1,
                                                                                                        tabPanel(title = icon("cloud"),
                                                                                                                 conditionalPanel("input.wc_cloudbar_pm == true",
                                                                                                                                  tagList(
                                                                                                                                    wordcloud2::wordcloud2Output("wordcloud_desc_pm"),
                                                                                                                                    tags$script(HTML(
                                                                                                                                      "$(document).on('click', '#canvas', function() {",
                                                                                                                                      'word = document.getElementById("wcSpan").innerHTML;',
                                                                                                                                      "Shiny.onInputChange('selected_word_desc_pm', word);",
                                                                                                                                      "});"
                                                                                                                                    ))
                                                                                                                                  )),
                                                                                                                 conditionalPanel("input.wc_cloudbar_pm == false",
                                                                                                                                  plotlyOutput("wordbar_desc_pm")
                                                                                                                 ),
                                                                                                                 sliderInput("wc_topn_pm", "Top words shown:", min = 1, 
                                                                                                                             max = 100, step = 1,value = 30, width = "60%"),
                                                                                                                 switchButton("wc_cloudbar_pm", label = "", col = "BW", type = "CLBR",value = F)
                                                                                                        ),
                                                                                                        tabPanel(title = icon("table"),
                                                                                                                 div(DT::dataTableOutput('pm_tab', width="100%"),
                                                                                                                     style='font-size:80%')
                                                                                                        )
                                                                                            )
                                                                                            
                                                                                   )
                                                                       )
                                                 )
                                   )
                 )),
      tabPanel("settings",  icon = icon("cog",class = "outlined"), value="options",
               navbarPage(inverse=TRUE,"Settings", id="tab_settings",
                          tabPanel("Mode", icon=icon("box-open"),
                                   switchButton(inputId = "db_only", label = "Run in database-only mode?",
                                                value = switch(opts$mode, dbonly=T, complete=F), 
                                                col = "BW", type = "YN")
                          ),
                          tabPanel("Project", icon=icon("gift"),
                                   #textInput(inputId="proj_name", label="Project name", value = ''),
                                   selectizeInput(inputId="proj_name",
                                                  label="Project name",
                                                  choices=lcl$vectors$project_names, # existing projects in user folder (generated in 'global')
                                                  selected = opts$proj_name,
                                                  options=list(create = TRUE)), # let users add new names
                                   actionButton("set_proj_name", label="Apply"),
                                   helpText("This name will be used in all save files."),
                                   textOutput("proj_name")
                          ),
                          # user directory picking
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
                          # change list of adducts used, or add your own
                          # TODO: fix, i think this is currently non-functional
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
                          # change toolbar colour, text font and size
                          tabPanel("Aesthetic", icon=icon("child"),
                                   h3("Change app settings"),
                                   hr(),
                                   h2("Navigation bar colours"),
                                   colourpicker::colourInput(inputId = "bar.col.1",
                                                             label = paste("Active background"),
                                                             value = opts$col1,
                                                             allowTransparent = FALSE),
                                   colourpicker::colourInput(inputId = "bar.col.2",
                                                             label = paste("Inactive background"),
                                                             value = opts$col2,
                                                             allowTransparent = FALSE),
                                   colourpicker::colourInput(inputId = "bar.col.3",
                                                             label = paste("Active tab"),
                                                             value = opts$col3,
                                                             allowTransparent = FALSE),
                                   colourpicker::colourInput(inputId = "bar.col.4",
                                                             label = paste("Inactive tab"),
                                                             value = opts$col4,
                                                             allowTransparent = FALSE),
                                   br(),
                                   h2("Fonts (Google fonts)"),
                                   textInput(inputId="font.1", label="h1", value = opts$font1),
                                   textInput(inputId="font.2", label="h2", value = opts$font2),
                                   textInput(inputId="font.3", label="h3", value = opts$font3),
                                   textInput(inputId="font.4", label="body", value = opts$font4),
                                   br(), # TODO: font size modifier slider
                                   h2("Font size"),
                                   sliderInput("size.1", label="h1", value=as.numeric(opts$size1),min = 5, max=50),
                                   sliderInput("size.2", label="h2", value=as.numeric(opts$size2),min = 5, max=50),
                                   sliderInput("size.3", label="h3", value=as.numeric(opts$size3),min = 5, max=50),
                                   sliderInput("size.4", label="body", value=as.numeric(opts$size4),min = 5, max=50),
                                   br(),
                                   h3("Taskbar image"),
                                   div(imageOutput("taskbar_image",inline = T)),
                                   shinyFilesButton('taskbar_image_path',
                                                    'Select image',
                                                    'Please select an image file',
                                                    FALSE),
                                   hr(),
                                   actionButton("change_css", "Save settings (restart to apply)") # need to reload CSS to enable new settings
                          )
               )
      )))
  }
})
