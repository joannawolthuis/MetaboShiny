library(shiny)
library(shinyBS)

shinyServer(function(input, output, session) {
  
# ================================= DEFAULTS ===================================

tbl <<- NA
tables <- list()
db_search_list <- c()
color.function <- rainbow
patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
mainmode <<- "stat"
shinyOptions(progress.style="old")

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

default.text <- list(list(name='options$work_dir',text=options$work_dir),
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


# -------------------
time.anal.ui <- reactive({
  # -------------
  navbarPage("Time Series", id="nav_time",
             tabPanel("iPCA", value = "ipca", 
                      plotlyOutput("plot_ipca",width="100%"),
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
                      fluidRow(plotlyOutput('meba_specific_plot'),width="100%"),
                      fluidRow(div(DT::dataTableOutput('meba_tab', width="100%"),style='font-size:80%'))
                      ),
             # =================================================================================
             tabPanel("ASCA", value="asca",
                      fluidRow(plotlyOutput('asca_specific_plot', width="100%")),
                      fluidRow(div(DT::dataTableOutput('asca_tab',width="100%"),style='font-size:80%'))
             )
             # =================================================================================
  )
})

stat.anal.ui <- reactive({
  navbarPage("Standard analysis", id="tab_stat",
             tabPanel("", value = "intro", icon=icon("comment-o"),
                      helpText("Info text here")
             ),
             tabPanel("PCA", value = "pca", #icon=icon("cube"),
                      plotlyOutput("plot_pca",width="100%"),
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
                                          plotlyOutput("plot_plsda",width="100%"),
                                          fluidRow(column(3,
                                                          selectInput("plsda_x", label = "X axis:", choices = paste0("PC",1:8),selected = "PC1",width="100%"),
                                                          selectInput("plsda_y", label = "Y axis:", choices = paste0("PC",1:8),selected = "PC2",width="100%"),
                                                          selectInput("plsda_z", label = "Z axis:", choices = paste0("PC",1:8),selected = "PC3",width="100%")
                                          ), 
                                          column(9, 
                                                 div(DT::dataTableOutput('plsda_tab',width="100%"),style='font-size:80%')
                                          ))),
                                 tabPanel("", icon=icon("star-o"), 
                                          plotlyOutput("plsda_vip_specific_plot",width="100%"),
                                          fluidRow(column(3,
                                                          selectInput("plsda_vip_cmp", label = "Compounds from:", choices = paste0("PC",1:5),selected = "PC1",width="100%")
                                          ), 
                                          column(9, 
                                                 div(DT::dataTableOutput('plsda_vip_tab',width="100%"),style='font-size:80%')
                                          ))
                                          
                                 )
                                 )
                                          
             ),
             # =================================================================================
             tabPanel("T-test", value="tt", 
                      fluidRow(plotlyOutput('tt_specific_plot',width="100%")),
                      navbarPage("",
                                 tabPanel("", icon=icon("table"),
                                          div(DT::dataTableOutput('tt_tab',width="100%"),style='font-size:80%'))
                                 ,tabPanel("", icon=icon("area-chart"),
                                           plotlyOutput('tt_overview_plot',height="300px")
                                 )
                      )),
             tabPanel("Fold-change", value="fc",
                      fluidRow(plotlyOutput('fc_specific_plot',width="100%")),
                      navbarPage("",
                                 tabPanel("", icon=icon("table"),
                                          div(DT::dataTableOutput('fc_tab',width="100%"),style='font-size:80%'))
                                 ,tabPanel("", icon=icon("area-chart"),
                                           plotlyOutput('fc_overview_plot',height="300px")
                                 )
                      )),
             tabPanel("Heatmap", value="heatmap",
                      plotlyOutput("heatmap",width="110%",height="700px"),
                      br(),
                      fluidRow(column(align="center",
                                      width=12,switchButton(inputId = "heatmode",
                                            label = "Use data from:", 
                                            value = TRUE, col = "BW", type = "TTFC"))
                      )
             ),
             # =================================================================================
             tabPanel("Volcano", value="volc",
                      fluidRow(plotlyOutput('volc_plot',width="100%",height="600px"))
                      #,fluidRow(div(DT::dataTableOutput('volc_tab'),style='font-size:80%'))
                      )
  )
})

  
observe({
  curr_mode <<- mainmode		
  if(exists("dataSet")){
    if("shinymode" %in% names(dataSet)){
      mainmode <<- dataSet$shinymode
    }
  } 
  whichUI <- switch(mainmode,
                    time={time.anal.ui()},
                    stat={stat.anal.ui()})
  output$analUI <- renderUI({whichUI})
})

# -----------------
optUI <- function(){		
  selectInput('your.time', 'What end time do you want to pick?', choices = get_times(patdb))
}

observeEvent(input$exp_type,{
  req(input$exp_type)
  modes <- strsplit(x=input$exp_type, split = '\\.')[[1]]
  mainmode <<- modes[[1]]
  submode <<- modes[[2]]
  # ---------------------
  if(mainmode == "time" & submode != "standard"){
    output$exp_opt <- renderUI({optUI()})
  }
})


output$pos_add_tab <- DT::renderDataTable({
    # -------------
    datatable(pos_adducts,
              selection = list(mode = 'multiple', selected = c(1, 2, nrow(pos_adducts)), target = 'row'),
              options = list(pageLength = 5, dom = 'tp')
              , rownames = F)
  })

output$neg_add_tab <- DT::renderDataTable({
  # -------------
  datatable(neg_adducts,
            selection = list(mode = 'multiple', selected = c(1, 2, nrow(neg_adducts)), target = 'row'),
            options = list(pageLength = 5, dom = 'tp')
            , rownames = F)
})

# ===== PACKAGE LOADING ====

observeEvent(input$nav_general, {
  load.necessarities(input$nav_general)
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
    output$package_tab <- DT::renderDataTable({
      # -------------
      datatable(get.package.table(input$nav_general),
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
  
# ================ observers? =================

observe({  
  shinyDirChoose(input, "get_db_dir", roots = getVolumes(), session = session)
  given_dir <- input$get_db_dir$path
  if(is.null(given_dir)) return()
  options$db_dir <<- paste(c(given_dir), collapse="/")
  output$curr_db_dir <- renderText(options$db_dir)
  # --- connect ---
  setOption(".conf", "db_dir", options$db_dir)
  options <<- getOptions(".conf")
})

observe({  
  shinyDirChoose(input, "get_work_dir", roots = getVolumes(), session = session)
  given_dir <- input$get_work_dir$path
  if(is.null(given_dir)) return()
  options$work_dir <<- paste(c(given_dir), collapse="/")
  output$options$work_dir <- renderText(options$work_dir)
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
  output$ramp_plot <- renderPlotly({
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
                              "ryw"=heat.colors)
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
    # --- render ---
    plot_ly(z = volcano, 
          colors = color.function(100), 
          type = "heatmap",
          showscale=FALSE)  %>%
    layout(xaxis = ax, yaxis = ax)
  })
})


color.pickers <- reactive({
  req(dataSet)
  # -------------
  switch(mainmode,
         time = {
           lbl.fac <- if(dataSet$facA.lbl == "Time") "facB" else "facA"
           facs <- levels(dataSet[[lbl.fac]])
           },
         stat = {
           facs <- levels(dataSet$filt.cls)
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

output$colourPickers <- renderUI({color.pickers()})

color.vec <- reactive({
  req(dataSet)
  # ----------
  switch(mainmode,
         time = {
           if("facA" %not in% names(dataSet) & "facB" %not in% names(dataSet)) return(c("Red", "Green"))
           lbl.fac <- if(dataSet$facA.lbl == "Time") "facB" else "facA"
           facs <- dataSet[[lbl.fac]]
         },
         stat = {
           facs <- dataSet$filt.cls
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
  output$adduct_tab <- renderRHandsontable({
    if (!is.null(DF))
      rhandsontable(DF, stretchH = "all", useTypes = TRUE)
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

output$adduct_tab <- renderRHandsontable({
  DF = adduct_tab_data()
  if (!is.null(DF))
    rhandsontable(DF, stretchH = "all", useTypes = TRUE)
})

observe({
  DF = adduct_tab_data()
  shinyFileSave(input, "save_adducts", roots=getVolumes(), session=session)
  fileinfo <- parseSavePath(getVolumes(), input$save_adducts)
  if (nrow(fileinfo) > 0) {
    switch(fileinfo$type,
           csv = fwrite(file = fileinfo$datapath, x = DF)
    )
    }
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
    load(input$outlist_neg$datapath, verbose=T)
    load(input$outlist_pos$datapath, verbose=T)
    setProgress(.50,message = "Creating experiment database file...")
    build.pat.db(patdb,
                 ppm = ppm,
                 poslist = outlist_pos_renamed,
                 neglist = outlist_neg_renamed,
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
  req(mainmode)
  req(submode)
  # ---------
  withProgress({
    setProgress(1/4)
    # create csv
    tbl = get.csv(patdb,
                 time.series = if(mainmode == "time") T else F,
                 exp.condition = input$exp_var,
                 group_adducts = T,
                 group_by = input$group_by,
                 which_dbs = add_search_list,
                 which_adducts = selected_adduct_list
                 )
    # --- experimental stuff ---
    switch(mainmode,
           time = {
             switch(submode, 
                    standard = { tbl.adj <- tbl }, #mode stays the same
                    custom = { 
                      tbl.adj <<- unique(tbl[Time == input$your.time, -"Time"])
                      mainmode <<- "stat"
                      tbl.adj$Sample <- gsub(tbl.adj$Sample,pattern = "_T.*", replacement = "")
                    },
                    subtract = { 
                      uniq.samples <- unique(tbl$Sample)
                      table.base <- tbl[Time == input$your.time,c(1,3)][on=uniq.samples]
                      time.end <- tbl[Time == input$your.time,][,4:ncol(tbl),on=uniq.samples]
                      time.begin <- tbl[Time == min(as.numeric(tbl$Time)),][,4:ncol(tbl),on=uniq.samples]
                      tbl.adj <<- cbind(table.base, 
                                      time.end - time.begin) 
                      tbl.adj$Sample <- gsub(tbl.adj$Sample,pattern = "_T.*", replacement = "")
                      # --- change mode ---
                      mainmode <<- "stat"
                    })
             },
           stat = print("nothing to do for now"))
        # -------------------------
    # save csv
    setProgress(2/4)
    csv_loc <<- file.path(options$work_dir, paste0(options$proj_name,".csv"))
    fwrite(tbl.adj, csv_loc, sep="\t")
    # --- overview table ---
    setProgress(3/4)
    output$csv_tab <- DT::renderDataTable({
      overview_tab <- if(mainmode == "time"){
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
      datatable(overview_tab, 
                selection = 'single',
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE))
    })
  })
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
})
  
# ===================== METABOSTART ========================

observeEvent(input$check_excel, {
  # get excel table stuff.
  updateSelectInput(session, "exp_var",
                    choices = get_exp_vars()
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
  switch(mainmode,
         time = {
           InitDataObjects("pktable",
                           "ts",
                           FALSE)
           SetDesignType("time")
           Read.TextData(csv_loc,
                         "rowts",
                         "disc")
           },
         stat = {
           InitDataObjects("pktable",
                           "stat",
                           FALSE)
           Read.TextData(csv_loc, 
                         "rowu", 
                         "disc")
           }
         )
  SanityCheckData()
  RemoveMissingPercent(percent = 0.5)
  ImputeVar(method = "min")
  ReplaceMin()
  FilterVariable(filter = input$filt_type,
                 qcFilter = "F", 
                 rsd = 25)
  # # ---- here facA/facB disappears?? ---
  GetPrenormSmplNms()
  GetPrenormFeatureNms()
  GetPrenormClsNms()
  UpdateGroupItems()
  UpdateSampleItems()
  UpdateFeatureItems()
  # # ------------------------------------
  setProgress(.2)
  Normalization(rowNorm = input$norm_type,
                transNorm = input$trans_type,
                scaleNorm = input$scale_type,
                ref = input$ref_var)
  dataSet$grouping <<- input$group_by
  dataSet$shinymode <<- mainmode
  if(dataSet$grouping == "pathway"){dataSet$norm <- abs(dataSet$norm)}
  setProgress(.3)
  output$var_norm_plot <- renderPlot(PlotNormSummary())
  setProgress(.4)
  output$samp_norm_plot <- renderPlot(PlotSampleNormSummary())
  dataSet$adducts <<- selected_adduct_list
  # --- ABS-ify if grouped by pathway --- (use glmnet ish thing later?) ---
})
})

observeEvent(input$nav_time, {
  if(input$nav_general != "analysis") return(NULL)
    # get excel table stuff.
  switch(input$nav_time,
         ipca = {
           # --------------------
           output$plot_ipca <- renderPlotly({
             if("ipca" %not in% names(analSet)){
               iPCA.Anal(file.path(options$work_dir, "ipca_3d_0_.json"))
             }
             fac.lvls <- unique(analSet$ipca$score[[input$ipca_factor]])
             chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
             # ---------------
             df <- t(as.data.frame(analSet$ipca$score$xyz))
             x <- gsub(input$ipca_x,pattern = "\\(.*$", replacement = "")
             y <- gsub(input$ipca_y,pattern = "\\(.*$", replacement = "")
             z <- gsub(input$ipca_z,pattern = "\\(.*$", replacement = "")
             x.num <- as.numeric(gsub(x, pattern = "PC", replacement = ""))
             y.num <- as.numeric(gsub(y, pattern = "PC", replacement = ""))
             z.num <- as.numeric(gsub(z, pattern = "PC", replacement = ""))
             # --- render! ---
             plot_ly(hoverinfo = 'text',
                     text = analSet$ipca$score$name ) %>%
               add_trace(
                 x = df[,1], 
                 y = df[,2], 
                 z = df[,3], 
                 type = "scatter3d",
                 color= analSet$ipca$score[[input$ipca_factor]], colors=chosen.colors
               ) %>%  layout(scene = list(
                 xaxis = list(
                   title = analSet$ipca$score$axis[x.num]),
                 yaxis = list(
                   title = analSet$ipca$score$axis[y.num]),
                 zaxis = list(
                   title = analSet$ipca$score$axis[z.num])))
           })
          
           # but perform first...
           updateSelectInput(session, "ipca_factor",
                             choices = grep(names(json_pca$score), pattern = "^fac[A-Z]", value = T))
           # ------------------
           
           ipca.table <- as.data.table(data.table(gsub(analSet$ipca$score$axis, pattern = "\\(.*$| ", replacement = ""),
                                                 gsub(analSet$ipca$score$axis, pattern = "PC\\d*| |\\(|%\\)", replacement = "")),
                                      keep.rownames = T)
           colnames(ipca.table) <- c("Principal Component", "% variance")
           output$ipca_tab <- DT::renderDataTable({
             req(ipca.table)
             # -------------
             datatable(ipca.table, 
                       selection = 'single',
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           })
         },
         meba = {
           if("MB" %not in% names(analSet)){
            performMB(10, dir=options$work_dir)
           }
           output$meba_tab <- DT::renderDataTable({
             # -------------
             datatable(analSet$MB$stats, 
                       selection = 'single',
                       colnames = c("Compound", "Hotelling/T2 score"),
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           })
         },
         asca = {
           if("asca" %not in% names(analSet)){
             Perform.ASCA(1, 1, 2, 2)
             CalculateImpVarCutoff(0.05, 0.9, dir=options$work_dir)
           }
           output$asca_tab <- DT::renderDataTable({
             # -------------
             datatable(analSet$asca$sig.list$Model.ab, 
                       selection = 'single',
                       colnames = c("Compound", "Leverage", "SPE"),
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           })

         })
  })

observeEvent(input$tab_stat, {
  if(input$nav_general != "analysis") return(NULL)
  # get excel table stuff.
  switch(input$tab_stat,
         pca = {
           PCA.Anal()
           output$plot_pca <- renderPlotly({
             df <- analSet$pca$x
             x <- input$pca_x
             y <- input$pca_y
             z <- input$pca_z
             x.var <- round(analSet$pca$variance[x] * 100.00, digits=1)
             y.var <- round(analSet$pca$variance[y] * 100.00, digits=1)
             z.var <- round(analSet$pca$variance[z] * 100.00, digits=1)
             fac.lvls <- unique(dataSet$filt.cls)
             chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
             # ---------------
             plot_ly(hoverinfo = 'text',
                     text = rownames(df) ) %>%
               add_trace(
                 x = analSet$pca$x[,x], 
                 y = analSet$pca$x[,y], 
                 z = analSet$pca$x[,z], 
                 type = "scatter3d",
                 color= dataSet$filt.cls, colors=chosen.colors
               ) %>%  layout(scene = list(
                 xaxis = list(
                   title = fn$paste("$x ($x.var %)")),
                 yaxis = list(
                   title = fn$paste("$y ($y.var %)")),
                 zaxis = list(
                   title = fn$paste("$z ($z.var %)"))))
           })
           pca.table <- as.data.table(round(analSet$pca$variance * 100.00,
                                            digits = 2),
                                      keep.rownames = T)
           colnames(pca.table) <- c("Principal Component", "% variance")
           output$pca_tab <- DT::renderDataTable({
             req(pca.table)
             # -------------
             datatable(pca.table, 
                       selection = 'single',
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           })},
         plsda = {
           PLSR.Anal()
           PLSDA.CV("L",5, "Q2")
           plsda.table <- as.data.table(round(analSet$plsr$Xvar 
                                              / analSet$plsr$Xtotvar 
                                              * 100.0,
                                              digits = 2),
                                        keep.rownames = T)
           colnames(plsda.table) <- c("Principal Component", "% variance")
           plsda.table[, "Principal Component"] <- paste0("PC", 1:nrow(plsda.table))
           # -----------
           output$plot_plsda <- renderPlotly({
             x <- input$plsda_x
             y <- input$plsda_y
             z <- input$plsda_z
             x.var <- plsda.table[`Principal Component` == x, 
                                  `% variance`]
             y.var <- plsda.table[`Principal Component` == y, 
                                  `% variance`]
             z.var <- plsda.table[`Principal Component` == z, 
                                  `% variance`]
             fac.lvls <- unique(dataSet$filt.cls)
             colnames(analSet$plsr$Yscores) <- paste0("PC", 1:ncol(analSet$plsr$Yscores))
             chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
             # ---------------
             plot_ly(hoverinfo = 'text',
                     text = rownames(dataSet$norm)) %>%
               add_trace(
                 x = analSet$plsr$Yscores[,x], 
                 y = analSet$plsr$Yscores[,y], 
                 z = analSet$plsr$Yscores[,z], 
                 type = "scatter3d",
                 color = dataSet$filt.cls, colors=chosen.colors
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
  output$plsda_tab <- DT::renderDataTable({
    req(plsda.table)
    # -------------
    datatable(plsda.table, 
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
  })
  # -------------
  output$plsda_vip_tab <- DT::renderDataTable({
    colnames(analSet$plsda$vip.mat) <- paste0("PC", 1:ncol(analSet$plsda$vip.mat))
    compounds_pc <- as.data.table(analSet$plsda$vip.mat,keep.rownames = T)
    ordered_pc <- setorderv(compounds_pc, input$plsda_vip_cmp, -1)
    plsda_tab <<- cbind(ordered_pc[1:50, c("rn")], 
                        Rank=c(1:50))
    # -------------
    datatable(plsda_tab,
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
    })
         },
         heatmap = {
           req(input$heatmode)
           req(input$color_ramp)
           req(analSet$tt)
           req(analSet$fc)
           # --- render ---
           output$heatmap <- renderPlotly({
             x <- if(input$heatmode == TRUE){
               dataSet$norm[,names(analSet$tt$inx.imp[analSet$tt$inx.imp == TRUE])]
             } else{dataSet$norm[,names(analSet$fc$inx.imp[analSet$fc$inx.imp == TRUE])]}
             final_matrix <- t(x)
             translator <- data.table(Sample=rownames(dataSet$norm),Group=dataSet$prenorm.cls)
             group_assignments <- translator[,"Group",on=colnames(final_matrix)]$Group
             hm_matrix <<- heatmapr(final_matrix, 
                                    Colv=T, 
                                    Rowv=T,
                                    col_side_colors = group_assignments,
                                    k_row = NA)
             heatmaply(hm_matrix,
                       Colv=F, Rowv=F,
                       branches_lwd = 0.3,
                       margins = c(60,0,NA,50),
                       colors = color.function(256),
                       col_side_palette = function(n){
                         if(n == length(color.vec())) color.vec() else rainbow(n)
                       },
                       subplot_widths = c(.9,.1),
                       subplot_heights =  c(.1,.05,.85),
                       column_text_angle = 90)
           })
           },
         tt = {
           Ttests.Anal(F, 0.05, FALSE, TRUE)
           output$tt_tab <- DT::renderDataTable({
             # -------------
             datatable(analSet$tt$sig.mat, 
                       selection = 'single',
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             
           })
           output$tt_overview_plot <- renderPlotly({
             # --- ggplot ---
             ggPlotTT(color.function, 20)
           })
           },
         fc = {
           FC.Anal.unpaired(2.0, 1)
           output$fc_tab <- DT::renderDataTable({
             # -------------
             datatable(analSet$fc$sig.mat, 
                       selection = 'single',
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             
           })
           output$fc_overview_plot <- renderPlotly({
             # --- ggplot ---
             ggPlotFC(color.function, 20)
           })
           },
         volc = {
           Volcano.Anal(FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
           output$volc_tab <- DT::renderDataTable({
             # -------------
             datatable(analSet$volc$sig.mat, 
                       selection = 'single',
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             
           })
           output$volc_plot <- renderPlotly({
             # --- ggplot ---
             ggPlotVolc(color.function, 20)
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
  print(add_search_list)
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
                        "asca", 
                        "meba",
                        "plsda_vip")

lapply(res.update.tables, FUN=function(table){
  observeEvent(input[[paste0(table, "_tab_rows_selected")]], {
    curr_row = input[[paste0(table, "_tab_rows_selected")]]
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(curr_row)) return()
    rownames(analSet$tt$sig.mat)
    curr_cpd <<- data.table(switch(table,
                                     tt = analSet$tt$sig.mat,
                                     fc = analSet$fc$sig.mat,
                                     asca = analSet$asca$sig.list$Model.ab,
                                     meba = analSet$MB$stats,
                                     plsda_vip = plsda_tab)
                              , keep.rownames = T)[curr_row, rn]
    output$curr_cpd <- renderText(curr_cpd)
    output[[paste0(table, "_specific_plot")]] <- renderPlotly({
      # --- ggplot ---
      if(table == 'meba'){
        ggplotMeba(curr_cpd, draw.average, cols = color.vec())
      }else{
        ggplotSummary(curr_cpd, cols = color.vec())
      }
    })
    if(input$autosearch & length(db_search_list > 0)){
      match_list <- lapply(db_search_list, FUN=function(match.table){
        get_matches(curr_cpd, match.table, searchid=dataSet$grouping)
      })
      match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
      output$match_tab <- DT::renderDataTable({
        datatable(if(dataSet$grouping == 'pathway') match_table else{match_table[,-"Description"]},
                  selection = 'single',
                  autoHideNavigation = T,
                  options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
      })  
    }
  })
})

observe({
  d <- event_data("plotly_click")
  req(d)
  # -----------------------------
  switch(input$tab_stat,
         tt = {
           if('key' %not in% colnames(d)) return(NULL)
           if(d$key %not in% names(analSet$tt$p.value)) return(NULL)
           curr_cpd <<- d$key
           output$tt_specific_plot <- renderPlotly({
             # --- ggplot ---
             ggplotSummary(curr_cpd, cols = color.vec())
           })
           },
         fc = {
           if('key' %not in% colnames(d)) return(NULL)
           if(d$key %not in% names(analSet$fc$fc.log)) return(NULL)
           curr_cpd <<- d$key
           output$fc_specific_plot <- renderPlotly({
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
           if(d$key %not in% rownames(analSet$volcano$sig.mat)) return(NULL)
           curr_cpd <<- d$key
         })
  # ----------------------------
  output$curr_cpd <- renderText(curr_cpd)
  if(input$autosearch & length(db_search_list) > 0){
    match_list <- lapply(db_search_list, FUN=function(match.table){
      get_matches(curr_cpd, match.table, searchid=input$group_by)
    })
    match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
    output$match_tab <- DT::renderDataTable({
      datatable(if(dataSet$grouping == 'pathway') match_table else{match_table[,-"Description"]},
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
    match_table <<- unique(as.data.table(rbindlist(match_list))[Name != ""])
    output$match_tab <- DT::renderDataTable({
      datatable(if(dataSet$grouping == 'pathway') match_table else{match_table[,-"Description"]},
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
  curr_def <<- if(dataSet$grouping == 'pathway') "Unavailable" else(match_table[curr_row,'Description'])
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
  output$browse_tab <- DT::renderDataTable({
    datatable(browse_table[,-c("Description", "Charge")],
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
  output$hits_tab <- DT::renderDataTable({
    datatable(hits_table,
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
  output$meba_specific_plot <- renderPlotly({ggplotMeba(curr_cpd, draw.average, cols = color.vec())})
  output$asca_specific_plot <- renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
  output$fc_specific_plot <- renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
  output$tt_specific_plot <- renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
  output$plsda_specific_plot <- renderPlotly({ggplotSummary(curr_cpd, cols = color.vec())})
  })

# --- ON CLOSE ---
session$onSessionEnded(function() {
  if(any(!is.na(session_cl))) stopCluster(session_cl)
  #save(dataSet, analSet, file = file.path(options$work_dir, paste0(options$proj_name, ".RData")))
})
})
