function(input, output, session) {

library(shiny)
  
# ===== defaults =====


mz <<- NA
session_cl <<- NA
patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))

packages <<- c("data.table", "DBI", "RSQLite", "ggplot2", "minval", "enviPat",
               "plotly", "parallel", "shinyFiles", "curl", "httr", "pbapply", 
               "sqldf", "plyr", "ChemmineR", "gsubfn", "stringr", "plotly", "heatmaply",
               "reshape2", "XML", "xlsx", "colourpicker", "DT","Rserve", "ellipse", 
               "scatterplot3d","pls", "caret", "lattice", "compiler",
               "Cairo", "randomForest", "e1071","gplots", "som", "xtable",
               "RColorBrewer", "xcms","impute", "pcaMethods","siggenes",
               "globaltest", "GlobalAncova", "Rgraphviz","KEGGgraph",
               "preprocessCore", "genefilter", "pheatmap", "igraph",
               "RJSONIO", "SSPA", "caTools", "ROCR", "pROC", "sva", "rJava")

options <- getOptions(".conf")

default.text <- list(list(name='options$work_dir',text=options$work_dir),
                     list(name='curr_db_dir',text=options$db_dir),
                     list(name='ppm',text=options$ppm),
                     list(name='analUI',text="Please choose an analysis mode!"),
                     list(name='proj_name',text=options$proj_name)
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
                      plotlyOutput("plot_ipca" ),
                      selectInput("ipca_factor", label = "Color based on:", choices =list("Time"="facA",
                                                                                          "Experimental group"="facB"))
             ),
             # =================================================================================
             tabPanel("MEBA", value="meba", 
                      fluidRow(plotlyOutput('meba_plot')),
                      fluidRow(div(DT::dataTableOutput('meba_tab'),style='font-size:80%'))),
             # =================================================================================
             tabPanel("ASCA", value="asca",
                      fluidRow(plotlyOutput('asca_plot')),
                      fluidRow(div(DT::dataTableOutput('asca_tab'),style='font-size:80%'))
             )
             # =================================================================================
  )
})

stat.anal.ui <- reactive({
  navbarPage("Standard analysis", id="tab_stat",
             tabPanel("PCA", value = "pca", 
                      navbarPage("Explore", id="tab_pca",
                                 tabPanel("PCA", icon=icon("eye"), helpText(plotlyOutput("plot_pca"),
                                                                                 selectInput("pca_x", label = "X axis:", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "Other")),
                                                                                 selectInput("pca_y", label = "Y axis:", choices =c("PC1", "PC2", "PC3", "PC4", "PC5","Other")),
                                                                                 selectInput("pca_z", label = "Z axis:", choices =c("PC1", "PC2", "PC3", "PC4", "PC5","Other")))),
                                 
                                 tabPanel("PLS-DA", icon=icon("bar-chart-o"),
                                          helpText("placeholder")
                                 )
                      )
             ),
             # =================================================================================
             tabPanel("Heatmap",
                      plotlyOutput("heatmap", height='600px', width='900px')
             ),
             # =================================================================================
             tabPanel("T-test", value="tt", 
                      fluidRow(plotOutput('tt_plot')),
                      fluidRow(div(DT::dataTableOutput('tt_tab'),style='font-size:80%'))),
             tabPanel("Fold-change", value="fc",
                      fluidRow(plotOutput('fc_plot')),
                      fluidRow(div(DT::dataTableOutput('fc_tab'),style='font-size:80%'))),
             # =================================================================================
             tabPanel("Volcano", value="volc",
                      fluidRow(plotOutput('volc_plot')),
                      fluidRow(div(DT::dataTableOutput('volc_tab'),style='font-size:80%')))
  )
})

observeEvent(input$exp_type,{
  req(input$exp_type)
  modes <- strsplit(x=input$exp_type, split = '\\.')[[1]]
  print(modes)
  mainmode <<- modes[[1]]
  submode <<- modes[[2]]
  # ---------------------
  optUI <- function(){
    selectInput('your.time', 'What end time do you want to pick?', choices = get_times(patdb))
  }
  if(mainmode == "time" & submode != "standard"){
    output$exp_opt <- renderUI({optUI()})
  }
})

# ===== PACKAGE LOADING ====

observeEvent(input$nav_general, {
  load.necessarities(input$nav_general)
})

# =================================================

images <- list(list(name = 'cute_package', path = 'backend/img/new-product.png', dimensions = c(80, 80)),
               list(name = 'umc_logo', path = 'backend/img/umclogo.jpg', dimensions = c(120, 100)),
               list(name = 'hmdb_logo', path = 'backend/img/hmdblogo.png', dimensions = c(150, 100)),
               list(name = 'chebi_logo', path = 'backend/img/chebilogo.png', dimensions = c(100, 100)),
               list(name = 'pubchem_logo', path = 'backend/img/pubchemlogo.png', dimensions = c(140, 100)),
               list(name = 'pos_icon', path = 'backend/img/handpos.png', dimensions = c(120, 120)),
               list(name = 'neg_icon', path = 'backend/img/handneg.png', dimensions = c(120, 120)),
               list(name = 'excel_icon', path = 'backend/img/excel.png', dimensions = c(120, 120)),
               list(name = 'db_icon', path = 'backend/img/office.png', dimensions = c(100, 100))
               )

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
    filename <- normalizePath(file.path('backend/img/yes.png'))
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
  default.colors <- rainbow(length(facs))
  # -------------
  lapply(seq_along(facs), function(i) {
    colourInput(inputId = paste("col", i, sep="_"), 
                label = paste("Choose colour for", facs[i]), 
                value = default.colours[i]) 
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
    print(i)
    input[[paste("col", i, sep="_")]]
  }))
})

# ======================== DB CHECK ============================



# --- check for db files ---

observeEvent(input$check_umc,{
  db_folder_files <- list.files(options$db_dir)
  is.present <- "internal.full.db" %in% db_folder_files | "noise.full.db" %in% db_folder_files
  check_pic <- if(is.present) "yes.png" else "no.png"
  output$umc_check <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('backend/img', check_pic))
    # Return a list containing the filename and alt text
    list(src = filename, width = 70,
         height = 70)
  }, deleteFile = FALSE)
})

observeEvent(input$check_hmdb,{
  db_folder_files <- list.files(options$db_dir)
  is.present <- "hmdb.full.db" %in% db_folder_files
  check_pic <- if(is.present) "yes.png" else "no.png"
  output$hmdb_check <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('backend/img', check_pic))
    # Return a list containing the filename and alt text
    list(src = filename, width = 70,
         height = 70)
  }, deleteFile = FALSE)
})

observeEvent(input$check_chebi,{
  db_folder_files <- list.files(options$db_dir)
  is.present <- "chebi.full.db" %in% db_folder_files
  check_pic <- if(is.present) "yes.png" else "no.png"
  output$chebi_check <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('backend/img', check_pic))
    # Return a list containing the filename and alt text
    list(src = filename, width = 70,
         height = 70)
  }, deleteFile = FALSE)
})

observeEvent(input$check_pubchem,{
  db_folder_files <- list.files(options$db_dir)
  is.present <- "pubchem.full.db" %in% db_folder_files
  check_pic <- if(is.present) "yes.png" else "no.png"
  output$pubchem_check <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('backend/img', check_pic))
    # Return a list containing the filename and alt text
    list(src = filename, width = 70,
         height = 70)
  }, deleteFile = FALSE)
})

# --- build db ---

observeEvent(input$build_umc, {
  # ---------------------------
  withProgress({
    setProgress(message = "Working...")
    build.base.db("internal", outfolder=options$db_dir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("internal", 
                      outfolder=options$db_dir,
                      adduct.table = wkz.adduct.confirmed, 
                      cl=session_cl, 
                      fetch.limit=10)
    build.base.db("noise", outfolder=options$db_dir) # does both because its a special db
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_hmdb,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("hmdb", outfolder=options$db_dir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("hmdb", 
                      outfolder=options$db_dir, 
                      adduct.table = wkz.adduct.confirmed, 
                      cl=session_cl, 
                      fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_chebi,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("chebi", outfolder=options$db_dir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("chebi", outfolder=options$db_dir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_pubchem,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("pubchem", outfolder=options$db_dir, cl = session_cl)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("pubchem", outfolder=options$db_dir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

# ================== DATA IMPORT ===========================

observeEvent(input$create_db,{
  # --------------------
  print(options$work_dir)
  patdb <<- file.path(options$work_dir, paste0(options$proj_name, ".db"))
  print(patdb)
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
    filename <- normalizePath('backend/img/yes.png')
    # Return a list containing the filename and alt text
    list(src = filename, width = 20,
         height = 20)
  }, deleteFile = FALSE)
})

# ==================== CREATE CSV =======================

observeEvent(input$import_csv, {
  req(input$pat_csv)
  # -----------------------------------
  csv_loc <<- input$pat_csv$datapath
  output$csv_upload_check <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath('backend/img/yes.png')
    # Return a list containing the filename and alt text
    list(src = filename, width = 20,
         height = 20)
  }, deleteFile = FALSE)
  mz <<- fread(csv_loc, sep="\t", header = T)
  # --- detect mode from csv ---
  mainmode <<- if("Time" %not in% names(mz)[1:5]) "stat" else "time"
  # ----------------------------
  overview_tab <- if(mainmode == "time"){
    t(data.table(keep.rownames = F,
                 Mz = ncol(mz) - 3,
                 Samples = nrow(mz),
                 Times = length(unique(mz$Time)),
                 Groups = length(unique(mz$Label))
    ))
  } else{
    t(data.table(keep.rownames = F,
                 Mz = ncol(mz) - 3,
                 Samples = nrow(mz),
                 Groups = length(unique(mz$Label))
    )) 
  }
  output$csv_tab <- DT::renderDataTable({
    datatable(overview_tab, 
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE))
    })
  whichUI <- switch(mainmode,
                    time={time.anal.ui()},
                    stat={stat.anal.ui()})
  output$analUI <- renderUI({whichUI})
  })

output$csv_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/office.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=100,
       height=100)
}, deleteFile = FALSE)
  
observeEvent(input$create_csv, {
  req(options$proj_name)
  req(options$work_dir)
  req(mainmode)
  req(submode)
  req(input$your.time)
  # ---------
  withProgress({
    setProgress(1/4, "Creating csv file for MetaboAnalyst...")
    # create csv
    patdb
    mz = get.csv(patdb,
                 time.series = if(mainmode == "time") T else F,
                 exp.condition = input$exp_var)
    # --- experimental stuff ---
    switch(mainmode,
           time = {
             switch(submode, 
                    standard = { mz.adj <- mz }, #mode stays the same
                    custom = { 
                      mz.adj <- unique(mz[Time == input$your.time, -"Time"])
                      mainmode <<- "stat"
                    },
                    subtract = { 
                      uniq.samples <- unique(mz$Sample)
                      table.base <- mz[Time == input$your.time,c(1,3)][on=uniq.samples]
                      time.end <- mz[Time == input$your.time,][,4:ncol(mz),on=uniq.samples]
                      time.begin <- mz[Time == min(as.numeric(mz$Time)),][,4:ncol(mz),on=uniq.samples]
                      mz.adj <- cbind(table.base, 
                                      time.end - time.begin) 
                      # --- change mode ---
                      mainmode <<- "stat"
                    })
             },
           stat = print("nothing to do for now"))
        # -------------------------
    # save csv
    setProgress(2/4, "Writing csv file...")
    csv_loc <<- file.path(options$work_dir, paste0(options$proj_name,".csv"))
    fwrite(mz.adj, csv_loc, sep="\t")
    # --- overview table ---
    setProgress(3/4, "Creating overview table...")
    
    overview_tab <- if(mode == "time"){
      t(data.table(keep.rownames = F,
                   Mz = ncol(mz.adj) - 3,
                   Samples = nrow(mz.adj),
                   Times = length(unique(mz.adj$Time)),
                   Groups = length(unique(mz.adj$Label))
      ))
    } else{
      t(data.table(keep.rownames = F,
                   Mz = ncol(mz.adj) - 3,
                   Samples = nrow(mz.adj),
                   Groups = length(unique(mz.adj$Label))
      )) 
      }
    output$csv_tab <- DT::renderDataTable({
      datatable(overview_tab, 
                selection = 'single',
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE))
    })
    whichUI <- switch(mainmode,
                      time={time.anal.ui()},
                      stat={stat.anal.ui()})
    output$analUI <- renderUI({whichUI})
  })
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
  if(input$norm_type == "CompNorm"){
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
  # get excel table stuff.
  updateSelectInput(session, "ref_var",
                    choices = get_ref_vars(input$exp_var)
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
  sourceAll(file.path("backend", 
                      "scripts", 
                      "metaboanalyst"))
  setProgress(.1, "Applying your settings...")
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
  FilterVariable(input$filt_type, 
                 "F", 
                 25)
  # # ---- here facA/facB disappears?? ---
  GetPrenormSmplNms()
  GetPrenormFeatureNms()
  GetPrenormClsNms()
  UpdateGroupItems()
  UpdateSampleItems()
  UpdateFeatureItems()
  # # ------------------------------------
  setProgress(.2, "Normalizing data...")
  Normalization(rowNorm = input$norm_type,
                transNorm = input$trans_type,
                scaleNorm = input$scale_type,
                ref = input$ref_var)
  setProgress(.3, "Plotting variable overview...")
  output$var_norm_plot <- renderPlot(PlotNormSummary())
  setProgress(.4, "Plotting sample overview...")
  output$samp_norm_plot <- renderPlot(PlotSampleNormSummary())
})
})

observeEvent(input$nav_time, {
  if(input$nav_general != "analysis") return(NULL)
    # get excel table stuff.
  switch(input$nav_time,
         ipca = {
           if("ipca" %not in% names(analSet)){
             iPCA.Anal(file.path(options$work_dir, "ipca_3d_0_.json"))
             analSet[["ipca"]] <- "Done!"             
           }
           json_pca <- fromJSON(file.path(options$work_dir, "ipca_3d_0_.json"))
           req(json_pca)
           # --------------------
           output$plot_ipca <- renderPlotly({
             req(json_pca)
             fac.lvls <- unique(json_pca$score[[input$ipca_factor]])
             chosen.colors <- if(length(fac.lvls) == length(color.vec())) color.vec() else rainbow(length(fac.lvls))
             # ---------------
             df <- t(as.data.frame(json_pca$score$xyz))
             plot_ly(hoverinfo = 'text',
                     text = json_pca$score$name ) %>%
               add_trace(
                 x = df[,1], 
                 y = df[,2], 
                 z = df[,3], 
                 type = "scatter3d",
                 color= json_pca$score[[input$ipca_factor]], colors=chosen.colors
               ) %>%  layout(scene = list(
                 xaxis = list(
                   title = json_pca$score$axis[1]),
                 yaxis = list(
                   title = json_pca$score$axis[2]),
                 zaxis = list(
                   title = json_pca$score$axis[3])))
           })
           # but perform first...
           updateSelectInput(session, "ipca_factor",
                             choices = grep(names(json_pca$score), pattern = "^fac[A-Z]", value = T))
         },
         meba = {
           if("MB" %not in% names(analSet)){
            performMB(10, dir=options$work_dir)
           }
           meba.table <<- read.csv(file.path(options$work_dir, 'meba_sig_features.csv'))
           output$meba_tab <- DT::renderDataTable({
             req(meba.table)
             # -------------
             datatable(meba.table, 
                       selection = 'single',
                       colnames = c("Mass/charge", "Hotelling/T2 score"),
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           })
         },
         asca = {
           if("asca" %not in% names(analSet)){
             Perform.ASCA(1, 1, 2, 2)
              CalculateImpVarCutoff(0.05, 0.9, dir=options$work_dir)
            }
           asca.table <<- read.csv(file.path(options$work_dir,'Sig_features_Model_ab.csv'))
           output$asca_tab <- DT::renderDataTable({
             req(asca.table)
             # -------------
             datatable(asca.table, 
                       selection = 'single',
                       colnames = c("Mass/charge", "Leverage", "SPE"),
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
           pca.table
           output$pca_tab <- DT::renderDataTable({
             req(pca.table)
             # -------------
             datatable(pca.table, 
                       selection = 'single',
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
           })},
         heatmap = {   NULL
           },
         tt = {
           Ttests.Anal(F, 0.05, FALSE, TRUE)
           tt.table <<- as.data.table(read.csv(file.path(options$work_dir,'t_test.csv')),keep.rownames = F)
           output$tt_tab <- DT::renderDataTable({
             req(tt.table)
             # -------------
             datatable(tt.table[,c(1,3)], 
                       selection = 'single',
                       colnames = c("Mass/charge", "p-value"),
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             
           })
           output$tt_overview_plot <- renderPlotly({
             # --- ggplot ---
             ggPlotTT()
           })
           },
         fc = {
           FC.Anal.unpaired(2.0, 1)
           fc.table <<- as.data.table(read.csv(file.path(options$work_dir,'fold_change.csv')),keep.rownames = F)
           fc.table
           output$fc_tab <- DT::renderDataTable({
             # -------------
             datatable(fc.table, 
                       selection = 'single',
                       colnames = c("Mass/charge", "FC", "log2(FC)"),
                       autoHideNavigation = T,
                       options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
             
           })
           output$fc_overview_plot <- renderPlotly({
             # --- ggplot ---
             ggPlotFC()
           })
           })
})

observeEvent(input$heatmode,{
  print(input$heatmode)
  x <- switch(input$heatmode, 
              tt = dataSet$norm[,names(analSet$tt$inx.imp[analSet$tt$inx.imp == TRUE])],
              fc = dataSet$norm[,names(analSet$fc$inx.imp[analSet$fc$inx.imp == TRUE])]
  )
  final_matrix <- t(x)
  translator <- data.table(Sample=rownames(dataSet$norm),Group=dataSet$prenorm.cls)
  group_assignments <- translator[,"Group",on=colnames(final_matrix)]$Group
  assigned_colours <- sapply(group_assignments, FUN=function(group){
    rc = rainbow(2)
    as.factor(switch(as.character(group),
           Safe = rc[1],
           Challenge = rc[2])
           )
  })
  hm_matrix <<- heatmapr(final_matrix, 
                        Colv=T, 
                        Rowv=T,
                        col_side_colors = group_assignments,
                        k_row = NA)
  output$heatmap <- renderPlotly({
    heatmaply(hm_matrix,
              Colv=F, Rowv=F,
              branches_lwd = 0.3,
              margins = c(60,0,NA,50),
              colors = rainbow(256),
              col_side_palette = RdYlGn,
              subplot_widths = c(.9,.1),
              subplot_heights =  c(.1,.05,.85),
              column_text_angle = 90)
    
      })
})


observeEvent(input$tt_tab_rows_selected,{
  curr_row = input$tt_tab_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- tt.table[curr_row,X]
  print(curr_mz)
  output$tt_specific_plot <- renderPlotly({
    # --- ggplot ---
    ggplotSummary(curr_mz, cols = color.vec())
  })
})

observeEvent(input$fc_tab_rows_selected,{
  curr_row = input$fc_tab_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- fc.table[curr_row,X]
  print(curr_mz)
  output$fc_specific_plot <- renderPlotly({
    # --- ggplot ---
    ggplotSummary(curr_mz, cols = color.vec())
  })
})

observe({
  d <- event_data("plotly_click")
  if(length(d) == 0) return(NULL)
  switch(input$tab_stat,
         tt = {
           if(d$key %not in% names(analSet$tt$p.value)) return(NULL)
           curr_mz <<- d$key
           output$tt_specific_plot <- renderPlotly({
             # --- ggplot ---
             ggplotSummary(curr_mz, cols = color.vec())
           })
           },
         fc = {
           if(d$key %not in% names(analSet$fc$fc.log)) return(NULL)
           curr_mz <<- d$key
           output$fc_specific_plot <- renderPlotly({
             # --- ggplot ---
             ggplotSummary(curr_mz, cols = color.vec())
           })
         },
         heat = {
           if(d$y > length(hm_matrix$matrix$rows)) return(NULL)
           curr_mz <<- hm_matrix$matrix$rows[d$y]
             })
})
         


# ================ MEBA ========================

observeEvent(input$meba_tab_rows_selected,{
  curr_row = input$meba_tab_rows_selected
  draw.average = T
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- meba.table[curr_row,'X']
  output$meba_plot <- renderPlotly({
    # --- ggplot ---
    ggplotMeba(curr_mz, draw.average, cols = color.vec())
  })
})

# check for selected mz row
observeEvent(input$asca_tab_rows_selected,{
  curr_row = input$asca_tab_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- asca.table[curr_row,'X']
  output$asca_plot <- renderPlotly({ggplotSummary(curr_mz, cols = color.vec())})})

# --- find matches ---

output$find_mol_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/search.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=70,
       height=70)
}, deleteFile = FALSE)

observeEvent(input$search_mz,{
  req(input$checkGroup)
  req(curr_mz)
  # -------------------
  match_list <- lapply(input$checkGroup, FUN=function(match.table){
    get_matches(curr_mz, match.table)})
  match_table <<- unique(as.data.table(rbindlist(match_list))[Compound != ""])
  output$match_tab <- DT::renderDataTable({
    datatable(match_table[,-"Description"],
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
  })
})

observeEvent(input$match_tab_rows_selected,{
  curr_row = input$match_tab_rows_selected
  curr_row <<- input$match_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_def <<- match_table[curr_row,'Description']
  output$curr_definition <- renderText(curr_def$Description)
  })


observeEvent(input$browse_db,{
  req(input$checkGroup)
  # -------------------
  cpd_list <- lapply(input$checkGroup, FUN=function(match.table){
  browse_db(match.table)})
  # ------------------
  browse_table <<- unique(as.data.table(rbindlist(cpd_list)))
  output$browse_tab <- DT::renderDataTable({
    datatable(browse_table[,-c("Description", "Charge")],
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 10, 20), pageLength = 5))
  })
})

observeEvent(input$search_cpd,{
  req(input$checkGroup)
  req(input$browse_tab_rows_selected)
  # -------------------
  search_cmd <- browse_table[curr_row,c('Formula', 'Charge')]
  print(search_cmd)
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
  curr_mz <<- hits_table[curr_row, mzmed.pgrp]
  output$meba_plot <- renderPlotly({ggplotMeba(curr_mz, draw.average, cols = color.vec() )})
  output$asca_plot <- renderPlotly({ggplotSummary(curr_mz, cols = color.vec())})
  })

# --- ON CLOSE ---
session$onSessionEnded(function() {
  if(any(!is.na(session_cl))) stopCluster(session_cl)
})
}
