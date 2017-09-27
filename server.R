function(input, output, session) {

library(shiny)
  
# ===== defaults =====


output$exp_dir <- renderText(exp_dir)
output$proj_name <- renderText(proj_name)
output$curr_db_dir <- renderText(dbDir)
output$ppm <- renderText(ppm)
output$analUI <- renderUI({helpText("Please choose a mode")})
session_cl <<- NA
patdb <<- file.path(exp_dir, paste0(proj_name, ".db"))


time.anal.ui <- reactive({
  # -------------
  navbarPage("Time Series", id="nav_time",
             tabPanel("iPCA", value = "ipca", 
                      plotlyOutput("plot_ipca" ),
                      selectInput("ipca_factor", label = "Color based on:", choices =list("Time"="facA",
                                                                                          "Experimental group"="facB"))
             ),
             # =================================================================================
             #tabPanel("Heatmap",
             #         plotOutput("time_heat_plot", height='600px', width='600px')
             #        ),
             #tabPanel("MANOVA",
             #        ),
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
                      plotOutput("heatmap", height='600px', width='600px')
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
  whichUI <- switch(input$exp_type,
                       time={time.anal.ui()},
                       stat={stat.anal.ui()})
  print(whichUI)
  output$analUI <- renderUI({whichUI})
})

# ===== PACKAGE LOADING ====

observeEvent(input$nav_general, {
  print(input$nav_general)
  switch(input$nav_general,
         setup = {
           library(pacman)
           library(data.table)},
         database = {
           library(RSQLite)
           library(gsubfn)
           library(DBI)
           library(parallel)
           library(XML)
           library(minval)
           library(curl)
           library(enviPat)
           data(isotopes, package = "enviPat")
           if(any(is.na(session_cl))){
             session_cl <<- makeCluster(detectCores())
             clusterExport(session_cl, envir = .GlobalEnv, varlist = list(
               "isotopes",
               "subform.joanna", 
               "mergeform.joanna",
               "multiform.joanna",
               "check.ded.joanna",
               "data.table",
               "rbindlist",
               "isopattern"
             ))
             }
         },
         upload = {
           library(RSQLite)
           library(DBI)
           library(reshape2)
           library(data.table)
           library(xlsx)
         },
         document = {
           library(plotly)
           library(data.table)
           library(reshape2)},
         filter = {
           library(plotly)
           library(data.table)
           library(Cairo)
           library(preprocessCore)},
         analysis = {
           library(plotly)
           library(data.table)},
         options = {
           library(shinyFiles)
           library(colourpicker)
         })
})
  

# =================================================

output$cute_package <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath(file.path('backend/img/new-product.png'))
    # Return a list containing the filename and alt text
    list(src = filename, width = 80,
         height = 80)
  }, deleteFile = FALSE)

observeEvent(input$nav_general, {
  if(input$nav_general == "setup"){
    status <- sapply(packages, FUN=function(package){
      if(package %in% rownames(installed.packages())) "Yes" else "No"
    })
    version <- sapply(packages, FUN=function(package){
      if(package %in% rownames(installed.packages())){packageDescription(package)$Version} else ""
    })
    # --------
    output$package_tab <- DT::renderDataTable({
      # -------------
      datatable(data.table(
        Package = packages,
        Installed = status,
        Version = version
      ),
      selection = 'none',
      autoHideNavigation = T,
      options = list(lengthMenu = c(10, 20, 30), pageLength = 10)
      , rownames = F)
    })
  }
})

observeEvent(input$install_packages, {
  req(packages)
  # ----------
  p_load(char = packages, character.only = T)
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
})

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
  # --- restart ---
  stopApp()
  runApp(".")
})
  

if(packages_installed == "Y"){
  
# ====================== SETTINGS =================

observe({  
  shinyDirChoose(input, "get_db_dir", roots = getVolumes(), session = session)
  given_dir <- input$get_db_dir$path
  if(is.null(given_dir)) return()
  dbDir <<- paste(c(given_dir), collapse="/")
  output$curr_db_dir <- renderText(dbDir)
  # --- connect ---
  setOption(".conf", "db_dir", dbDir)
})

observe({  
  shinyDirChoose(input, "get_work_dir", roots = getVolumes(), session = session)
  given_dir <- input$get_work_dir$path
  if(is.null(given_dir)) return()
  exp_dir <<- paste(c(given_dir), collapse="/")
  output$exp_dir <- renderText(exp_dir)
  if(!dir.exists(exp_dir)) dir.create(exp_dir)
  # --- connect ---
  setOption(".conf", "work_dir", exp_dir)
})

observeEvent(input$set_proj_name, {
  proj_name <<- input$proj_name
  patdb <<- file.path(exp_dir, paste0(proj_name,".db", sep=""))
  output$proj_name <<- renderText(proj_name)
  # --- connect ---
  setOption(".conf", "proj_name", proj_name)
})

observeEvent(input$set_ppm, {
  ppm <<- input$ppm
  output$ppm <<- renderText(ppm)
  # --- connect ---
  setOption(".conf", "ppm", ppm)
})


color.pickers <- reactive({
  req(dataSet)
  # -------------
  if("facA" %not in% names(dataSet) & "facB" %not in% names(dataSet)){print('aah')}
  lbl.fac <- if(dataSet$facA.lbl == "Time") "facB" else "facA"
  default.colours <- rainbow(length(levels(dataSet[[lbl.fac]])))
  facs <- levels(dataSet[[lbl.fac]])
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
  if("facA" %not in% names(dataSet) & "facB" %not in% names(dataSet)) return(c("Red", "Green"))
  lbl.fac <- if(dataSet$facA.lbl == "Time") "facB" else "facA"
  default.colours <- rainbow(length(levels(dataSet[[lbl.fac]])))
  facs <- levels(dataSet[[lbl.fac]])
  # -------------
  unlist(lapply(seq_along(facs), function(i) {
    input[[paste("col", i, sep="_")]]
  }))
})

# ======================== DB CHECK ============================

output$umc_logo <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/umclogo.jpg'))
# Return a list containing the filename and alt text
  list(src = filename,
       alt = "UMC Utrecht",
       width=120,
       height=100)
  }, deleteFile = FALSE)

output$hmdb_logo <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/hmdblogo.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       alt = "The Human Metabolome DataBase",
       width = 150,
       height = 100)
}, deleteFile = FALSE)

output$chebi_logo <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/chebilogo.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       alt = "Chemical Entities of Biological Interest",
       width=100,
       height=100)
  
}, deleteFile = FALSE)

output$pubchem_logo <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/pubchemlogo.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       alt = "PubChem",
       width=140,
       height=100)
}, deleteFile = FALSE)

# --- check for db files ---

observeEvent(input$check_umc,{
  db_folder_files <- list.files(dbDir)
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
  db_folder_files <- list.files(dbDir)
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
  db_folder_files <- list.files(dbDir)
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
  db_folder_files <- list.files(dbDir)
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
    build.base.db("internal", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("internal", 
                      outfolder=dbDir,
                      adduct.table = wkz.adduct.confirmed, 
                      cl=session_cl, 
                      fetch.limit=10)
    build.base.db("noise", outfolder=dbDir) # does both because its a special db
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_hmdb,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("hmdb", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("hmdb", 
                      outfolder=dbDir, 
                      adduct.table = wkz.adduct.confirmed, 
                      cl=session_cl, 
                      fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_chebi,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("chebi", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("chebi", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_pubchem,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("pubchem", outfolder=dbDir, cl = session_cl)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("pubchem", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

# ================== DATA IMPORT ===========================

output$pos_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/handpos.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=120,
       height=120)
}, deleteFile = FALSE)

output$neg_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/handneg.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=120,
       height=120)
}, deleteFile = FALSE)

output$excel_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/excel.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=120,
       height=120)
}, deleteFile = FALSE)

output$db_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/office.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=100,
       height=100)
}, deleteFile = FALSE)

observeEvent(input$create_db,{
  # --------------------
  print(exp_dir)
  patdb <<- file.path(exp_dir, paste0(proj_name, ".db"))
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
  mode <<- if("Time" %not in% names(mz)[1:5]) "stat" else "time"
  # ----------------------------
  overview_tab <- if(mode == "time"){
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
  })

output$csv_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/office.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=100,
       height=100)
}, deleteFile = FALSE)

observeEvent(input$exp_type, {
 modes <- strsplit(input$exp_type, pattern = '//.')
 mainmode <<- modes[[1]]
 submode <<- modes[[2]]
 })
  
observeEvent(input$create_csv, {
  req(proj_name)
  req(exp_dir)
  # ---------
  withProgress({
    setProgress(1/4, "Creating csv file for MetaboAnalyst...")
    # create csv
    patdb
    mz = get.csv(patdb,
                 time.series = if(input$exp_type == "time") T else F,
                 exp.condition = input$exp_var)
    mz = get.csv(patdb,
                 time.series = T,
                 exp.condition = "diet")
    # --- experimental stuff ---
    yourTime <- 3
    submode="custom"
    switch(submode, 
           standard = { mz.adj <- mz }, 
           custom = { mz.adj <- unique(mz[Time == yourTime, -"Time"])},
           subtract = { 
             uniq.samples <- unique(mz$Sample)
             uniq.samples
             table.base <- mz[Time == yourTime,c(1,3)][on=uniq.samples]
             time.end <- mz[Time == yourTime,][,4:ncol(mz),on=uniq.samples]
             time.begin <- mz[Time == min(as.numeric(mz$Time)),][,4:ncol(mz),on=uniq.samples]
             mz.adj <- cbind(table.base, 
                              time.end - time.begin) 
           })
    unique(mz[Sample %in% uniq.samples, 1:3])
    View(unique(mz[,1:3]))
    View(mz[Time == yourTime,4:ncol(mz)] -mz [Time == 1,4:ncol(mz)])
    mz
    # -------------------------
    # save csv
    setProgress(2/4, "Writing csv file...")
    csv_loc <<- file.path(exp_dir, paste0(proj_name,".csv"))
    fwrite(mz, csv_loc, sep="\t")
    # --- overview table ---
    setProgress(3/4, "Creating overview table...")
    
    overview_tab <- if(input$exp_type == "time"){
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
  switch(main_mode,
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
           Read.TextData(csv_loc, "row", "disc")
           print(names(dataSet))
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
             iPCA.Anal(file.path(exp_dir, "ipca_3d_0_.json"))
             analSet[["ipca"]] <- "Done!"             
           }
           json_pca <- fromJSON(file.path(exp_dir, "ipca_3d_0_.json"))
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
            performMB(10, dir=exp_dir)
           }
           meba.table <<- read.csv(file.path(exp_dir, 'meba_sig_features.csv'))
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
              CalculateImpVarCutoff(0.05, 0.9, dir=exp_dir)
            }
           asca.table <<- read.csv(file.path(exp_dir,'Sig_features_Model_ab.csv'))
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
             req(analSet$pca)
             req(input$pca_x)
             req(input$pca_y)
             req(input$pca_z)
             fac.lvls <- unique(dataSet$facB.lbl)
             dataSet$facA.lbl
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
             
           })},
         heatmap = {NULL},
         tt = {Ttests.Anal(F, 0.05, FALSE, TRUE)},
         fc = {FC.Anal.unpaired(2.0, 1)})
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
}