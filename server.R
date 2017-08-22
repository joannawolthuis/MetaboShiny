function(input, output, session) {
  
# ===== defaults =====
output$exp_dir <- renderText(exp_dir)
output$proj_name <- renderText(proj_name)
  
# ====================== SETTINGS =================

volumes = getVolumes()
observe({  
  shinyDirChoose(input, "get_work_dir", roots = volumes, session = session)
  given_dir <- input$get_work_dir$path
  if(is.null(given_dir)) return()
  exp_dir <<- paste(c(given_dir), collapse="/")
  output$exp_dir <- renderText(exp_dir)
  if(!dir.exists(exp_dir)) dir.create(exp_dir)
  options_raw[[1]] <- paste("work_dir = '",exp_dir,"'", sep="")
  # --- connect ---
  opt_conn <- file(".conf")
  writeLines(opt_conn, text = options_raw)
  close(opt_conn)
})

observeEvent(input$set_proj_name, {
  patdb <<- file.path(exp_dir, paste0(input$proj_name,".db", sep=""))
  print(patdb)
  output$proj_name <- renderText(proj_name)
  options_raw[[2]] <- paste("proj_name = '", exp_dir,"'", sep="")
  # --- connect ---
  opt_conn <- file(".conf")
  writeLines(opt_conn)
  close(opt_conn)
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

# --- check for db files ---

db_folder_files <- list.files("./backend/db")

observeEvent(input$check_umc,{
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

# --- build db ---

observeEvent(input$build_umc,{
  session_cl <<- if(is.na(session_cl)) makeCluster(4, type="FORK") else session_cl
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
  session_cl <<- if(is.na(session_cl)) makeCluster(4, type="FORK") else session_cl
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
  session_cl <<- if(is.na(session_cl)) makeCluster(4, type="FORK") else session_cl
  withProgress({
    setProgress(message = "Working...")
    build.base.db("chebi", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("chebi", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
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
  print(head(mz[,1:5]))
  overview_tab <- t(data.table(
    Mz = ncol(mz) - 3,
    Samples = nrow(mz),
    Groups = length(unique(mz$Label)),
    Timepoints = length(unique(mz$Time))
  ))
  output$csv_tab <- DT::renderDataTable({
    datatable(overview_tab, 
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE))
  })  })

output$csv_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/office.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       width=100,
       height=100)
}, deleteFile = FALSE)

observeEvent(input$create_csv, {
  req(proj_name)
  req(exp_dir)
  # ---------
  withProgress({
    setProgress(1/5, "Finding and filtering matches for mz values...") # move somewhere else later??
    total.matches <- 0
    for(db in c("internal", "noise", "hmdb", "chebi")){
      dbfile <- paste0(tolower(db), ".full.db")
      match.count <- iden.code.binned(patdb, file.path("./backend/db", dbfile), isofilt=T)
      print(match.count)
      total.matches <- total.matches + match.count
    }
    setProgress(2/5, "Creating csv file for MetaboAnalyst...")
    # create csv
    mz = get.csv(patdb,
                 time.series = T,
                 group.by.adduct = F,
                 exp.condition = input$exp_var,
                 chosen.display = "mz")
    # save csv
    setProgress(3/5, "Writing csv file...")
    csv_loc <<- file.path(exp_dir, paste0(proj_name,".csv"))
    fwrite(mz, csv_loc, sep="\t")
    # --- overview table ---
    setProgress(4/5, "Creating overview table...")
    overview_tab <- t(data.table(keep.rownames = F,
      Mz = ncol(mz) - 3,
      Samples = nrow(mz),
      Timepoints = length(unique(mz$Time)),
      Groups = length(unique(mz$Label)),
      Matches = total.matches
    ))
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

observeEvent(input$initialize,{
  # match
  withProgress({
  # get curr values from: input$ exp_type, filt_type, norm_type, scale_type, trans_type (later)
  sourceAll(file.path("backend", 
                      "scripts", 
                      "metaboanalyst"))
  setProgress(.1, "Applying your settings...")
  #Below is your R command history: 
  InitDataObjects("pktable", "ts", FALSE)
  SetDesignType("time")
  Read.TextData(csv_loc, "rowts", "disc")
  SanityCheckData()
  RemoveMissingPercent(percent=0.5)
  ImputeVar(method="min")
  ReplaceMin()
  FilterVariable("iqr", "F", 25)
  # # ---- here facA/facB disappears?? ---
  GetPrenormSmplNms()
  GetPrenormFeatureNms()
  GetPrenormClsNms()
  UpdateGroupItems()
  UpdateSampleItems()
  UpdateFeatureItems()
  # # ------------------------------------
  setProgress(.2, "Normalizing data...")
  Normalization("QuantileNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
  setProgress(.3, "Plotting variable overview...")
  output$var_norm_plot <- renderPlot(PlotNormSummary())
  setProgress(.4, "Plotting sample overview...")
  output$samp_norm_plot <- renderPlot(PlotSampleNormSummary())
  # ------ start matching analyses ---------
  setProgress(.5, "Starting analyses...")
  if(mode == "time"){
    setProgress(.6, "Running iPCA...")
    iPCA.Anal(file.path(exp_dir, "ipca_3d_0_.json"))
    # --- asca stuff ---
    setProgress(.7, "Running ASCA")
    Perform.ASCA(1, 1, 2, 2)
    CalculateImpVarCutoff(0.05, 0.9, dir=exp_dir)
    asca.table <<- read.csv(file.path(exp_dir,'Sig_features_Model_ab.csv'))
    # --- do anova only if balanced ---
    NULL
    # --- meba stuff ---
    setProgress(.8, "Running MEBA...")
    performMB(10, dir=exp_dir)
    meba.table <<- read.csv(file.path(exp_dir, 'meba_sig_features.csv'))
  } else{NULL}
  json_pca <<- fromJSON(file.path(exp_dir, "ipca_3d_0_.json")) # but perform first...
  updateSelectInput(session, "ipca_factor",
                    choices = grep(names(json_pca$score), pattern = "^fac[A-Z]", value = T))
  setProgress(.9, "Saving results to R-loadable file...")
  save(dataSet, analSet, file=file.path(exp_dir, paste0(proj_name, "analysed.RData")))
  })
})

# ======================== IPCA ==========================

  output$plot_ipca <- renderPlotly({
    req(json_pca)
    # ---------------
    df <- t(as.data.frame(json_pca$score$xyz))
    plot_ly() %>%
        add_trace(
          x = df[,1], 
          y = df[,2], 
          z = df[,3], 
          type = "scatter3d",
          color= json_pca$score[[input$ipca_factor]], colors=c("pink", "skyblue")
        ) %>%  layout(scene = list(xaxis = list(
          title = json_pca$score$axis[1]),
          yaxis = list(
            title = json_pca$score$axis[2]),
          zaxis = list(
            title = json_pca$score$axis[3])))
    })
  
# =================== HEATMAP =====================

#output$time_heat_plot <- renderPlot(
  #PlotHeatMap2("euclidean","ward.D","bwm","overview", F, 1, F, F)
#)

# =================== MEBA ========================

output$meba_tab <- DT::renderDataTable({
  req(meba.table)
  # -------------
  datatable(meba.table, 
            selection = 'single',
            colnames = c("Mass/charge", "Hotelling/T2 score"),
            autoHideNavigation = T,
            options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
})

observeEvent(input$meba_tab_rows_selected,{
  curr_row = input$meba_tab_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- meba.table[curr_row,'X']
  output$meba_plot <- renderPlot(PlotMBTimeProfile(curr_mz)
)
})

# ======================== ASCA ============================

output$asca_tab <- DT::renderDataTable({
  req(asca.table)
  # -------------
  datatable(asca.table, 
            selection = 'single',
            colnames = c("Mass/charge", "Leverage", "SPE"),
            autoHideNavigation = T,
            options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
})

# check for selected mz row
observeEvent(input$asca_tab_rows_selected,{
  curr_row = input$asca_tab_rows_selected
  # do nothing if not clicked yet, or the clicked cell is not in the 1st column
  if (is.null(curr_row)) return()
  curr_mz <<- asca.table[curr_row,'X']
  output$asca_plot <- renderPlot(PlotCmpdSummary(curr_mz))
})

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
  result_table <<- unique(as.data.table(rbindlist(match_list))[Compound != ""])
  output$match_tab <- DT::renderDataTable({
    datatable(result_table[,-"Description"],
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 10, 15), pageLength = 5))
  })
})

observeEvent(input$match_tab_rows_selected,{
  curr_row = input$match_tab_rows_selected
  if (is.null(curr_row)) return()
  # -----------------------------
  curr_def <<- result_table[curr_row,'Description']
  output$curr_definition <- renderText(curr_def$Description)
  })

# --- ON CLOSE ---
session$onSessionEnded(function() {
  if(!is.na(session_cl)) stopCluster(session_cl)
})
}