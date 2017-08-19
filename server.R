function(input, output, session) {
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

output$db_icon <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/database.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       alt = "Your own db",
       width = 100,
       height = 100)
}, deleteFile = FALSE)

output$db_icon2 <- renderImage({
  # When input$n is 3, filename is ./images/image3.jpeg
  filename <- normalizePath(file.path('./backend/img/database.png'))
  # Return a list containing the filename and alt text
  list(src = filename,
       alt = "Your own db",
       width = 100,
       height = 100)
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
  withProgress({
    setProgress(message = "Working...")
    #build.base.db("internal", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("internal", 
                      outfolder=dbDir,
                      adduct.table = wkz.adduct.confirmed, 
                      cl=cl, 
                      fetch.limit=10)
    build.base.db("noise", outfolder=dbDir) # does both because its a special db
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_hmdb,{
  withProgress({
    setProgress(message = "Working...")
    #build.base.db("hmdb", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("hmdb", 
                      outfolder=dbDir, 
                      adduct.table = wkz.adduct.confirmed, 
                      cl=cl, 
                      fetch.limit=100)
    setProgress(message = "Ok!")
  })
})

observeEvent(input$build_chebi,{
  withProgress({
    setProgress(message = "Working...")
    build.base.db("chebi", outfolder=dbDir)
    setProgress(0.5,message = "Halfway there...")
    build.extended.db("chebi", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=cl, fetch.limit=100)
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

observeEvent(input$import_files,{
  withProgress({
    patdb <<- "./backend/appdata/experiment.db"
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

# ====================== MATCHING ==========================

observeEvent(input$match_start,{
  req(input$match_group)
  req(patdb)
  progress = 0.0
  withProgress({
    for(db in input$match_group){
      dbfile <- paste0(tolower(db), ".full.db")
      setProgress(progress ,message = fn$paste("Finding $db matches..."))
      iden.code.binned(patdb, file.path("./backend/db", dbfile), isofilt=T)
      progress <- progress + 1/length(input$match_group)
    }
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
  # requirements
  req(patdb)
  # create csv
mz = get.csv(patdb,
             time.series = T,
             group.by.adduct = F,
             exp.condition = input$exp_var,
             chosen.display = "mz")
# save csv
  metabo_csv <<- "./backend/appdata/mz.csv"
  #write.table(mz, metabo_csv, sep="\t", row.names=F)
  # get curr values from: input$ exp_type, filt_type, norm_type, scale_type, trans_type (later)
  sourceAll(file.path("backend", 
                      "scripts", 
                      "metaboanalyst"))
  #Below is your R command history: 
  InitDataObjects("pktable", "ts", FALSE)
  SetDesignType("time")
  Read.TextData(metabo_csv, "rowts", "disc")
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
  Normalization("QuantileNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
  output$var_norm_plot <- renderPlot(PlotNormSummary())
  output$samp_norm_plot <- renderPlot(PlotSampleNormSummary())
})

# ======================== ASCA ============================

output$asca_tab <- DT::renderDataTable({
  datatable(asca.table, 
            selection = 'single',
            colnames = c("Mass/charge", "Leverage", "SPE"),
            autoHideNavigation = T,
            options = list(lengthMenu = c(10, 30, 50), pageLength = 10))
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
  match_list <- lapply(input$checkGroup, FUN=function(match.table){
    get_matches(curr_mz, match.table)})
  result_table <- as.data.table(rbindlist(match_list))
  output$match_tab <- DT::renderDataTable({
    datatable(unique(result_table[Compound != "",]),
              selection = 'single',
              autoHideNavigation = T,
              options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })
})

# --- ON CLOSE ---
session$onSessionEnded(function() {
  stopCluster(cl)
})
}