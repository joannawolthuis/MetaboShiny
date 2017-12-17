shinyServer(function(input, output, session) {
  
  # dir
  shinyDirChoose(input, 'dir', roots = c(home = '~'), filetypes = c('', 'txt'))
  dir <- reactive(input$dir)
  output$dir <- renderPrint(dir())
  
  observe({
    shinyFileChoose(input, id = 'sample_sheet_loc', roots = c(home = '~'), filetypes = c('csv','xls', 'xlt', 'xlsx'))
    sample_sheet <- reactive(input$sample_sheet_loc)
    if(!is.null(sample_sheet())){
      sheetpath <<- as.character(parseFilePaths( c(home = '~'), sample_sheet())$datapath)
      sample_tab <<- as.data.table(read.xlsx(file = sheetpath, sheetIndex = 1))[,1:2]
      print(head(sample_tab))
      colnames(sample_tab) <<- c("File.Name", "Card.ID")
    }
    # -----------------------------------------------------------------
  })

  observe({
    shinyFileSave(input, id='report_save_loc', roots = c(home = '~'), filetype=c('csv'))
    report_loc <- reactive(input$report_save_loc)
    adj_loc <- report_loc()
    adj_loc$type <- "csv"
    if(!is.null(report_loc()) & exists("final_report") ){ 
      # -------------------
      savepath <- as.character(parseSavePath(roots =  c(home = '~'), selection = adj_loc )$datapath)
      print(savepath)
      write.csv(x = final_report, file = savepath)
    }
    # -----------------------------------------------------------------
  })
  # check
  
  observe({
    shinyFileSave(input, id='ticreport_save_loc', roots = c(home = '~'), filetype=c('csv'))
    tic_report_loc <- reactive(input$ticreport_save_loc)
    adj_loc <- tic_report_loc()
    adj_loc$type <- "csv"
    if(!is.null(tic_report_loc()) & 
       exists("ticSelection") ){ 
      # -------------------
      savepath <- as.character(parseSavePath(roots =  c(home = '~'), selection = adj_loc )$datapath)
      print(savepath)
      write.csv(x = ticSelection, file = savepath)
    }
    # -----------------------------------------------------------------
  })
  
  observeEvent(input$convert, {
    "
Requires the following:\
WINE - 32 BITS INSTALL (please adjust your path below to the correct executable)
     - PROPER WINEPREFIX (see below for example)
MSFileReader - Thermo program, requires registration. Installs through wine.
             - May need to 'install' the DLL's manually through wine.
READW - Nice converter that requires MSFileReader. Find a one-file executable please!
--- This may be really tricky to get wo work, I feel lucky it works currently lol :I ---
    "
    minFiles = 0
    progr = 0
    files = list.files(path(),pattern = "raw",full.names = T)
    done = list.files(path(),pattern = "mzXML",full.names = T)
    maxFiles <- length(files) - length(done)
    READW <- file.path(getwd(),"exec/ReAdW.2016010.msfilereader.exe")
    WINE <- "/Users/jwolthuis/brew/Cellar/wine/2.19_3/bin/wine"
    # ----------------------
    withProgress(min = 0, max = maxFiles, {
      for(INFILE in files){
        OUTFILE = gsub(INFILE, pattern = "\\.raw$|\\.RAW$", replacement = "\\.mzXML")
        if(file.exists(OUTFILE)) next
        system(fn$paste(
          "env WINEPREFIX=/Users/jwolthuis/.wine32 $WINE '$READW' '$INFILE' '$OUTFILE'"
          ))
        progr <- progr + 1
        setProgress(message = paste(basename(INFILE), "converted", sep=": "),value = progr)
      }
    })
  })
  
  
  observeEvent(input$check, {
    minFiles = 0
    progr = 0
    files = list.files(path(),pattern = "mzXML",full.names = T)
    maxFiles <- length(files)
    # ----------------------
    withProgress(min = 0, max = maxFiles, {
      for(f in files){
        result <- CheckFile(f)
        progr <- progr + 1
        setProgress(message = paste(basename(f), result, sep=": "),value = progr)
      }
    })
  })
  

  observeEvent(input$plot, {
    req(projname)
    # -----------
    ticDir <<- file.path(getwd(), "www", projname)
    if(!dir.exists(ticDir)) dir.create(ticDir)
    minFiles = 0
    progr = 0
    xmlfiles <<- list.files(path(),pattern = "mzXML",full.names = T)
    donefiles = list.files(ticDir, pattern="png")
    maxFiles <- length(xmlfiles) - length(donefiles)
    if(maxFiles != 0){
      withProgress(min = 0, max = maxFiles, {
        for(f in xmlfiles){
          imgName <- paste(ticDir, "/", gsub(basename(f), pattern = "\\.mzXML", replacement = "\\.png"), sep="")
          if(file.exists(imgName)){print("already plotted");next;}
          # ----------
          rawF <- NULL
          attempt <- 1
          while( is.null(rawF) && attempt <= 3 ) {
            attempt <- attempt + 1
            try(
              rawF <- xcmsRaw(f, profstep=0.1)
            )
          }
          # -----------
          png(filename=imgName, 320, 240)
          plotTIC(rawF, ident=FALSE, msident=FALSE)
          dev.off()
          progr <- progr + 1
          setProgress(message = paste("plotting",basename(f)),value = progr)
        }
      })
    }
    # ----------------------
    tics <<- list.files(ticDir, pattern = "png")
    #files <<- list.files(path(),pattern = "mzXML|raw")
    output$tics <- renderUI({
      req(projname)
      req(ticDir)
      # ---------------------------
      lapply(1:length(tics), FUN=function(x){
        tic <- tics[x]
        sardine(fadeImageButton(paste0("tic_", x), img.path = file.path(projname, tic),value = T))
      })
    })
    ticTable <<- data.table(filenames = gsub(tics, pattern = "\\.png", replacement = ".mzXML"),
                           judgement = c("OK"))
    output$ticTable <- DT::renderDataTable({
      datatable(ticTable,
      autoHideNavigation = T,
      options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE)
      ) %>% formatStyle(
        'judgement',
        target = 'row',
        backgroundColor = styleEqual(c('OK', "BAD"), c('lightgreen', "lightred"))
      )
    })
  })
  

  observeEvent(input$submit, {
    nTics <- length(tics)
    judgement <<- sapply(1:nTics, FUN=function(i){
      ticJudge <- !input[[paste0("tic_", i)]]
      # --- return ---
      if(ticJudge == TRUE){"OK"} else("BAD")
    })
    ticSelection <<- data.table(filenames = gsub(tics, pattern = "\\.png", replacement = ".mzXML"),
                                judgement = judgement)
    output$ticTable <- DT::renderDataTable({
      datatable(ticSelection,
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE)
      ) %>% formatStyle(
        'judgement',
        target = 'row',
        backgroundColor = styleEqual(c('OK', 'BAD'), c('lightgreen', 'red'))
      )
    })
  })
  

  observeEvent(input$make_report, {
    # create final report.
    # ------------------
    rawfiles <- data.table(RAW=list.files(path(),pattern = "raw"))
    rawfiles <- data.table(filenames=gsub(pattern = "\\.raw", 
                                replacement = "",
                                x = list.files("/Users/jwolthuis/PROCESSING_HPC/",
                                               pattern = "raw")))
    ticSelection_subbed <- ticSelection
    ticSelection_subbed$filenames <- gsub(ticSelection_subbed$filenames,
                                   pattern = "\\.mzXML", 
                                   replacement = "")
    
    setkey(rawfiles,filenames)
    setkey(ticSelection_subbed,filenames)
    
    res_1 <- merge(rawfiles, ticSelection_subbed, all=TRUE)
    res_2 <- merge(res_1, sample_tab, all=TRUE, by.x = "filenames", by.y='File.Name')    
    # -----------
    library(dplyr)
    set.seed(1)
    res_3 <- as.data.table(
      res_2 %>% 
      mutate(OK=(judgement=="OK")) %>% 
      group_by(Card.ID) %>%
      summarise(OK=sum(OK))
    )
    # ===========
    warning_samples <- res_3[OK < 2]
    warning_samples$result <- c("BAD")
    ok_samples <- res_3[OK >= 2]
    ok_samples$result <- c("OK")
    # -----------
    final_report <<- rbind(ok_samples, warning_samples)
    
    output$final_report <- DT::renderDataTable({
      datatable(final_report,
                autoHideNavigation = T,
                options = list(lengthMenu = c(10, 30, 50), pageLength = 30,scrollX=TRUE, scrollY=TRUE)
      ) %>% formatStyle(
        'result',
        target = 'row',
        backgroundColor = styleEqual(c('OK', 'BAD'), c('lightgreen', 'red'))
      )
    })
    })
  
  observeEvent(input$save_report,{
    req(input$report_title)
    # ---------------------
    info <- final_report
    textplot( info, valign="top",cex = 1)
    title(input$report_title)
  })
  
  # path
  path <- reactive({
    home <- normalizePath("~")
    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })
  
  # files
  observe({
    files <<- list.files(path(),pattern = "mzXML|raw")
    projname <<- gsub(basename(files[1]),pattern = ".*\\/|_\\d*.raw|_\\d*.mzXML", replacement = "")
    })
  
  output$files <- renderPrint(list.files(path(),pattern = "mzXML|raw"))
}) 

makeSampleNames <- function(dir, injfile){ # IMPLEMENT THIS W/ PEAKFINDING SOON
  library(data.table)
  
  dir <- "/Users/jwolthuis/Downloads/Data/Project 2017-025 DSM feed-9 (SRUC-Ayr) - Saskia v Mil/RES-2017-11-14_DSM DBS Ayr Scotland/"
  mzfiles <- list.files(file.path(dir, "MZXML"))
  injfile <- xlsx::read.xlsx(file = "/Users/jwolthuis/Downloads/Data/Project 2017-025 DSM feed-9 (SRUC-Ayr) - Saskia v Mil/Injection list Ayr (Scotland).xlsx",
                             sheetIndex = 1,colIndex = c(1:3),rowIndex = 1:length(mzfiles)/3 + 1)
  df <- as.data.frame(injfile)
  
  filepfx <- gsub(mzfiles[[1]], 
                  pattern = "_\\d*\\.mzXML", 
                  replacement = "" )
  
  df.expanded <- df[rep(1:nrow(df),each=3),] # replicados
  samplecount = nrow(df)
  filenames <- paste(filepfx, c(1:nrow(df.expanded)), sep = "_")
  filenames <- gsub(filenames,
                    pattern="(_)([0-9])$",
                    replacement="\\100\\2")
  filenames <- gsub(filenames,
                    pattern="(_)([0-9][0-9])$",
                    replacement="\\10\\2")
  df.w.filenames <- df.expanded
  df.w.filenames$filename <- filenames
  dt <- as.data.table(df.w.filenames)
  
  sampleNames <- dt[,c(3,2)]
  colnames(sampleNames) <- c("File_Name", "Sample_Name")
  
  fwrite(sampleNames,file = file.path(dir, "sampleNames_UK.txt"),sep = "\t")
}
