shinyServer(function(input, output, session) {
  
  # dir
  shinyDirChoose(input, 'dir', roots = c(home = '~'), filetypes = c('', 'txt'))
  dir <- reactive(input$dir)
  output$dir <- renderPrint(dir())
  
  # check
  
  
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
    maxFiles <- length(files)
    READW <- file.path(getwd(),"exec/ReAdW.2016010.msfilereader.exe")
    WINE <- "/Users/jwolthuis/brew/Cellar/wine/2.19_3/bin/wine"
    # ----------------------
    withProgress(min = 0, max = maxFiles, {
      for(INFILE in files){
        OUTFILE = gsub(INFILE, pattern = "\\.raw$|\\.RAW$", replacement = "\\.mzXML")
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
    ticDir <<- file.path(getwd(), "www", "TICs")
    if(!dir.exists(ticDir)) dir.create(ticDir)
    print(ticDir)
    minFiles = 0
    progr = 0
    files = list.files(path(),pattern = "mzXML",full.names = T)
    maxFiles <- length(files)
    # ----------------------
    withProgress(min = 0, max = maxFiles, {
      for(f in files){
        imgName <- paste(ticDir, "/", basename(f), ".png", sep="")
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
  })
  
  ticDir <<- file.path(getwd(), "www", "TICs")
  ticsAll <- list.files(ticDir, pattern = "png")

  output$tics <- renderUI({
    tics <<- ticsAll
    # ---------------------------
     lapply(1:length(tics), FUN=function(x){
      tic <- tics[x]
      sardine(fadeImageButton(paste0("tic_", x), img.path = file.path("TICs",tic),value = T))
    })
  })

  
  # path
  path <- reactive({
    home <- normalizePath("~")
    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })
  
  # files
  output$files <- renderPrint(list.files(path(),pattern = "mzXML|raw"))
  
}) 

