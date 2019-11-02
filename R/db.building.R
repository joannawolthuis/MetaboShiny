#' @export
build.pat.db <- function(db.name,
                         pospath,
                         negpath,
                         overwrite=FALSE,
                         rtree=TRUE,
                         make.full = TRUE,
                         ppm=2,
                         inshiny=F){

  ppm = as.numeric(ppm)

  
  poslist <- data.table::fread(pospath,header = T)
  neglist <- data.table::fread(negpath,header = T) 
  
  keepcols <- intersect(colnames(poslist), colnames(neglist))
  
  poslist <- poslist[,..keepcols]
  neglist <- neglist[,..keepcols]
  
  # replace commas with dots
  poslist <- as.data.table(sapply(poslist, gsub, pattern = ",", replacement= "."))
  neglist <- as.data.table(sapply(neglist, gsub, pattern = ",", replacement= "."))

  if(inshiny) setProgress(.20)
  
  gc()
  
  mzvals <- data.table::data.table(mzmed = c(as.numeric(poslist$mzmed), as.numeric(neglist$mzmed)),
                                   foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  
  mzvals$foundinmode <- trimws(mzvals$foundinmode)
  
  gc()
  
  if(inshiny) setProgress(.30)
  
  poslist <- data.table::melt(poslist,#[,(rmv.cols) := NULL],
                              id.vars="mzmed",
                              variable.name="filename",
                              value.name="intensity",
                              variable.factor=TRUE
  )
  neglist <- data.table::melt(neglist,#[,(rmv.cols) := NULL],
                              id.vars="mzmed",
                              variable.name="filename",
                              value.name="intensity",
                              variable.factor=TRUE
  )
  if(inshiny) setProgress(.40)
  

  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)
  
  sql.make.int <- strwrap("CREATE TABLE mzintensities(
                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
                          mzmed decimal(30,13),
                          filename text,
                          intensity float)", width=10000, simplify=TRUE)
  
  RSQLite::dbExecute(conn, sql.make.int)
  
  if(inshiny) setProgress(.60)
  
  # --- write intensities to table and index ---
  RSQLite::dbWriteTable(conn, "mzintensities", poslist, append=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "mzintensities", neglist, append=TRUE) # insert into
  
  RSQLite::dbExecute(conn, "CREATE INDEX intindex ON mzintensities(filename,'mzmed',intensity)")
  RSQLite::dbExecute(conn, "CREATE INDEX intindex2 ON mzintensities('mzmed')")
  
  # ------------------------
  
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, sql.make.meta)
  RSQLite::dbExecute(conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  
  if(inshiny) setProgress(.70)
  
  RSQLite::dbWriteTable(conn, "params", 
                        data.table::data.table(ppm=ppm), 
                        overwrite=T)
  
  # --- write vals to table ---
  RSQLite::dbWriteTable(conn, "mzvals", mzvals, append=TRUE) # insert into
  
  # --- make range table (choose if R*tree or not) ---
 
  if(inshiny) setProgress(.80)
  
  # --- cleanup ---
  RSQLite::dbExecute(conn, "VACUUM")
  
  if(inshiny) setProgress(.90)
  
  # ----------------
  RSQLite::dbDisconnect(conn)
}


#' @export
load.metadata.csv <- function(path.to.csv,
                              path.to.patdb){
  
  # --- connect to sqlite db ---
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  filenames = RSQLite::dbGetQuery(conn, "SELECT DISTINCT filename FROM mzintensities")
  
  tab <- data.table::fread(path.to.csv)
  colnames(tab) <- tolower(colnames(tab))
  
  if(any(colnames(tab)=="sampling_date")){
    if(any(is.na(as.numeric(tab$sampling_date)))){
      tab$sampling_date <- as.factor(as.Date(as.character(tab$sampling_date),
                                                         format = "%d-%m-%y"))
    }else{
      tab$sampling_date <- as.factor(as.Date(as.numeric(tab$sampling_date),
                                                         origin = "1899-12-30"))
    }
    if(is.na(tab$sampling_date[1])) levels(tab$sampling_date) <- factor(1)
  }
  
  if(length(intersect(filenames[,1], tab$sample)) == 0){
    stop("Complete mismatch between metadata sample names and file sample names. Aborting...")
  }
  
  colnames(tab) <- tolower(gsub(x=colnames(tab), pattern = "\\.$|\\.\\.$", replacement = ""))
  colnames(tab) <- gsub(x=colnames(tab), pattern = "\\.|\\.\\.| ", replacement = "_")
  
  RSQLite::dbWriteTable(conn, "individual_data", tab, overwrite=TRUE) # insert into
  
  RSQLite::dbDisconnect(conn)
}

#' @export
load.metadata.excel <- function(path.to.xlsx,
                                path.to.patdb,
                                tabs.to.read = c(
                                  "Setup",
                                  "Individual Data"
                                )){
  
  # --- connect to sqlite db ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  # -------------------------------
  data.store <- pbapply::pblapply(tabs.to.read, FUN=function(tab.name){
    #tab <- data.table::as.data.table(xlsx::read.xlsx(path.to.xlsx, sheetName = tab.name))
    tab <- data.table::as.data.table(openxlsx::read.xlsx(path.to.xlsx, sheet = tab.name))
    # --- reformat colnames ---
    colnames(tab) <- tolower(gsub(x=colnames(tab), pattern = "\\.$|\\.\\.$", replacement = ""))
    colnames(tab) <- gsub(x=colnames(tab), pattern = "\\.|\\.\\.", replacement = "_")
    #colnames(tab)[grep(x=colnames(tab), pattern= "*date*")] <- "sampling_date"
    # -----------------------------------------------------------------------------------------
    data.table::as.data.table(tab, keep.rownames=F)
  })
  
  # --- convert to data table --- ## make this nicer loooking in the future
  setup <- data.store[[1]]
  individual.data <- data.store[[2]]
  
  # --- fill empty cells w/ na ---
  
  indx <- which(sapply(setup, is.character))
  for (j in indx) set(setup, i = grep("^$|^ $", setup[[j]]), j = j, value = NA_character_)
  
  indx <- which(sapply(individual.data, is.character))
  for (j in indx) set(individual.data, i = grep("^$|^ $", individual.data[[j]]), j = j, value = NA_character_)
  
  # --- remove empty lines ---
  
  setup <- setup[rowSums(is.na(setup)) != ncol(setup),]
  individual.data <- individual.data[rowSums(is.na(individual.data)) != ncol(individual.data),]
  
  # --------------------------
  
  if(any(colnames(individual.data)=="sampling_date")){
    if(any(is.na(as.numeric(individual.data$sampling_date)))){
      individual.data$sampling_date <- as.factor(as.Date(as.character(individual.data$sampling_date),
                                                         format = "%d-%m-%y"))
    }else{
      individual.data$sampling_date <- as.factor(as.Date(as.numeric(individual.data$sampling_date),
                                                         origin = "1899-12-30"))
    }
    if(is.na(individual.data$sampling_date[1])) levels(individual.data$sampling_date) <- factor(1)
  }
  
  individual.data$sample <- as.character(individual.data$sample)
  individual.data$individual <- as.character(individual.data$individual)
 
  setup <- data.table::as.data.table(apply(setup, MARGIN=2, trimws))
  individual.data <- data.table::as.data.table(apply(individual.data, MARGIN=2, trimws))

  # --- import to patient sql file ---
  RSQLite::dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "individual_data", individual.data, overwrite=TRUE) # insert into
 # --- disconnect ---
  RSQLite::dbDisconnect(conn)
}

