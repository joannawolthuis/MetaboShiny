#' @export
build.pat.db <- function(db.name,
                         pospath,
                         negpath,
                         overwrite=FALSE,
                         rtree=TRUE,
                         make.full = TRUE,
                         ppm=2,
                         inshiny=F){

  # pospath="/Users/jwolthuis/Documents/Documents/xls/example_pos.csv"
  # negpath="/Users/jwolthuis/Documents/Documents/xls/example_neg.csv"
  
  ppm = as.numeric(ppm)

  
  poslist <- data.table::fread(pospath,header = T)
  neglist <- data.table::fread(negpath,header = T) 
  
  # ------------------------------------
  
  keepcols <- intersect(colnames(poslist), colnames(neglist))
  
  poslist <- poslist[,..keepcols]
  neglist <- neglist[,..keepcols]
  
  # replace commas with dots
  poslist <- as.data.table(sapply(poslist, gsub, pattern = ",", replacement= "."))
  neglist <- as.data.table(sapply(neglist, gsub, pattern = ",", replacement= "."))
  
  # - - - fix QCs - - -
  
  which.qc <- grep(colnames(poslist), pattern = "^QC")
  qc.i = 1
  
  for(qc in which.qc){
    new.qc.name <- paste0("QC", qc.i)
    new.qc.name <- gsub(colnames(poslist)[qc], pattern = "(^QC[\\d|\\d\\d])", replacement = new.qc.name,perl = T)
    colnames(poslist)[qc] <- new.qc.name
    colnames(neglist)[qc] <- new.qc.name
    # - - -
    qc.i = qc.i + 1
  }
  # - - - - - - - - - -
  if(inshiny) setProgress(.20)
  
  gc()
  
  mzvals <- data.table::data.table(mzmed = c(as.numeric(poslist$mzmed), as.numeric(neglist$mzmed)),
                                   foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  
  mzvals$foundinmode <- trimws(mzvals$foundinmode)
  
  # --- SAVE BATCH INFO (kinda ugly...  ; _;") ---
  
  if(any(grepl("\\*", x = colnames(poslist)))){
    samp_split = strsplit(colnames(poslist)[2:ncol(poslist)], "\\*")
    batch_split = strsplit(unlist(lapply(samp_split, function(x) x[2])), "\\_")
    batch_info = data.table::data.table(sample = sapply(samp_split, function(x) x[1]),
                                        batch = sapply(samp_split, function(x) x[2]),
                                        injection = sapply(samp_split, function(x) x[3]))
    colnames(poslist) = gsub(colnames(poslist), pattern = "(\\*.*$)", replacement = "")
    colnames(neglist) = gsub(colnames(poslist), pattern = "(\\*.*$)", replacement = "")
  } else{
    batch_info = NULL
  }
  
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
  
  # --- reconnect / remake ---
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)
  
  # ------------------------s
  
  if(!is.null(batch_info)){
    RSQLite::dbWriteTable(conn, "batchinfo", batch_info, overwrite=T) # insert into
  }
  
  # ------------------------
  
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
  
  #path.to.csv = "~/Downloads/maria_meta.csv"
  #path.to.patdb = "~/Downloads/maria_3ppm.db"
  
  
  # --- connect to sqlite db ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  
  tab <- data.table::fread(path.to.csv)
  colnames(tab) <- tolower(colnames(tab))
  
  colnames(tab) <- tolower(gsub(x=colnames(tab), pattern = "\\.$|\\.\\.$", replacement = ""))
  colnames(tab) <- gsub(x=colnames(tab), pattern = "\\.|\\.\\.| ", replacement = "_")
  colnames(tab)[grep(x=colnames(tab), pattern= "*date*")] <- "sampling_date"
  
  setup <- data.table(group = as.character(unique(tab[,c("group")][[1]])))
  
  RSQLite::dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "individual_data", tab, overwrite=TRUE) # insert into
  
  RSQLite::dbDisconnect(conn)
}

#' @export
load.metadata.excel <- function(path.to.xlsx,
                                path.to.patdb,
                                tabs.to.read = c(
                                  #"General",
                                  "Setup",
                                  "Individual Data"
                                  #,"Pen Data",
                                  #"Admin"
                                ),ppm){
  
  # --- connect to sqlite db ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  # -------------------------------
  data.store <- pbapply::pblapply(tabs.to.read, FUN=function(tab.name){
    #tab <- data.table::as.data.table(xlsx::read.xlsx(path.to.xlsx, sheetName = tab.name))
    tab <- data.table::as.data.table(openxlsx::read.xlsx(path.to.xlsx, sheet = tab.name))
    # --- reformat colnames ---
    colnames(tab) <- tolower(gsub(x=colnames(tab), pattern = "\\.$|\\.\\.$", replacement = ""))
    colnames(tab) <- gsub(x=colnames(tab), pattern = "\\.|\\.\\.", replacement = "_")
    colnames(tab)[grep(x=colnames(tab), pattern= "*date*")] <- "sampling_date"
    # -----------------------------------------------------------------------------------------
    data.table::as.data.table(tab, keep.rownames=F)
  })
  
  # --- convert to data table --- ## make this nicer loooking in the future
  #general <- data.store[[1]]
  setup <- data.store[[1]]
  individual.data <- data.store[[2]]
  
  # --- fill empty cells w/ na ---
  
  indx <- which(sapply(setup, is.character))
  for (j in indx) set(setup, i = grep("^$|^ $", setup[[j]]), j = j, value = NA_character_)
  
  # indx <- which(sapply(general, is.character))
  # for (j in indx) set(general, i = grep("^$|^ $", general[[j]]), j = j, value = NA_character_)
  
  indx <- which(sapply(individual.data, is.character))
  for (j in indx) set(individual.data, i = grep("^$|^ $", individual.data[[j]]), j = j, value = NA_character_)
  
  # --- remove empty lines ---
  
  #general <- general[rowSums(is.na(general)) != ncol(general),]
  setup <- setup[rowSums(is.na(setup)) != ncol(setup),]
  individual.data <- individual.data[rowSums(is.na(individual.data)) != ncol(individual.data),]
  
  # --------------------------
  
  if(any(is.na(as.numeric(individual.data$sampling_date)))){
    individual.data$sampling_date <- as.factor(as.Date(as.character(individual.data$sampling_date),
                                                       format = "%d-%m-%y"))
  }else{
    individual.data$sampling_date <- as.factor(as.Date(as.numeric(individual.data$sampling_date),
                                                       origin = "1899-12-30"))
    
  }
  
  individual.data$card_id <- as.character(individual.data$card_id)
  individual.data$animal_internal_id <- as.character(individual.data$animal_internal_id)
  
  if(is.na(individual.data$sampling_date[1])) levels(individual.data$sampling_date) <- factor(1)
  
  setup <- data.table::as.data.table(apply(setup, MARGIN=2, trimws))
  individual.data <- data.table::as.data.table(apply(individual.data, MARGIN=2, trimws))
  #general <- data.table::as.data.table(apply(general, MARGIN=2, trimws))
  
  # --- add the QC samples ---
  
  qc_samps = RSQLite::dbGetQuery(conn, "SELECT * FROM batchinfo WHERE sample LIKE '%QC%'")
  
  placeholder_date <- individual.data$sampling_date[[1]]
  
  qc_ind_data <- lapply(qc_samps$sample, function(qc) {
    data.table(label = c(1),
               card_id = qc,
               animal_internal_id = qc,
               sampling_date = placeholder_date,
               sex = "qc",
               group = "qc",
               farm = "QcLand")
  })
  
  qc_tab_setup = data.table(group = "qc",
                            stool_condition = "qc")
  qc_tab_ind = unique(rbindlist(qc_ind_data))
  
  # --- join to existing ---
  
  setup <- rbind(setup, qc_tab_setup, fill=TRUE)
  individual.data <- rbindlist(list(individual.data, qc_tab_ind), fill=TRUE)
  individual.data$label <- 1:nrow(individual.data)
  
  #pen.data <- data.table::as.data.table(apply(pen.data, MARGIN=2, trimws))
  #admin <- data.table::as.data.table(apply(admin, MARGIN=2, trimws))
  
  RSQLite::dbWriteTable(conn, "params", 
                        data.table::data.table(ppm=ppm), 
                        overwrite=T)
  
  # --- import to patient sql file ---
  #RSQLite::dbWriteTable(conn, "general", general, overwrite=TRUE) # insert into BUGGED FIX LATER
  RSQLite::dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "individual_data", individual.data, overwrite=TRUE) # insert into
  #RSQLite::dbWriteTable(conn, "pen_data", pen.data, overwrite=TRUE) # insert into
  #RSQLite::dbWriteTable(conn, "admin", admin, overwrite=TRUE) # insert into
  # --- disconnect ---
  RSQLite::dbDisconnect(conn)
}

db.build.custom <- function(db.name = "MyDb",
                            db.short = "mydb",
                            db.description = "Personal custom database.",
                            db.icon = "www/questionmark.png",
                            outfolder = getOptions(lcl$paths$opt.loc)$db_dir,
                            csv){
  
  db.base = data.table::fread(csv)
  
  columns =  c("compoundname",
               "description",
               "baseformula",
               "identifier",
               "charge",
               "structure")
  
  keep.columns <- intersect(columns,colnames(db.base))
  
  db.formatted <- db.base[, ..keep.columns]
  
  if(all(is.na(db.formatted$identifier))){
    db.formatted$identifier <- c(1:nrow(db.formatted))
  }
  
  # check the formulas
  checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                             db.formatted$baseformula))
  db.formatted$baseformula <- checked$new_formula
  keep <- checked[warning == FALSE, which = TRUE]
  db.formatted <- db.formatted[keep]
  
  # open db
  outfolder <- file.path(outfolder, "custom")
  if(!dir.exists(outfolder)) dir.create(outfolder)
  
  db <- file.path(outfolder, paste0(db.short, ".base.db"))
  if(file.exists(db)) file.remove(db)
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  RSQLite::dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text, structure text)")
  
  # create folder for db
  RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
  RSQLite::dbDisconnect(conn)
  
  # write metadata to file.. json?
  meta.dbpage =
    list(title = db.name,
         description = db.description,
         image_id = paste0(db.short, "_icon"))
  
  meta.img =
    list(name = paste0(db.short, "_icon"), path = db.icon, dimensions = c(200, 200))
  
  save(list = c("meta.img", "meta.dbpage"), file = file.path(outfolder, paste0(db.short, ".RData")))
}
