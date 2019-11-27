#' @export
build.pat.db <- function(db.name,
                         metapath,
                         pospath,
                         negpath,
                         overwrite=FALSE,
                         rtree=TRUE,
                         make.full = TRUE,
                         ppm=2,
                         inshiny=F,
                         wipe.regex = ".*_(?>POS|NEG)_[0+]*"){

  ppm = as.numeric(ppm)

  # ==== TESTING ====
  
  #metapath = "~/Documents/Documents/code/metshi_paper/MTBLS28_20181107_105647/s_mtbls28_v2.txt"
  #pospath = "~/Documents/Documents/code/metshi_paper/MTBLS28_20181107_105647/m_mtbls28_POS_v2_maf.tsv"
  #negpath = "~/Documents/Documents/code/metshi_paper/MTBLS28_20181107_105647/m_mtbls28_NEG_v2_maf.tsv"
  
  # =================
  
  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)
  
  metadata <- data.table::fread(metapath)
  colnames(metadata) <- tolower(colnames(metadata))
  
  colnames(metadata) <- tolower(gsub(x=colnames(metadata), pattern = "\\.$|\\.\\.$|\\]", replacement = ""))
  colnames(metadata) <- gsub(x=colnames(metadata), pattern = "\\[|\\.|\\.\\.| ", replacement = "_")

  keep.cols = sapply(colnames(metadata), function(x) length(unique(metadata[,..x][[1]])) > 1)  
  
  metadata.filt = metadata[, ..keep.cols]
  
  colnames(metadata.filt) <- gsub(colnames(metadata.filt), pattern = "characteristics_|factor_value_", replacement="")
  setnames(metadata.filt, old = c("sample_name", "source_name"), new = c("sample", "individual"), skip_absent = T)
  if(!("individual" %in% colnames(metadata.filt))) metadata.filt$individual <- metadata.filt$sample
  
  RSQLite::dbWriteTable(conn, "individual_data", metadata.filt, overwrite=TRUE) # insert into
  
  # =================
  
  poslist <- data.table::fread(pospath,header = T)
  neglist <- data.table::fread(negpath,header = T) 
  colnames(poslist) <- gsub(colnames(poslist), pattern = wipe.regex, replacement = "", perl=T)
  colnames(neglist) <- gsub(colnames(neglist), pattern = wipe.regex, replacement = "", perl=T)
  
  unique.samples = unique(metadata.filt$sample)
  
  keep.samples.pos <- intersect(colnames(poslist), unique.samples)
  keep.samples.neg <- intersect(colnames(neglist), unique.samples)
  
  missing.samples = setdiff(unique.samples, unique(c(keep.samples.pos, keep.samples.neg)))
  if(length(missing.samples) > 0 ){
    shiny::showNotification(paste0("Missing samples: ", paste0(missing.samples, collapse = ","), ". Please check your peak table < > metadata for naming mismatches!"))
  }
  
  setnames(poslist, "mass_to_charge", "mzmed", skip_absent = T)
  setnames(neglist, "mass_to_charge", "mzmed", skip_absent = T)
  
  keepcols.pos <- c("mzmed", keep.samples.pos)
  keepcols.neg <- c("mzmed", keep.samples.neg)
  
  poslist <- poslist[, ..keepcols.pos]
  neglist <- neglist[, ..keepcols.neg]
  
  # replace commas with dots
  poslist <- data.table::as.data.table(sapply(poslist, gsub, pattern = ",", replacement= "."))
  neglist <- data.table::as.data.table(sapply(neglist, gsub, pattern = ",", replacement= "."))
  
 if(inshiny) setProgress(.20)
  
  mzvals <- data.table::data.table(mzmed = c(as.numeric(poslist$mzmed), as.numeric(neglist$mzmed)),
                                   foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  
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

