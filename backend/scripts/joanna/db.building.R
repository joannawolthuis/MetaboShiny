#' @export
build.base.db <- function(dbname=NA, 
                          outfolder="/Users/jwolthuis/Google Drive/Metabolomics/Databases", 
                          cl=FALSE){
  # --- check if user chose something ---
  if(is.na(dbname)) return("~ Please choose one of the following options: HMDB, ChEBI, PubChem! d(>..w.<)b ~")
  # -------------------------------------
  library(data.table)
  library(pbapply)
  library(RSQLite)
  library(DBI)
  library(sqldf)
  # --- make dir if not exists ---
  if(!dir.exists(outfolder)) dir.create(outfolder)
  # --- create .db file ---
  db <- file.path(outfolder, paste0(dbname, ".base.db"))
  if(file.exists(db)) file.remove(db)
  conn <- dbConnect(RSQLite::SQLite(), db)
  dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text)")
  # --- create base ---
  function.of.choice <- switch(tolower(dbname),
                               internal = function(dbname){ # BOTH NOISE AND NORMAL
                                 # --- uses csv package ---
                                 int.loc <- file.path(wd, "backend","umcfiles", "internal")
                                 # --- non-noise ---
                                 internal.base.db <- read.csv(file.path(int.loc, 
                                                                     "TheoreticalMZ_NegPos_noNoise.txt"),
                                                           sep="\t")
                                 db.formatted <- data.table(compoundname = internal.base.db[,1],
                                                            description = c("Internal"),
                                                            baseformula = internal.base.db[,2], 
                                                            identifier=c("Internal"),
                                                            charge=c(0))
                                 # fix some unicode stuff
                                 db.formatted <- db.formatted[baseformula != ""]
                                 db.formatted$baseformula <-  gsub(x=db.formatted$baseformula, 
                                                                   pattern="( )|(\xa0)|(\xe1)|( .*$)", 
                                                                   replacement="")
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # --- write ---
                                 dbWriteTable(conn, "base", db.formatted, append=TRUE)},
                               noise = function(dbname){
                                 int.loc <- file.path(wd, "backend","umcfiles", "internal")
                                 # --- noise ---
                                 noise.base.db <- read.csv(file.path(int.loc, 
                                                                  "TheoreticalMZ_NegPos_yesNoise.txt"), 
                                                        sep="\t")
                                 db.formatted <- data.table(compoundname = noise.base.db[,1],
                                                            description = str_match(noise.base.db[,1], "(?<=\\().+?(?=\\))"),
                                                            baseformula = noise.base.db[,2], 
                                                            identifier=c("Noise"),
                                                            charge=c(0))
                                 dbWriteTable(conn, "base", db.formatted, append=TRUE)
                                 ### do extended table and write to additional db (exception...)
                                db.full <- file.path(outfolder, paste0(dbname, ".full.db"))
                                if(file.exists(db.full)) file.remove(db.full)
                                full.conn <- dbConnect(RSQLite::SQLite(), db.full)
                                # --- create base table here too ---
                                dbExecute(full.conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text)")
                                # --- create extended table ---
                                sql.make.meta <- strwrap("create table extended(
                                                         baseformula text,
                                                         fullformula text, 
                                                         basemz decimal(30,13), 
                                                         fullmz decimal(30,13), 
                                                         adduct text,
                                                         basecharge int,
                                                         totalcharge int,
                                                         isoprevalence float,
                                                         foundinmode text)", width=10000, simplify=TRUE)
                                dbExecute(full.conn, sql.make.meta)
                                # -- reformat noise table ---
                                db.formatted.list <- pblapply(1:nrow(noise.base.db), FUN=function(x){
                                  row <- c(noise.base.db[x,])
                                  print(row[1])
                                  # --- find adduct ---
                                  split.col <- strsplit(unlist(row[1]), split = "\\[|\\]")
                                  compoundname <- split.col[[1]][1]
                                  adduct <- split.col[[1]][2]
                                  # --------------
                                  mz.neg <- as.numeric(row[3])
                                  mz.pos <- as.numeric(row[4])
                                  
                                  # --------------
                                  result <- data.table(
                                      baseformula = as.character(row[2]),
                                      fullformula = as.character(row[2]), 
                                      basemz = c(mz.pos, mz.neg), 
                                      fullmz = c(mz.pos, mz.neg), 
                                      adduct = if(is.na(adduct)) c("M+H", "M-H") else c(adduct),
                                      basecharge = c(0),
                                      totalcharge = c(1, -1),
                                      isoprevalence = c(100),
                                      foundinmode = c("positive", "negative")
                                    )
                                  result
                                })
                                db.formatted.full <- rbindlist(db.formatted.list[!is.na(db.formatted.list)])
                                db.formatted.full <- db.formatted.full[basemz > 0] # remove not-found ions
                                # --- write to db ---
                                db.formatted$compoundname <- as.character(strsplit(db.formatted$compoundname, 
                                                                      split=", \\[.*\\][\\+-]"))
                                dbWriteTable(full.conn, 
                                             "base", 
                                             db.formatted, 
                                             append=TRUE)
                                dbWriteTable(full.conn, 
                                             "extended", 
                                             db.formatted.full, 
                                             append=TRUE)
                                # --- index ---
                                dbExecute(full.conn, "create index b_idx1 on base(baseformula, charge)")
                                dbExecute(full.conn, "create index e_idx1 on extended(baseformula, basecharge)")
                                dbExecute(full.conn, "create index e_idx2 on extended(fullmz, foundinmode)")# -------------
                                dbDisconnect(full.conn)
                                },
                               hmdb = function(dbname){
                                 library(XML)
                                 library(plyr)
                                 print("Downloading XML database...")
                                 file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
                                 # ----
                                 base.loc <- file.path(dbDir, "hmdb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc)
                                 zip.file <- file.path(base.loc, "HMDB.zip")
                                 download.file(file.url, zip.file)
                                 unzip(zip.file, exdir = base.loc)
                                 # --- go through xml ---
                                 print("Parsing XML...")
                                 data <- xmlParse(file.path(base.loc,"hmdb_metabolites.xml"), useInternalNodes = T)
                                 # --- xpath magic! :-) --- [NOTE: leaves out anything that doens't have formal charge available, another option would be to default to zero]
                                 compoundnames <- getNodeSet(data, 
                                                             "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:name", 
                                                             c(pf = "http://www.hmdb.ca"))
                                 formulae <- getNodeSet(data, 
                                                        "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:chemical_formula", 
                                                        c(pf = "http://www.hmdb.ca"))
                                 identifiers <- getNodeSet(data, 
                                                           "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:accession", 
                                                           c(pf = "http://www.hmdb.ca"))
                                 charges <- getNodeSet(data, 
                                                       "/*/pf:metabolite/pf:predicted_properties/pf:property[pf:kind='formal_charge']/pf:value", 
                                                       c(pf = "http://www.hmdb.ca"))
                                 description <- getNodeSet(data, 
                                                           "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:description", 
                                                           c(pf = "http://www.hmdb.ca"))
                                 # --- make nice big table ---
                                 db.formatted <- data.table(
                                   compoundname = sapply(compoundnames, xmlValue),
                                   description = sapply(description, xmlValue),
                                   baseformula = sapply(formulae, xmlValue),
                                   identifier = sapply(identifiers, xmlValue),
                                   charge = sapply(charges, xmlValue)
                                 )
                                 # --- check formulae ---
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                          db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # ----------------------
                                 dbWriteTable(conn, "base", db.formatted, append=TRUE)
                               },
                               chebi = function(dbname){
                                 library(minval)
                                 print("you chose chebi")
                                 db.full <- as.data.table(download.chebi.joanna(release = "latest", woAssociations = FALSE))
                                 db.formatted <- unique(db.full[, list(compoundname = ChEBI, 
                                                                       description = DEFINITION,
                                                                       baseformula = FORMULA,
                                                                       identifier = ID, 
                                                                       charge = gsub(CHARGE,pattern = "$\\+\\d", replacement = "")
                                                                       )
                                                                ]
                                                        )
                                 # --- check formulae ---
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # ----------------------
                                 dbWriteTable(conn, "base", db.formatted, append=TRUE)
                               },
                               pubchem = function(dbname, ...){
                                 library(stringr)
                                 library(ChemmineR)
                                 library(parallel)
                                 library(curl)
                                 library(httr)
                                 # --- create working space ---
                                 baseLoc <- file.path(dbDir, "pubchem_source")
                                 sdf.loc <- file.path(baseLoc, "sdf")
                                 csv.loc <- file.path(baseLoc, "csv")
                                 # --- check the user's determination ---
                                 # continue <- readline("You chose PubChem, this will take a while and requires at least 30gb of space.Continue? (yes/no): ")
                                 # if(continue == "no" | continue == "n") return(NA)
                                 # --- download the 1 million sdf files ---
                                 if(!dir.exists(sdf.loc)) dir.create(sdf.loc)
                                 if(!dir.exists(csv.loc)) dir.create(csv.loc)
                                 print("Downloading SDF files...")
                                 folder.url = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
                                 ftp.handle = new_handle(dirlistonly=TRUE)
                                 ftp.conn = curl(folder.url, "r", ftp.handle)
                                 file.table = read.table(ftp.conn, stringsAsFactors=TRUE, fill=TRUE)
                                 close(ftp.conn)
                                 file.urls <- paste0(folder.url, file.table[,1])
                                 # ---------------------------------------
                                 cb <- function(req){
                                   cat("done:", req$url, ": HTTP:", req$status, "\n")
                                   file.handle <- file(file.path(sdf.loc,basename(req$url)), open = "wb")
                                   counter <<- counter + 1
                                   print(paste(counter, max_counter, sep="//"))
                                   writeBin(req$content, file.handle)
                                   close(file.handle)
                                   # --- download ---
                                 }
                                 counter  = 0
                                 max_counter = length(file.urls)
                                 # --- start downloady ---
                                 pool <- new_pool(host_con = length(cl))
                                 sapply(file.urls, FUN=function(url){
                                   #if(file.exists(file.path(sdf.loc, basename(url)))) return(NA)
                                   curl_fetch_multi(url, done = cb, pool = pool)
                                 })
                                 # lotsa files, need tiem.
                                 
                                 out <- multi_run(pool=pool)
                                 # ------------------------------
                                 print("Converting SDF files to tab delimited matrices...")
                                 sdf.files <- list.files(path = sdf.loc, pattern = "\\.sdf\\.gz$")
                                 # -----------------------
                                 pbsapply(cl=cl, sdf.files, FUN=function(sdf.file){
                                   input <- file.path(sdf.loc, sdf.file)
                                   output <- file.path(csv.loc, gsub("\\.sdf.gz$", "\\.csv", sdf.file))
                                   sdfStream.joanna(input=input,
                                             output=output,
                                             fct = function(sdfset, test){
                                               valid <- validSDF(sdfset)
                                               sdfset <- sdfset[valid]
                                               print(head(datablock(sdfset)))
                                               blockmatrix <- datablock2ma(datablock(sdfset)) # Converts data block to matrix
                                               # --------------
                                               db.formatted <- data.table(
                                                 compoundname = blockmatrix[, if("PUBCHEM_IUPAC_TRADITIONAL_NAME" %not in% colnames(blockmatrix)) "PUBCHEM_MOLECULAR_FORMULA" else("PUBCHEM_IUPAC_TRADITIONAL_NAME")],
                                                 description = c("PubChem"),
                                                 baseformula = gsub(x = blockmatrix[, "PUBCHEM_MOLECULAR_FORMULA"], 
                                                                    pattern="[\\+\\-]\\d*$", 
                                                                    replacement=""),
                                                 identifier = as.numeric(blockmatrix[, "PUBCHEM_COMPOUND_CID"]),
                                                 charge = blockmatrix[,"PUBCHEM_TOTAL_CHARGE"]
                                                )
                                               # --- check formulae ---
                                               checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                              db.formatted$baseformula))
                                               db.formatted$baseformula <- checked$new_formula
                                               keep <- checked[warning == FALSE, which = TRUE]
                                               db.keep <- db.formatted[keep]
                                               db.keep},
                                             append = FALSE,
                                             silent = FALSE,
                                             Nlines = 100000 )
                                 })
                                 # --- assemble ---
                                 csv.files <- list.files(path = csv.loc, pattern = "\\.csv$", full.names = TRUE)
                                 print("Assembling and putting in db file...")
                                 pbsapply(cl=cl, csv.files, FUN=function(file){
                                   print(file)
                                   first.row <- read.csv(file, nrows=3, header=TRUE, sep="\t")
                                   if("description" %not in% colnames(first.row)){print("NOPE!"); file.remove(file); return(NULL)}
                                   read.csv.sql(file, sep="\t", sql = c(fn$paste("insert into base select compoundname, description, baseformula, identifier, charge from file")), dbname = db)
                                 })
                                 dbDisconnect(conn)
                                 print("Done!")
                                 # --- cleanup --- downloading takes forever, would not recommend unless compressing... ---
                                 # remove.residuals <- readline("Done! Remove downloaded SDF / XLS files? (yes/no): ")
                                 # if(remove.residuals == "yes" | remove.residuals == "y") unlink(sdf.loc, recursive=TRUE) else return("Alrighty, files kept.")
                               })
  # --- execute ;) ---
  function.of.choice(dbname)
}

#' @export
build.extended.db <- function(dbname, 
                              outfolder, 
                              adduct.table, 
                              charge.mode, 
                              continue=FALSE, 
                              reset.base=FALSE, 
                              cl=FALSE, 
                              fetch.limit=-1, 
                              test.mode="OFF",
                              cpd.limit=-1){
  # --- GET BASE DATA ---
  library(enviPat)
  library(pbapply)
  library(data.table)
  library(RSQLite)
  library(DBI)
  library(sqldf)
  # ------------------------
  data(isotopes)
  base.db <- file.path(outfolder, paste0(dbname, ".base.db"))
  full.db <- file.path(outfolder, paste0(dbname, ".full.db"))
  # ------------------------
  if(file.exists(full.db)) file.remove(full.db)
  full.conn <- dbConnect(RSQLite::SQLite(), full.db)
  base.conn <- dbConnect(RSQLite::SQLite(), base.db)
  # ------------------------
  sql.make.meta <- strwrap("CREATE TABLE IF NOT EXISTS extended(
                           baseformula text,
                           fullformula text, 
                           basemz decimal(30,13), 
                           fullmz decimal(30,13), 
                           adduct text,
                           basecharge int,
                           totalcharge int,
                           isoprevalence float,
                           foundinmode text)", width=10000, simplify=TRUE)
  dbExecute(full.conn, sql.make.meta)
  # ------------------------
  limit.query <- if(cpd.limit == -1) "" else fn$paste("LIMIT $cpd.limit")
  # ------------------------
  total.formulae <- dbGetQuery(base.conn, fn$paste("SELECT Count(*)
                               FROM (SELECT DISTINCT
                               baseformula, charge
                               FROM base $limit.query)"))
  print(total.formulae)
  results <- dbSendQuery(base.conn, fn$paste("SELECT DISTINCT baseformula, charge FROM base $limit.query"))
  formula.count <- total.formulae[1,]
  # --- start pb ---
  pb <- startpb(0, formula.count)
  # --- waow, my first while in R ---
  while(!dbHasCompleted(results)){
    # --- fetch part of results ---
    partial.results <- as.data.table(dbFetch(results, fetch.limit))
    if(length(partial.results$baseformula) == 0) next
    # -----------------------
    checked.formulae <- as.data.table(check.chemform.joanna(isotopes, 
                                                            partial.results$baseformula))
    keep.rows <- checked.formulae[warning == FALSE & monoisotopic_mass %between% c(60, 600), which=TRUE]
    # -----------------------------
    if(length(keep.rows) == 0) next
    backtrack <- data.table(baseformula = checked.formulae[keep.rows, new_formula],
                            basemz = checked.formulae[keep.rows, monoisotopic_mass],
                            charge = c(partial.results$charge)[keep.rows])
    # -- go through each adduct --
    do.calc <- function(x){
      row <- adduct.table[x,]
      name <- row$Name
      if(test.mode == "ON") print(name)
      adduct <- row$Formula_add
      deduct <- row$Formula_ded
      # --- fix booleans ---
      adduct <- if(adduct == "FALSE") FALSE else adduct
      deduct <- if(deduct == "FALSE") FALSE else deduct
      # --------------------
      mode <- row$Ion_mode
      multiplier <- as.numeric(row$Multi)
      # --- multiplication ---
      formulae.mult <- multiform.joanna(backtrack$baseformula, multiplier)
      backtrack$multiform <- formulae.mult
      backtrack$adducted <- backtrack$multiform
      # --- is adduction necessary? ---
      if(adduct != FALSE){
        formulae.add <- mergeform(formula1 = backtrack$baseformula, 
                                  formula2 = adduct)
        backtrack$adducted <- formulae.add
      }
      # --- placeholder ---
      backtrack$final <- backtrack$adducted
      # --- is deduction necessary? ---
      if(deduct != FALSE){
        formulae <- backtrack$final
        can.deduct <- which(!check.ded.joanna(formulas = backtrack$adducted,
                                              deduct = deduct))
        if(length(can.deduct) == 0) return(NA)
        deductibles <- formulae[can.deduct]
        formulae.ded <- subform(deductibles, deduct)
        backtrack$final[can.deduct] <- formulae.ded
        backtrack <- backtrack[can.deduct]
      }
      if(nrow(backtrack) == 0) return(NA)
      backtrack$final.charge <- c(as.numeric(backtrack$charge)) + c(as.numeric(row$Charge))
      backtrack <- backtrack[final.charge != 0]
      if(nrow(backtrack) == 0) return(NA)
      # --- get isotopes ---
      isotopes <- isopattern(
        isotopes,
        backtrack$final,
        threshold = 0.1,
        plotit = FALSE,
        charge = backtrack$final.charge,
        algo = 2,
        verbose = FALSE
      )
      # -------------------------
      isolist <- lapply(isotopes, function(isotable){
        if(test.mode == "ON") print(isotable[[1]])
        if(isotable[[1]] == "error") return(NA)
        iso.dt <- data.table(isotable, fill=TRUE)
        result <- iso.dt[,1:2]
        names(result) <- c("fullmz", "isoprevalence")
        # --- return ---
        result
      })
      isolist.nonas <- isolist[!is.na(isolist)]
      isotable <- rbindlist(isolist.nonas)
      keep.isos <- names(isolist.nonas)
      # --- remove 'backtrack' rows that couldn't be calculated ---
      backtrack.final <- backtrack[final %in% keep.isos]
      # --------------------
      repeat.times <- c(unlist(lapply(isolist.nonas, FUN=function(list) nrow(list)), use.names = FALSE))
      meta.table <- data.table(baseformula = rep(backtrack.final$baseformula, repeat.times),
                               fullformula = rep(backtrack.final$final, repeat.times),
                               basemz = rep(backtrack.final$basemz, repeat.times),
                               fullmz =isotable$fullmz,
                               adduct = c(name),
                               basecharge = rep(backtrack.final$charge, repeat.times),
                               totalcharge = rep(backtrack.final$final.charge, repeat.times),
                               isoprevalence = isotable$isoprevalence,
                               foundinmode = c(mode)
      )
      # --- return ---
      meta.table
    }
    tab.list <- switch(test.mode,
                       ON = pblapply(cl=FALSE, 1:nrow(adduct.table), FUN=function(x) do.calc(x)),
                       OFF = parLapply(cl=cl, 1:nrow(adduct.table), fun=function(x) do.calc(x))
    )
    # --- progress bar... ---
    setpb(pb, dbGetRowCount(results))
    total.table <- rbindlist(tab.list[!is.na(tab.list)])
    dbWriteTable(full.conn, "extended", total.table, append=TRUE)
  }
  dbClearResult(results)
  # --- add base db to the new one ---
  dbExecute(full.conn, fn$paste("ATTACH '$base.db' AS tmp"))
  dbExecute(full.conn, fn$paste("CREATE TABLE base AS SELECT * FROM tmp.base"))
  # --- indexy ---
  dbExecute(full.conn, "create index b_idx1 on base(baseformula, charge)")
  dbExecute(full.conn, "create index e_idx1 on extended(baseformula, basecharge)")
  dbExecute(full.conn, "create index e_idx2 on extended(fullmz, foundinmode)")
  # --- vacuum --- !! needs altered settings
  dbGetQuery(full.conn, "VACUUM")
  # --------------
  dbDisconnect(base.conn)
  dbDisconnect(full.conn)
  on.exit(closepb(pb))
  return("Done! :-)")
}

#' @export
build.pat.db <- function(db.name, 
                      poslist, 
                      neglist, 
                      overwrite=FALSE,
                      rtree=TRUE,
                      rmv.cols=c("mzmin.pgrp", 
                                 "mzmax.pgrp",
                                 "fq.best", 
                                 "fq.worst",
                                 "nrsamples",
                                 "avg.int")){
  library(reshape2)
  library(DBI)
  library(RSQLite)
  # ------------------------
  mzvals <- data.table(mzmed.pgrp = c(poslist$mzmed.pgrp, neglist$mzmed.pgrp),
                       foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  mzranges <- data.table(mzmin = c(poslist$mzmin.pgrp, neglist$mzmin.pgrp),
                         mzmax = c(poslist$mzmax.pgrp, neglist$mzmax.pgrp))
  mzintensities <- melt(as.data.table(rbind(poslist, neglist))[,(rmv.cols) := NULL], 
                        id="mzmed.pgrp", 
                        variable="filename",
                        value="intensity")
  print("here")
  filenames <- gsub(x=as.character(mzintensities$filename), "", pattern="-\\d*$", perl=T)
  replicates <- gsub(x=as.character(mzintensities$filename), "", pattern=".*-(?=\\d*$)", perl=T)
  print("here")
  # ------------------------
  mzintensities$filename <- filenames
  mzintensities$replicate <- replicates
  # ------------------------
  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)
  # --- reconnect / remake ---
  conn <- dbConnect(RSQLite::SQLite(), db.name)
  print(unique(mzintensities$filename))
  print(unique(mzintensities$replicate))
  # ------------------------
  sql.make.int <- strwrap("CREATE TABLE mzintensities(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           'mzmed.pgrp' decimal(30,13),
                           filename text,
                           intensity float,
                            replicate int)", width=10000, simplify=TRUE)
  dbExecute(conn, sql.make.int)
  # --- write intensities to table and index ---
  dbWriteTable(conn, "mzintensities", mzintensities, append=TRUE) # insert into
  dbExecute(conn, "CREATE INDEX intindex ON mzintensities(filename,'mzmed.pgrp',intensity)")
  # ------------------------
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           'mzmed.pgrp' decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  dbExecute(conn, sql.make.meta)
  dbExecute(conn, "create index mzfind on mzvals([mzmed.pgrp], foundinmode);")
  # --- write vals to table ---
  dbWriteTable(conn, "mzvals", mzvals, append=TRUE) # insert into
  # --- make range table (choose if R*tree or not) ---
  sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges USING rtree(
                            ID INTEGER PRIMARY KEY AUTOINCREMENT,
                            mzmin decimal(30,13),
                            mzmax decimal(30,13));"
                            , width=10000, simplify=TRUE)
  sql.make.normal <- strwrap("CREATE TABLE mzranges(
                             ID INTEGER PRIMARY KEY AUTOINCREMENT,
                             mzmin decimal(30,13),
                             mzmax decimal(30,13));", width=10000, simplify=TRUE)
  dbExecute(conn, if(rtree) sql.make.rtree else sql.make.normal)
  # --- write ranges to table ---
  dbWriteTable(conn, "mzranges", mzranges, append=TRUE) # insert into
  # --- cleanup ---
  dbGetQuery(conn, "VACUUM")
  # ----------------
  dbDisconnect(conn)
  print("Made!")}

#' @export
load.excel <- function(path.to.xlsx, 
                       path.to.patdb = patdb,
                       tabs.to.read = c("General",
                                        "Setup",
                                        "Individual Data",
                                        "Pen Data",
                                        "Admin")){
  # --- should put in patdb sql file ---
  library(xlsx)
  library(DBI)
  library(RSQLite)
  # --- connect to sqlite db ---
  conn <- dbConnect(RSQLite::SQLite(), path.to.patdb)
  # -------------------------------
  tab.store <- pblapply(tabs.to.read, FUN=function(tab.name){
    tab <- as.data.table(read.xlsx(path.to.xlsx, sheetName = tab.name))
    print(tab)
    # --- reformat colnames ---
    split.cols <- strsplit(colnames(tab), split = "\\.\\.")
    reformat.cols <- lapply(split.cols, FUN=function(col.split){
        # --- remove trailing dot ---
        col.split <- gsub(col.split, pattern="\\.$", replacement="", perl=T)
        # ---------------------------
        col.type <- gsub(tolower(col.split[1]), pattern="\\.", replacement="_", perl=T)
        col.unit <- gsub(tolower(col.split[2]), pattern="\\.", replacement="/", perl=T)
        # --- return ---
        data.table(table=c(tab.name),
                   column = as.character(col.type), 
                   unit = as.character(col.unit))})
    unit.table <- rbindlist(reformat.cols)
    print(head(tab))
    colnames(tab) <- unit.table$column
    # --- return both ---
    result <- list(units = unit.table, data = tab)
  })
  # --- create unit table ---
  unit.store <- list()
  data.store <- list()
  for(x in 1:length(tab.store)){
    curr <- tab.store[[x]]
    unit.store[[x]] <- curr$units
    data.store[[x]] <- curr$data
    }
  units <- rbindlist(unit.store)
  # --- convert to data table --- ## make this nicer loooking in the future
  general <- data.store[[1]]
  setup <- data.store[[2]]
  individual.data <- data.store[[3]]
  individual.data$sampling_date <- as.factor(as.Date(individual.data$sampling_date, format = "%d-%m-%y"))
  #individual.data$sampling_date <- as.numeric(as.factor(as.Date(individual.data$sampling_date, format = "%d-%m-%y")))
  pen.data <- data.store[[4]]
  admin <- data.store[[5]]
  print("here")
  print(individual.data)
  # --- import to patient sql file ---
  #dbWriteTable(conn, "general", general, overwrite=TRUE) # insert into BUGGED FIX LATER
  dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  dbWriteTable(conn, "units", units, overwrite=TRUE) # insert into
  dbWriteTable(conn, "individual_data", individual.data, overwrite=TRUE, 
               field.types=list(label="integer",
                                card_id="text",
                                animal_internal_id="text",
                                sampling_date="integer")) # insert into
  dbWriteTable(conn, "pen_data", pen.data, overwrite=TRUE) # insert into
  dbWriteTable(conn, "admin", admin, overwrite=TRUE) # insert into
  # --- disconnect ---
  dbDisconnect(conn)
  return(colnames(setup))
}
