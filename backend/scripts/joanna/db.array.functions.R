# example commands

# echo "/home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save < /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/db.make.hpc.R" | qsub -N make_pubchem -l h_vmem=20g -l h_rt=06:00:00 -pe threaded 40
# /home/cog/jwolthuis/R-3.4.0/bin/Rscript --vanilla < /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/db.make.hpc.R" | qsub -l h_vmem=20g -l h_rt=12:00:00 -pe threaded 80


# libraries

library(ggplot2)
library(DT)
library(DBI)
library(RSQLite)
library(gsubfn)
library(data.table)
library(parallel)
library(pbapply)
library(enviPat)
library(plotly)
library(jsonlite)
library(shinyFiles)
library(stringr)
library(ChemmineR)
library(curl)
library(httr)

# clone metaboshiny first...

dbDir <- "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db"

if(!dir.exists(dbDir)) dir.create(dbDir)

# init

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#' @export
build.extended.db.send <- function(dbname,
                                   out_folder, 
                                   fetch.limit=10){
  # ------------------------
  library(RSQLite)
  library(DBI)
  library(gsubfn)
  # ------------------------
  base.db <- file.path(out_folder, paste0(dbname, ".base.db"))
  # ------------------------
  base.conn <- dbConnect(RSQLite::SQLite(), base.db)
  # --------------------
  csv_folder <- file.path(out_folder, paste0(dbname, "_csv"))
  if(!dir.exists(csv_folder)) dir.create(csv_folder)
  
  get.query <- fn$paste("SELECT DISTINCT b.baseformula, b.charge FROM base b")
  total.formulae <- dbGetQuery(base.conn, fn$paste("SELECT Count(*)
                                                   FROM ($get.query)"))
  formula.count <- total.formulae[1,]
  results <- dbGetQuery(base.conn, get.query)
  # --- split and write batches? ----
  split_res <- split(results,(seq(nrow(results))-1) %/% fetch.limit)
  # write
  for(i in seq_along(split_res)){
    tab <- split_res[[i]]
    csv_name <- paste0(dbname, "_", i, ".csv")
    print(csv_name)
    print(tab)
    # ---
    write.csv(x = tab, file = file.path(csv_folder, csv_name))
  }
}

build.extended.db.send("hmdb",out_folder = dbDir,fetch.limit = 100)

# --- SEPERATE JOBS ---

build.extended.db.batch <- function(dbname,
                                    out_folder,
                                    batchnum,
                                    script_dir,
                                    adduct.table){
  # load necessary libraries (this should be a seperate script on its own)
  library(enviPat)
  library(data.table)
  data(isotopes, package = "enviPat")
  # ------------------------
  sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
      if(trace) cat(nm,":")           
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
  }
  sourceDir(script_dir)
  csv_out <- file.path(out_folder, paste0(dbname, "_csv_ext"))
  if(!dir.exists(csv_out)) dir.create(csv_out)
  # load in file with correct num
  csv_folder <- file.path(out_folder, paste0(dbname, "_csv"))
  csv_name <- paste0(dbname, "_", batchnum, ".csv")
  results <- read.csv(file.path(csv_folder,csv_name),header = T)
  # do the adduct thing
  checked.formulae <- as.data.table(check.chemform.joanna(isotopes, 
                                                          results$baseformula))
  # -----------------------------
  keep.rows <- checked.formulae[warning == FALSE, which=TRUE]
  # -----------------------------
  if(length(keep.rows) == 0) next
  backtrack <- data.table(baseformula = checked.formulae[keep.rows, new_formula],
                          basemz = checked.formulae[keep.rows, monoisotopic_mass],
                          charge = c(results$charge)[keep.rows])
  # -- go through each adduct --
  tab.list <- lapply(1:nrow(adduct.table), FUN=function(x){
    row <- adduct.table[x,]
    name <- row$Name
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
    extended.table <- data.table(baseformula = rep(backtrack.final$baseformula, repeat.times),
                                 fullformula = rep(backtrack.final$final, repeat.times),
                                 basemz = rep(backtrack.final$basemz, repeat.times),
                                 fullmz =isotable$fullmz,
                                 adduct = c(name),
                                 basecharge = rep(backtrack.final$charge, repeat.times),
                                 totalcharge = rep(backtrack.final$final.charge, repeat.times),
                                 isoprevalence = isotable$isoprevalence,
                                 foundinmode = c(mode)
    )
    # ------
    extended.table
  })
  print(tab.list[[1]])
  total.table <- rbindlist(tab.list[!is.na(tab.list)])
  write.table(total.table, file.path(csv_out, csv_name), row.names = F)
  # --- return ---
  # write result table to folder
}

for(i in 1:10){
  build.extended.db.batch("hmdb", dbDir, batchnum = i, script_dir = "./backend/scripts/joanna", wkz.adduct.confirmed)
}

build.extended.db.collect <- function(dbname,
                                      out_folder){
  library(RSQLite)
  library(DBI)
  library(gsubfn)
  # -----------------------
  csv_folder <- file.path(out_folder, paste0(dbname, "_csv_ext"))
  files <- list.files(csv_folder, pattern = ".csv$",full.names = T)
  # -----------------------
  full.db <- file.path(out_folder, paste0(dbname, ".full.db"))
  base.db <- file.path(out_folder, paste0(dbname, ".base.db"))
  
  if(file.exists(full.db)) file.remove(full.db)
  full.conn <- dbConnect(RSQLite::SQLite(), full.db)
  
  print("Attaching base...")
  dbExecute(full.conn, fn$paste("ATTACH '$base.db' AS tmp"))
  dbExecute(full.conn, fn$paste("CREATE TABLE IF NOT EXISTS done(baseformula text, basecharge text)"))
  dbExecute(full.conn, fn$paste("CREATE TABLE IF NOT EXISTS base AS SELECT * FROM tmp.base"))
  print("Indexing base...")
  dbExecute(full.conn, "CREATE INDEX IF NOT EXISTS b_idx1 on base(baseformula, charge)")
  
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
  # --- progress bar... ---
  sapply(files, FUN=function(f){
    print(f)
    result <- read.table(f, header = T)
    print(result)
    dbWriteTable(full.conn, "extended", result, append=T)
  })
  # --- indexy ---
  print("Indexing extended table...")
  dbExecute(full.conn, "CREATE INDEX IF NOT EXISTS b_idx1 on base(baseformula, charge)")
  dbExecute(full.conn, "create index e_idx1 on extended(baseformula, basecharge)")
  dbExecute(full.conn, "create index e_idx2 on extended(fullmz, foundinmode)")
  print("Disconnecting...")
  # --------------
  dbDisconnect(full.conn)
  print("Done! :-)")
}

build.extended.db.collect("hmdb",dbDir)
# === GET OPTIONS ===

wd <- "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny"

# --- laod adduct table for general use ---

load(file.path(wd, "backend/umcfiles/adducts/AdductTableWKZ.RData"))
sourceDir(file.path(wd, "backend/scripts/joanna"))
data(isotopes)

# get amount of cores

nslots <- Sys.getenv( "NSLOTS" )
print(nslots)

# --- wipe session log to save file size ---
session_cl <- makeCluster(nslots, type="FORK")

#build.base.db("pubchem", outfolder=dbDir, cl = session_cl)

build.extended.db("pubchem", 
                  continue=TRUE, 
                  outfolder=dbDir, 
                  adduct.table = wkz.adduct.confirmed, 
                  cl=session_cl, fetch.limit=1000)
