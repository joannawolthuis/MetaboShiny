#$ -S /home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

args = commandArgs(trailingOnly=TRUE)

i <- args[1]

load("/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/umcfiles/adducts/AdductTableWKZ.RData")

dbname <- "pubchem"
dbdir <- "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db"
script_dir <-  "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/joanna"

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

build.extended.db.batch(dbname, dbdir, batchnum = i, script_dir = script_dir, wkz.adduct.confirmed)
