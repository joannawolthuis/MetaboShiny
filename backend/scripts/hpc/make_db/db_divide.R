dbname <- "pubchem"
dbdir <- "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db"
script_dir <-  "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/scripts/joanna"

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
    # ---
    write.csv(x = tab, file = file.path(csv_folder, csv_name))
  }
}

build.extended.db.send(dbname, dbdir, fetch.limit = 1000)

