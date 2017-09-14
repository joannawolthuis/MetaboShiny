build.extended.db.collect <- function(dbname,
                                      out_folder){
  library(RSQLite)
  library(DBI)
  library(gsubfn)
  library(sqldf)
  library(data.table)
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
  sql.make.meta <- strwrap("DROP TABLE IF EXISTS extended")
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
  for(f in files){
    print(f)
    res <- fread(f, sep = " ", quote = '"')
    dbWriteTable(full.conn, "extended", res, append=T)
  }
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

build.extended.db.collect("pubchem", "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db")


