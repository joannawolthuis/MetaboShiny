library(RSQLite)
library(DBI)
library(data.table)

dbFolder <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db"
dbs <- list.files(dbFolder, 
                  pattern = "\\.full\\.db$",
                  full.names = T)

dbFrags <- lapply(dbs, FUN=function(db){
  conn <- dbConnect(RSQLite::SQLite(), db)
  print(db)
  # -----------------
  contents <- dbGetQuery(conn, 
                         "SELECT DISTINCT fullmz as mz, baseformula as formula, foundinmode as mode FROM extended")
  # -----------------
  dbDisconnect(conn)
  contents$source <- c(gsub(basename(db), 
                            pattern = "\\.full\\.db", 
                            replacement = ""))
  return(contents)
})

bigDB <- unique(rbindlist(dbFrags))

bigDB_pathways <- bigDB[which(bigDB$source %in% c("smpdb", "kegg", "wikipathways")),]
bigDB_pathways <- bigDB_pathways[order(mz),]

references <- bigDB[mode == "positive"]
save(references, file = "positive.RData")
head(references)

references <- bigDB[mode == "negative"]
save(references, file = "negative.RData")