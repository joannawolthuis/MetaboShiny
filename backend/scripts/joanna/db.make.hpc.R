# clone metaboshiny first...

dbDir <- "/hpc/cog_bioinf/ridder/users/jwolthuis/shinyDB"

if(!dir.exists(dbDir)) dir.create(dbDir)

# init

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# === GET OPTIONS ===

wd <- "/Users/jwolthuis/Google Drive/MetaboShiny"

# --- laod adduct table for general use ---

load(file.path(wd, "backend/umcfiles/adducts/AdductTableWKZ.RData"))
sourceDir(file.path(wd, "backend/scripts/joanna"))
data(isotopes)

# get amount of cores
session_cl <<- if(is.na(session_cl)) makeCluster(4, type="FORK")




build.base.db("pubchem", outfolder=dbDir, cl = session_cl)
build.extended.db("pubchem", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
