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

# === GET OPTIONS ===

wd <- "/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny"

# --- laod adduct table for general use ---

load(file.path(wd, "backend/umcfiles/adducts/AdductTableWKZ.RData"))
sourceDir(file.path(wd, "backend/scripts/joanna"))
data(isotopes)

# get amount of cores
nslots <- Sys.getenv( "NSLOTS" )
print( nslots )

session_cl <<- if(is.na(session_cl)) makeCluster(nslots, type="FORK")

build.base.db("pubchem", outfolder=dbDir, cl = session_cl)
build.extended.db("pubchem", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=100)
