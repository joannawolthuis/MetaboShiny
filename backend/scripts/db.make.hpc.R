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
if(file.exists("/home/cog/jwolthuis/shinyClusterLog.txt")) file.remove("/home/cog/jwolthuis/shinyClusterLog.txt")

session_cl <- makeCluster(nslots, type="FORK")

#build.base.db("pubchem", outfolder=dbDir, cl = session_cl)

build.extended.db("pubchem", outfolder=dbDir, adduct.table = wkz.adduct.confirmed, cl=session_cl, fetch.limit=1000)
