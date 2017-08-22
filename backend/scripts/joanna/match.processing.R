#' @export
iden.code.binned <- function(outlist.path,
                             db.path,
                             isofilt=FALSE,
                             excl.adducts=c("PLACEHOLDER")){
  # --------------------------------- -------
  library(pbapply)
  library(parallel)
  library(data.table)
  library(gsubfn)
  library(DBI)
  library(RSQLite)
  # - connect -
  conn <- dbConnect(RSQLite::SQLite(), outlist.path)
  # --- skip if already exists ---
  result.table <- paste("matches", gsub("\\.db|\\.|full", "", basename(db.path)), sep="_")
  #if(dbExistsTable(conn, result.table)) return("Already exists!")
  # --- attach patient outlist and get mzmed pgrp values ---
  prep.query <- strwrap(fn$paste("ATTACH '$db.path' AS db"),width=10000, simplify=TRUE)
  dbExecute(conn, prep.query)
  # --------------------
  get.query <- strwrap(fn$paste("CREATE TEMPORARY TABLE patresults AS
                                SELECT DISTINCT base.compoundname, base.identifier, base.description, 
			                          cpd.baseformula, cpd.adduct, cpd.isoprevalence, mz.[mzmed.pgrp],
                                cpd.basemz, cpd.fullmz, cpd.basecharge, cpd.totalcharge
                                FROM mzvals mz INDEXED BY valindex
                                JOIN mzranges rng ON rng.ID = mz.ID  
                                LEFT JOIN db.extended cpd INDEXED BY cpdindex
                                ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
                                LEFT JOIN db.base base
                                ON base.baseformula = cpd.baseformula AND
                                base.charge = cpd.basecharge
                                ORDER BY mz.'mzmed.pgrp' ASC"),width=10000, simplify=TRUE)
  dbExecute(conn, get.query)
  print(dbGetQuery(conn, "select * from base limit 10"))
  dbExecute(conn, "DROP INDEX IF EXISTS pat_idx")
  dbExecute(conn, "CREATE INDEX pat_idx ON patresults(isoprevalence, adduct)")
  # --- isofilt? ---
  adduct.string <- gsub(" ", "", paste("'", excl.adducts, "'", collapse=","))
  filt.query <- strwrap(fn$paste("CREATE TEMPORARY TABLE isocount AS
                                  SELECT DISTINCT *
                                  FROM patresults INDEXED BY pat_idx
                                  WHERE (isoprevalence = 100
                                  AND adduct NOT IN ($adduct.string))"), width=10000, simplify=TRUE)
  dbExecute(conn, filt.query)
  dbExecute(conn, "DROP INDEX IF EXISTS iso_idx")
  dbExecute(conn, "CREATE INDEX iso_idx ON isocount(baseformula, adduct)")
  match.sql <- strwrap(fn$paste("SELECT DISTINCT pat.compoundname, pat.identifier, pat.description, 
                                                pat.baseformula, pat.adduct, 
                                                pat.[mzmed.pgrp], pat.isoprevalence, 
                                                pat.basecharge, pat.totalcharge 
                                      FROM patresults pat INDEXED BY pat_idx
                                      JOIN isocount iso INDEXED BY iso_idx
                                      ON pat.baseformula = iso.baseformula AND
                                      pat.adduct = iso.adduct"), width=10000, simplify=TRUE)
  dbExecute(conn, "DROP INDEX IF EXISTS pat_idx2")
  dbExecute(conn, "CREATE INDEX pat_idx2 ON patresults(baseformula)")
  nomatch.sql <- strwrap("SELECT DISTINCT compoundname, identifier, description,
                                          baseformula, adduct, 
                                          [mzmed.pgrp], isoprevalence, 
                                          basecharge, totalcharge
                         FROM patresults pat INDEXED BY pat_idx2
                         WHERE pat.baseformula ISNULL", width=10000, simplify=TRUE)
  # --- getto!! ---
  matches <- as.data.table(dbGetQuery(conn, match.sql))
  nomatches <- as.data.table(dbGetQuery(conn, nomatch.sql))
  # --- create result table
  dbExecute(conn, fn$paste("DROP TABLE IF EXISTS $result.table"))
  dbExecute(conn, fn$paste("CREATE TABLE $result.table(compoundname text, identifier text, description text, baseformula text, adduct text, [mzmed.pgrp] float, isoprevalence float, basecharge int, totalcharge int)"))
  # --- return ---
  results <- unique(rbind(matches, nomatches, fill=TRUE))
  # --- write to table ---
  dbWriteTable(conn, result.table, results, append=TRUE)
  # --- index ---
  dbExecute(conn, fn$paste("DROP INDEX IF EXISTS $result.table_idx1"))
  dbExecute(conn, fn$paste("CREATE INDEX $result.table_idx1 ON $result.table([mzmed.pgrp], baseformula, adduct)"))
  # --- add eventual extra indices here ---
  NULL
  # ---------------------------------------
  print("Done!")
  # --- disconnect ---
  dbDisconnect(conn)
  # --- AND return to R (can remove later) ---
  return(nrow(matches))
}

#' @export
find.formulas <- function(mzvals, cl=FALSE, ppm=3, charge=1, element.counts = list(c("C",0,50),c("H",0,50),
                                                                                   c("N",0,50),c("O",0,50),
                                                                                   c("S",0,50),c("Na", 0, 5),
                                                                                   c("Cl", 0, 5), c("P", 0,5))){
  require(rcdk)
  require(pbapply)
  require(data.table)
  # ------------------------------------------
  found.rows <- pblapply(mzvals,cl=cl, function(mz){
    window = mz * (ppm / 1e6)
    # --- generate molecular formulae ---
    found.mfs <- generate.formula(mz, window=0.3, element.counts, validation=TRUE, charge=charge)
    rows <- if(length(found.mfs) == 0) NA else(
      rows <- lapply(found.mfs, function(formula){
        # --- check for ppm range ---
        mz.found <- formula@mass
        within.ppm <- abs(mz - mz.found) < window
        if(within.ppm){
          data.table(origMZ = mz,
                     genMZ = mz.found,
                     BaseFormula = formula@string)
        } else(data.table(origMZ=mz, 
                          genMZ=NA, 
                          BaseFormula=NA))
      })
    )
    # --- return ---
    unique(rbindlist(rows[!is.na(rows)]))
  })
  rbindlist(found.rows[!is.na(found.rows)])
}
