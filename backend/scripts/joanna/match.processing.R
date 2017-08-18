#' @export
get.conservative.adducts <- function(matches,
                                     source.db,
                                     min.frac=0){
  conn <- dbConnect(RSQLite::SQLite(), source.db)
  # --- subset internal database ---
  internal.standards <- dbGetQuery(conn, 'SELECT DISTINCT compoundname, baseformula FROM base WHERE compoundname LIKE "%(IS)%"')
  # --- find them in matches ---
  standard.tabs <- pblapply(1:nrow(internal.standards), FUN=function(x){
    row <- internal.standards[x,]
    standard.matches <- matches[BaseFormula == row$BaseFormula,
                                nomatch = 0L]
    return(standard.matches)
  })
  counts <- adduct.counter(rbindlist(standard.tabs))
  print(counts)
  to.keep.adducts <- counts[Occurence >= (min.frac * max(Occurence)), Adduct]
  # --- return ---
  to.keep.adducts
}

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
  # --- attach patient outlist and get mzmed pgrp values ---
  prep.query <- strwrap(fn$paste("ATTACH '$db.path' AS db"),width=10000, simplify=TRUE)
  print(prep.query)
  dbExecute(conn, prep.query)
  print("here")
  # --------------------
  get.query <- strwrap(fn$paste("CREATE TEMPORARY TABLE patresults AS
                                SELECT DISTINCT base.compoundname, base.identifier, 
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
  print('here')
  print(dbGetQuery(conn, "select * from base limit 10"))
  dbExecute(conn, "DROP INDEX IF EXISTS pat_idx")
  dbExecute(conn, "CREATE INDEX pat_idx ON patresults(isoprevalence, adduct, baseformula, compoundname, identifier)")
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
  match.sql <- strwrap(fn$paste("SELECT DISTINCT pat.compoundname, pat.identifier, 
                                                pat.baseformula, pat.adduct, 
                                                pat.[mzmed.pgrp], pat.isoprevalence, 
                                                pat.basecharge, pat.totalcharge 
                                      FROM patresults pat INDEXED BY pat_idx
                                      JOIN isocount iso INDEXED BY iso_idx
                                      ON pat.baseformula = iso.baseformula AND
                                      pat.adduct = iso.adduct"), width=10000, simplify=TRUE)
  dbExecute(conn, "DROP INDEX IF EXISTS pat_idx2")
  dbExecute(conn, "CREATE INDEX pat_idx2 ON patresults(baseformula)")
  nomatch.sql <- strwrap("SELECT DISTINCT compoundname, identifier, 
                                          baseformula, adduct, 
                                          [mzmed.pgrp], isoprevalence, 
                                          basecharge, totalcharge
                         FROM patresults pat INDEXED BY pat_idx2
                         WHERE pat.baseformula ISNULL", width=10000, simplify=TRUE)
  # --- getto!! ---
  matches <- as.data.table(dbGetQuery(conn, match.sql))
  nomatches <- as.data.table(dbGetQuery(conn, nomatch.sql))
  # --- create result table
  result.table <- paste("matches", gsub("\\.db|\\.|full", "", basename(db.path)), sep="_")
  dbExecute(conn, fn$paste("DROP TABLE IF EXISTS $result.table"))
  dbExecute(conn, fn$paste("CREATE TABLE $result.table(compoundname text, identifier text, baseformula text, adduct text, [mzmed.pgrp] float, isoprevalence float, basecharge int, totalcharge int)"))
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
  return(results)
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
