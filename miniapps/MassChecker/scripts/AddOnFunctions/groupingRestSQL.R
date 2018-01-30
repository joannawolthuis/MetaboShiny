groupingRestSQL <- function(outdir, fileIn, cl=NULL, scanmode, ppm=2){
  # fileIn="./results/specpks_all_rest/negative_outlist_i_min_1.1.RData"
  # scanmode="negative"
  # outdir="./results"
  # ppm=2
  
  library(data.table)
  library(gsubfn)
  library(DBI)
  library(RSQLite)
  library(pbapply)
  
  options(digits=16)
  load(fileIn)
  print(fileIn)
  unident.peaks <- unident.peaks[,-c(1,5,6)]
  colnames(unident.peaks)[1] <- "mzmed"
  print(head(unident.peaks))
  
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  
  
  if (scanmode=="negative"){
    # repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
  } else {
    # repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  # load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  
  outpgrlist = NULL
  
  if (scanmode=="negative"){
    # repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
  } else {
    # repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  
  # Then group on highest peaks
  range = ppm*1e-06
  startcol=7
  
  # === SQL-IFY ===
  db.name <- gsub(fileIn, pattern = "\\.RData$", replacement = ".db")
  
  print(db.name)
  
  if(file.exists(db.name)) file.remove(db.name)
  
  conn <- dbConnect(RSQLite::SQLite(), db.name)
  
  # --- write intensities to table and index ---
  
  sql.make.meta <- strwrap("CREATE TABLE mzintensities(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           samplenr text,
                           fq decimal(30,13),
                           intensity decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  
  dbExecute(conn, sql.make.meta)   
  
  dbWriteTable(conn, "mzintensities", unique(unident.peaks), append=T) # insert into
  
  dbExecute(conn, "CREATE INDEX intindex ON mzintensities(mzmed)")
  
  mzvals <- data.table(mzmed = unique(unident.peaks$mzmed),
                       foundinmode = c(scanmode))
  
  mzranges <- data.table(mzmin = pbsapply(as.numeric(mzvals$mzmed),cl = cl,
                                          FUN=function(mz, ppm){
                                            mz - mz * (ppm / 1E6)}, ppm=ppm),
                         mzmax = pbsapply(as.numeric(mzvals$mzmed), cl = cl,
                                          FUN=function(mz, ppm){
                                            mz + mz * (ppm / 1E6)}, ppm=ppm))
  
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  
  dbExecute(conn, sql.make.meta)   
  dbExecute(conn, "create index mzfind on mzvals(mzmed);")
  
  # --- write vals to table ---
  dbWriteTable(conn, "mzvals", unique(mzvals), append=T) # insert into
  
  # --- make range table (choose if R*tree or not) ---
  
  sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges USING rtree(
                            ID INTEGER PRIMARY KEY AUTOINCREMENT,
                            mzmin decimal(30,13),
                            mzmax decimal(30,13));"
                            , width=10000, simplify=TRUE)
  
  dbExecute(conn, sql.make.rtree)
  
  # --- write ranges to table ---
  dbWriteTable(conn, "mzranges", unique(mzranges), append=T) # insert into
  
  # ------------------
  print("a")
  
  
  dbExecute(conn, "DROP TABLE IF EXISTS results")
  
  dbExecute(conn, "CREATE TABLE results AS
            SELECT DISTINCT int.mzmed as peakmz, mz.mzmed as groupmz, rng.mzmin as mzmin, rng.mzmax as mzmax, int.intensity as intensity, int.fq as fq, int.samplenr as samplenr
            FROM mzintensities int
            JOIN mzranges rng
            ON int.mzmed BETWEEN rng.mzmin AND rng.mzmax
            JOIN mzvals mz
            ON rng.ID = mz.ID")
  
  print("b")
  
  dbExecute(conn, "CREATE INDEX res1 ON results(groupmz, peakmz, intensity, mzmin, mzmax, fq, samplenr)")
  
  dbExecute(conn, "DROP TABLE IF EXISTS results_grouped")
  
  dbExecute(conn, "
            CREATE TABLE results_grouped AS
            SELECT DISTINCT groupmz, peakmz, samplenr, intensity, mzmin, mzmax, fq FROM results res WHERE groupmz IN(
            SELECT groupmz from results WHERE
            (groupmz,peakmz,intensity) IN 
            ( Select distinct groupmz, groupmz, MAX(intensity) from results
            Group by groupmz
            ))
            ORDER BY (SELECT MAX(intensity) FROM results WHERE groupmz = res.groupmz) DESC, intensity desc")
  
  print("c")
  
  dbExecute(conn, "
            CREATE INDEX res2 ON results_grouped(peakmz, samplenr)")
  
  print(paste("Used peaks removed:", dbGetQuery(conn, "
                                                SELECT count(rowid) FROM results_grouped
                                                WHERE rowid NOT IN
                                                (
                                                SELECT min(rowid)
                                                FROM results_grouped
                                                GROUP BY
                                                peakmz,
                                                samplenr
                                                );
                                                ")))
  
  dbExecute(conn, "
            DELETE FROM results_grouped
            WHERE rowid NOT IN
            (
            SELECT  min(rowid)
            FROM    results_grouped
            GROUP BY
            peakmz,
            samplenr
            )")
  
  print("d")
  
  
  # --- cleanup ---
  

  
  res <- as.data.table(dbGetQuery(conn, "SELECT * FROM results_grouped"))
  
  dbDisconnect(conn)
  
  # === IT WORKS UP UNTIL HERE!!! ===
  
  return(res)
  
  cl=parallel::makeCluster(3,"FORK")
  
  # --- get intensities ---
  
  outrows <- pblapply(unique(res$groupmz), cl=cl,FUN=function(mz){
    print(paste("---", mz, "---"))
    # FIND PEAKS IN THIS GROUP
    members <- res[groupmz == mz,]
    # --- get intensities ---
    mzmin = as.numeric(min(members$mzmin))
    mzmax = as.numeric(max(members$mzmax))
    fq.worst.pgrp = as.numeric(max(members$fq))
    fq.best.pgrp = as.numeric(min(members$fq))
    ints.allsamps = rep(0, length(groupNames))
    names(ints.allsamps) = groupNames # same order as sample list!!!
    # # Check for each sample if multiple peaks exists, if so take the sum!
    labels=unique(members$samplenr)
    nrsamples <- length(labels)
    ints.allsamps[labels] = as.vector(unlist(lapply(labels, 
                                                    function(x){
                                                      sum(members[samplenr==x, 'intensity'])
                                                    })))
    # --- make dt ---
    outpgrlist = rbind(cbind('mzmed.pgrp' = mz, 
                             "fq.best"=fq.best.pgrp, 
                             "fq.worst"=fq.worst.pgrp, 
                             nrsamples, 
                             'mzmin.pgrp' = mzmin, 
                             'mzmax.pgrp' = mzmax,          
                             t(as.matrix(ints.allsamps)
                             )))
    # -----------------------
    as.data.table(outpgrlist)
  })
  
  outpgrlist <- rbindlist(outrows)
  colnames(outpgrlist)[1:6] = c("mzmed.pgrp", "fq.best", "fq.worst", "nrsamples", "mzmin.pgrp", "mzmax.pgrp")
  
  # --- return ---
  dir.create(paste(outdir, "grouping_rest", sep="/"), showWarnings = FALSE)
  
  save(outpgrlist, file=paste(outdir, "grouping_rest", paste(scanmode,".RData", sep=""), sep="/"))
}
