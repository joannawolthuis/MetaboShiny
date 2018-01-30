groupingOnSQL <- function(f,
                  searchdb,
                  outdir,
                  mode,
                  overwrite=F,
                  ppm=2){
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(gsubfn)
  library(DBI)
  library(parallel)
  library(pbapply)
  # ------------------------
  ppm = as.numeric(ppm)
  load(f)
  load(paste(outdir, "repl.pattern.RData", sep="/"))
  load(paste(outdir, "breaks.fwhm.RData", sep="/"))
  if (mode=="negative"){
    repl.pattern=repl.pattern.neg
    groupNames=groupNames.neg
    
  } else {
    repl.pattern=repl.pattern.pos
    groupNames=groupNames.pos
  }
  options(digits=22)
  outlist.tot <- as.data.table(outlist.tot)
  print(colnames(outlist.tot))
  db.name <- gsub(f, pattern = "\\.RData",replacement = ".db")
  conn <- dbConnect(RSQLite::SQLite(), db.name)
  outlist.tot$foundinmode <- c(mode)
  # ------------------------
  dbWriteTable(conn,
               name = "peakvalues",
               outlist.tot,
               overwrite = T)
  # --- INDEX ---
  dbExecute(conn, "create index idx1 on peakvalues('mzmed.pkt', samplenr, 'height.pkt')")
  
  # === SEARCH ===
  query.zero <- fn$paste("ATTACH '$searchdb' AS db")
  print(query.zero)
  dbExecute(conn, query.zero)
  query.one <- fn$paste(strwrap(
    "SELECT mz.mzmed as groupmz, pk.'mzmed.pkt', pk.samplenr, pk.fq, rng.mzmin, rng.mzmax, sum(pk.'height.pkt') as intensity, pk.foundinmode
    FROM peakvalues pk
    LEFT JOIN db.mzranges rng
    ON pk.'mzmed.pkt' BETWEEN rng.mzmin AND rng.mzmax
    LEFT JOIN db.mzvals mz
    ON rng.ID = mz.ID AND mz.foundinmode = pk.foundinmode
    GROUP BY pk.'mzmed.pkt', pk.samplenr, groupmz
    ORDER BY groupmz ASC",
    width=10000,
    simplify=TRUE))
  
  res <- as.data.table(dbGetQuery(conn, query.one))

  unident.peaks <- unique(res[is.na(groupmz),])
  peakgroups <- unique(res[!is.na(groupmz),])
  outrows <- list()
  restrows <- list()
  used <- c()
  
  i = 1
  for(mz in unique(res$groupmz)){
    #print(paste("Used values:", paste(used, collapse=",")))
    # FIND PEAKS IN THIS GROUP
    allmembers <- peakgroups[groupmz == mz,]
    members <- allmembers[!(mzmed.pkt %in% used),]
    usedmembers <- allmembers[mzmed.pkt %in% used,]
    idkused <- nrow(allmembers)
    unused <- nrow(members)
    alreadyused <- nrow(usedmembers)
    if(!is.null(idkused)){
      allmembercount = idkused
    } else{ allmembercount = 0 }
    if(!is.null(unused)){
      membercount = unused
    } else{ membercount = 0 }
    if(!is.null(alreadyused)){
      usedcount = nrow(usedmembers)
    } else{ usedcount = 0 }
    print(paste("Members in this group:", allmembercount))
    print(paste("Already used, so omitting:", usedcount))
    #members <- members[!(members %in% used)]
    if(membercount > 0){
      if(!is.null(alreadyused)){
        restrows[[i]] <- as.data.table(usedmembers)
      }
      print(paste("---", mz, "---"))
      used <- unique(c(used, as.numeric(members[,mzmed.pkt])))
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
      outrows[[i]] <- as.data.table(outpgrlist)
      i = i + 1
    }else{print("NO UNUSED PEAKS LEFT")}
    # ========================
  }
  outpgrlist <- rbindlist(outrows)
  resttab <- rbindlist(restrows)
  unident.peaks <- rbind(unident.peaks, resttab)
  # --- return ---
  dir.create(paste(outdir, "grouping_db", sep="/"), showWarnings = FALSE)
  save(outpgrlist, file=paste(paste(outdir, "grouping_db", sep="/"), paste(mode, "RData", sep="."), sep="/"))
  dir.create(paste(outdir, "specpks_all_rest", sep="/"), showWarnings = FALSE)
  save(unident.peaks, file=paste(paste(outdir, "specpks_all_rest", sep="/"), paste(mode, "RData", sep="."), sep="/"))
}
